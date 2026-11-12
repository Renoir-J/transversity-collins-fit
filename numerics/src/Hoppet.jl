# Fixed linear maps through an unchanged retained backend. No AD through settings.
const HOPPET_NATIVE_LOCK=ReentrantLock()
const HOPPET_CONTRACT=normpath(joinpath(@__DIR__,"../../config/hoppet-contract.json"))

mutable struct HoppetState
    library::Ptr{Nothing}
    context::Int64
    functions::Dict{Symbol,Ptr{Nothing}}
end

struct HoppetPlan
    state::HoppetState
    grid::Vector{Float64}
    Q2::Vector{Float64}
    Q02::Float64
    tau::Vector{Float64}
    quantity::Symbol
    nf::Int
    coupling_model::Symbol
    lambda_GeV::Union{Nothing,Float64}
    layout::Symbol
    native_nodes::Int
    coefficient_payload_bytes::Int64
    signature::UInt
end

function verified_hoppet_build(path)
    receipt=JSON.parsefile(path)
    expected=JSON.parsefile(HOPPET_CONTRACT)["upstream"]
    receipt["status"]=="BUILT_NOT_VALIDATED" && receipt["scope"]=="reconstructed_theory" &&
        receipt["production_ready"]===false && receipt["ABI_version"]==1 &&
        receipt["upstream"]==expected && receipt["source_clean"]===true || error("invalid HOPPET build scope/version")
    for kind in ("source_hashes","output_hashes")
        !isempty(receipt[kind]) || error("missing HOPPET build provenance")
        for (file,digest) in receipt[kind]
            isfile(file) && bytes2hex(sha256(read(file)))==digest || error("stale HOPPET build: $file")
        end
    end
    for (file,target) in receipt["source_links"]
        islink(file) && readlink(file)==target || error("stale HOPPET source link: $file")
    end
    receipt
end

function HoppetPlan(buildfile,Q2,Q02,coupling;nf=5,quantity=:transversity,
                    ymax=nothing,xmin=nothing,dy=.025,order=-6,layout=:nested,eps=1e-10,du=.0125)
    quantity in (:transversity,:KPSY_Hhat3_homogeneous) || error("unsupported HOPPET quantity/convention")
    nf isa Integer && !(nf isa Bool) && 3<=nf<=6 || error("fixed integer nf from three to six required")
    layout in (:single,:nested) || error("choose single or nested HOPPET grid")
    order isa Integer && !(order isa Bool) && -10<=order<=-2 || error("negative order -2 through -10 required")
    qs=[fixed_scalar(q,"Q2") for q in Q2]; q0=fixed_scalar(Q02,"Q02")
    1<=length(qs)<=4096 && q0>0 && all(>(0),qs) || error("finite positive fixed scales required")
    ymax!==nothing && xmin!==nothing && error("specify xmin or ymax, not both")
    y,lower = if xmin!==nothing
        low=fixed_scalar(xmin,"xmin")
        0<low<1 || error("xmin must lie strictly between zero and one")
        (-log(low),low)
    elseif ymax!==nothing
        high=fixed_scalar(ymax,"ymax")
        (high,exp(-high))
    else
        (-log(1e-6),1e-6)
    end
    settings=[fixed_scalar(v,n) for (v,n) in zip((y,dy,eps,du),("ymax","dy","eps","du"))]
    y,d,e,u=settings
    0<y<=30 && d>0 && 0<e<=1e-4 && 0<u<=.2 && abs(order)+1<=round(y/d)<=20000 || error("unsupported grid/solver settings")
    layout==:nested && !(y>=2.5 && d<=.2) && error("nested grid requires ymax>=2.5 and dy<=0.2")
    times,model,lambda = if coupling isa AsymptoticTwoLoopCoupling
        coupling.nf==nf || error("coupling/nf mismatch")
        (evolution_time.(Ref(coupling),qs,q0),:asymptotic_two_loop,coupling.lambda_GeV)
    else
        a=fixed_scalar(coupling,"alpha0")
        ([-2log(alpha_lo(q,q0,a;nf=nf)/a)/(11-2nf/3) for q in qs],:one_loop,nothing)
    end
    all(t -> isfinite(t) && abs(t)<=1,times) || error("unsupported auxiliary evolution-time interval")
    receipt=verified_hoppet_build(buildfile)
    lock(HOPPET_NATIVE_LOCK) do
        library=Libdl.dlopen(receipt["library"])
        state=HoppetState(library,0,Dict(name=>Libdl.dlsym(library,name) for name in
            (:cf_hoppet_open,:cf_hoppet_close,:cf_hoppet_grid,:cf_hoppet_forward,
             :cf_hoppet_adjoint,:cf_hoppet_generator,:cf_hoppet_evaluate,:cf_hoppet_info)))
        try
            ccall(Libdl.dlsym(library,:cf_hoppet_abi),Cint,())==1 || error("wrong HOPPET ABI")
            n=Ref{Cint}(0); context=Ref{Int64}(0); status=Ref{Cint}(0)
            ccall(state.functions[:cf_hoppet_open],Cvoid,
                (Cdouble,Cdouble,Cdouble,Cint,Cint,Cdouble,Cdouble,Cint,Cint,Ptr{Cdouble},Ref{Cint},Ref{Int64},Ref{Cint}),
                y,lower,d,order,layout==:nested ? 1 : 0,e,u,nf,length(times),times,n,context,status)
            status[]==0 || error("HOPPET initialization rejected (status $(status[])); use one active context per isolated job")
            state.context=context[]
            x=zeros(n[]); native=Ref{Cint}(0); payload=Ref{Int64}(0)
            ccall(state.functions[:cf_hoppet_grid],Cvoid,(Int64,Cint,Ptr{Cdouble},Ref{Cint}),state.context,n[],x,status)
            status[]==0 && all(diff(x).>0) && last(x)==1 || error("invalid HOPPET grid")
            ccall(state.functions[:cf_hoppet_info],Cvoid,(Int64,Cint,Ref{Cint},Ref{Int64},Ref{Cint}),state.context,n[],native,payload,status)
            status[]==0 || error("missing HOPPET resource metadata")
            signature=hash((x,qs,times))
            finalizer(close,state)
            HoppetPlan(state,x,qs,q0,times,quantity,nf,model,lambda,layout,native[],payload[],signature)
        catch
            close(state)
            rethrow()
        end
    end
end

function Base.close(state::HoppetState)
    lock(HOPPET_NATIVE_LOCK) do
        if state.context!=0
            status=Ref{Cint}(0)
            ccall(state.functions[:cf_hoppet_close],Cvoid,(Int64,Ref{Cint}),state.context,status)
            status[]==0 || error("stale HOPPET context during close")
            state.context=0
        end
        if state.library!=C_NULL
            Libdl.dlclose(state.library)
            state.library=C_NULL
        end
    end
    nothing
end
Base.close(plan::HoppetPlan)=close(plan.state)
Base.isopen(plan::HoppetPlan)=plan.state.context!=0

function hoppet_valid(plan,initial;endpoint=true)
    isopen(plan) || error("closed HOPPET plan")
    hash((plan.grid,plan.Q2,plan.tau))==plan.signature || error("mutated HOPPET settings/grid")
    size(initial,1)==length(plan.grid) && size(initial,2)>0 && all(isfinite,initial) || error("invalid HOPPET initial shape/values")
    endpoint && !all(iszero,selectdim(initial,1,size(initial,1))) && error("h(1)=0 is required")
end

# Linearity transports every tangent, including nested ForwardDiff values. The
# components are concatenated along the batch dimension before crossing the ABI.
function hoppet_dual_batch(operation,initial::AbstractArray{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    pieces=(map(ForwardDiff.value,initial),ntuple(k -> map(v -> ForwardDiff.partials(v)[k],initial),N)...)
    combined=operation(cat(pieces...;dims=2))
    nb=size(initial,2)
    values=selectdim(combined,2,1:nb)
    partials=ntuple(k -> selectdim(combined,2,k*nb+1:(k+1)*nb),N)
    map((v,ds...) -> ForwardDiff.Dual{T}(v,ForwardDiff.Partials(ds)),values,partials...)
end

function apply_operator(plan::HoppetPlan,initial::AbstractMatrix{Float64})
    lock(HOPPET_NATIVE_LOCK) do
        hoppet_valid(plan,initial)
        result=Array{Float64}(undef,size(initial)...,length(plan.tau)); status=Ref{Cint}(0)
        GC.@preserve plan initial result ccall(plan.state.functions[:cf_hoppet_forward],Cvoid,
            (Int64,Cint,Cint,Cint,Ptr{Cdouble},Ptr{Cdouble},Ref{Cint}),
            plan.state.context,size(initial,1),size(initial,2),length(plan.tau),Matrix(initial),result,status)
        status[]==0 || error("HOPPET evolution rejected (status $(status[]))")
        result
    end
end
apply_operator(plan::HoppetPlan,initial::AbstractMatrix{<:ForwardDiff.Dual})=hoppet_dual_batch(x -> apply_operator(plan,x),initial)
apply_operator(plan::HoppetPlan,initial::AbstractVector)=dropdims(apply_operator(plan,reshape(initial,:,1));dims=2)
apply_operator(plan::HoppetPlan,initial::AbstractMatrix)=error("HOPPET supports Float64 and nested Float64 ForwardDiff values")

function operator_adjoint(plan::HoppetPlan,g::AbstractArray{Float64,3})
    lock(HOPPET_NATIVE_LOCK) do
        hoppet_valid(plan,g;endpoint=false)
        size(g,3)==length(plan.tau) || error("HOPPET cotangent scale count mismatch")
        result=zeros(size(g,1),size(g,2)); status=Ref{Cint}(0)
        ccall(plan.state.functions[:cf_hoppet_adjoint],Cvoid,
            (Int64,Cint,Cint,Cint,Ptr{Cdouble},Ptr{Cdouble},Ref{Cint}),
            plan.state.context,size(g,1),size(g,2),size(g,3),Array(g),result,status)
        status[]==0 || error("HOPPET adjoint rejected (status $(status[]))")
        result
    end
end
operator_adjoint(plan::HoppetPlan,g::AbstractArray{<:ForwardDiff.Dual,3})=hoppet_dual_batch(x -> operator_adjoint(plan,x),g)
operator_adjoint(plan::HoppetPlan,g::AbstractMatrix)=vec(operator_adjoint(plan,reshape(g,size(g,1),1,size(g,2))))

function hoppet_generator(plan::HoppetPlan,initial::AbstractMatrix{Float64})
    lock(HOPPET_NATIVE_LOCK) do
        hoppet_valid(plan,initial)
        result=zeros(size(initial)); status=Ref{Cint}(0)
        ccall(plan.state.functions[:cf_hoppet_generator],Cvoid,
            (Int64,Cint,Cint,Ptr{Cdouble},Ptr{Cdouble},Ref{Cint}),
            plan.state.context,size(initial,1),size(initial,2),Matrix(initial),result,status)
        status[]==0 || error("HOPPET generator rejected (status $(status[]))")
        result
    end
end
hoppet_generator(plan::HoppetPlan,initial::AbstractMatrix{<:ForwardDiff.Dual})=hoppet_dual_batch(x -> hoppet_generator(plan,x),initial)
hoppet_generator(plan::HoppetPlan,initial::AbstractVector)=vec(hoppet_generator(plan,reshape(initial,:,1)))

function hoppet_evaluate(plan::HoppetPlan,initial::AbstractMatrix{Float64},queries)
    q=[fixed_scalar(v,"HOPPET query fraction") for v in queries]
    !isempty(q) && all(v -> first(plan.grid)<=v<=1,q) || error("HOPPET query outside support")
    lock(HOPPET_NATIVE_LOCK) do
        hoppet_valid(plan,initial)
        result=zeros(length(q),size(initial,2)); status=Ref{Cint}(0)
        ccall(plan.state.functions[:cf_hoppet_evaluate],Cvoid,
            (Int64,Cint,Cint,Cint,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ref{Cint}),
            plan.state.context,size(initial,1),size(initial,2),length(q),Matrix(initial),q,result,status)
        status[]==0 || error("HOPPET interpolation rejected (status $(status[]))")
        result
    end
end
hoppet_evaluate(plan::HoppetPlan,initial::AbstractMatrix{<:ForwardDiff.Dual},q)=hoppet_dual_batch(x -> hoppet_evaluate(plan,x,q),initial)
hoppet_evaluate(plan::HoppetPlan,initial::AbstractVector,q)=vec(hoppet_evaluate(plan,reshape(initial,:,1),q))
function hoppet_evaluate(plan::HoppetPlan,initial::AbstractArray{T,3},q) where T
    evaluated=hoppet_evaluate(plan,reshape(initial,size(initial,1),:),q)
    reshape(evaluated,length(q),size(initial,2),size(initial,3))
end
