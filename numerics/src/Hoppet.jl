# Fixed linear maps through an unchanged retained backend. No AD through settings.
const HOPPET_NATIVE_LOCK=ReentrantLock()
const HOPPET_CONTRACT=normpath(joinpath(@__DIR__,"../../config/hoppet-contract.json"))

# An internal token for one serial point-plan preparation. Verify the complete
# build once before and once after assembling its rows, rather than per query.
struct HoppetBuildVerification
    receipt::Dict{String,Any}
end

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
    construction::Symbol
    native_nodes::Int
    coefficient_payload_bytes::Int64
    signature::String
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
                    ymax=nothing,xmin=nothing,dy=.025,order=-6,layout=:nested,eps=1e-10,du=.0125,
                    construction=:batched_probes)
    quantity in (:transversity,:KPSY_Hhat3_homogeneous) || error("unsupported HOPPET quantity/convention")
    nf isa Integer && !(nf isa Bool) && 3<=nf<=6 || error("fixed integer nf from three to six required")
    layout in (:single,:nested) || error("choose single or nested HOPPET grid")
    construction in (:independent,:batched_probes) || error("unsupported retained construction method")
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
    receipt=buildfile isa HoppetBuildVerification ? buildfile.receipt : verified_hoppet_build(buildfile)
    lock(HOPPET_NATIVE_LOCK) do
        library=Libdl.dlopen(receipt["library"])
        state=HoppetState(library,0,Dict(name=>Libdl.dlsym(library,name) for name in
            (:cf_hoppet_open,:cf_hoppet_open_batched,:cf_hoppet_close,:cf_hoppet_grid,:cf_hoppet_forward,
             :cf_hoppet_adjoint,:cf_hoppet_generator,:cf_hoppet_evaluate,:cf_hoppet_info,:cf_hoppet_rows,
             :cf_hoppet_collins,:cf_hoppet_collins_rows)))
        try
            ccall(Libdl.dlsym(library,:cf_hoppet_abi),Cint,())==1 || error("wrong HOPPET ABI")
            n=Ref{Cint}(0); context=Ref{Int64}(0); status=Ref{Cint}(0)
            operation=construction==:batched_probes ? :cf_hoppet_open_batched : :cf_hoppet_open
            ccall(state.functions[operation],Cvoid,
                (Cdouble,Cdouble,Cdouble,Cint,Cint,Cdouble,Cdouble,Cint,Cint,Ptr{Cdouble},Ref{Cint},Ref{Int64},Ref{Cint}),
                y,lower,d,order,layout==:nested ? 1 : 0,e,u,nf,length(times),times,n,context,status)
            status[]==0 || error("HOPPET initialization rejected (status $(status[])); use one active context per isolated job")
            state.context=context[]
            x=zeros(n[]); native=Ref{Cint}(0); payload=Ref{Int64}(0)
            ccall(state.functions[:cf_hoppet_grid],Cvoid,(Int64,Cint,Ptr{Cdouble},Ref{Cint}),state.context,n[],x,status)
            status[]==0 && all(diff(x).>0) && last(x)==1 || error("invalid HOPPET grid")
            ccall(state.functions[:cf_hoppet_info],Cvoid,(Int64,Cint,Ref{Cint},Ref{Int64},Ref{Cint}),state.context,n[],native,payload,status)
            status[]==0 || error("missing HOPPET resource metadata")
            signature=numerical_snapshot((x,qs,times))
            finalizer(close,state)
            HoppetPlan(state,x,qs,q0,times,quantity,nf,model,lambda,layout,construction,native[],payload[],signature)
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
    numerical_snapshot((plan.grid,plan.Q2,plan.tau))==plan.signature || error("mutated HOPPET settings/grid")
    size(initial,1)==length(plan.grid) && size(initial,2)>0 && all(isfinite,initial) || error("invalid HOPPET initial shape/values")
    endpoint && !all(iszero,selectdim(initial,1,size(initial,1))) && error("h(1)=0 is required")
end

# Linearity transports every tangent, including nested ForwardDiff values.
# Components share a batch through either the Float64 ABI or detached Julia rows.
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

"Fixed selected-node rows with O(nodes*scales) output and O(native nodes) workspace."
function hoppet_rows(plan::HoppetPlan,node;collins=false)
    node isa Integer && !(node isa Bool) && 1<=node<=length(plan.grid) || error("fixed valid grid node required")
    collins isa Bool || error("fixed Collins row selection required")
    collins && plan.layout!=:single && error("regular matching rows require a single query-bound grid")
    lock(HOPPET_NATIVE_LOCK) do
        isopen(plan) || error("closed HOPPET plan")
        numerical_snapshot((plan.grid,plan.Q2,plan.tau))==plan.signature || error("mutated HOPPET settings/grid")
        rows=zeros(length(plan.tau),length(plan.grid));status=Ref{Cint}(0)
        operation=collins ? :cf_hoppet_collins_rows : :cf_hoppet_rows
        ccall(plan.state.functions[operation],Cvoid,
            (Int64,Cint,Cint,Cint,Ptr{Cdouble},Ref{Cint}),
            plan.state.context,length(plan.grid),length(plan.tau),node,rows,status)
        status[]==0 || error("HOPPET row export rejected (status $(status[]))")
        rows
    end
end

"Full retained regular matching, used to check detached row construction."
function hoppet_collins_regular(plan::HoppetPlan,initial::AbstractMatrix{Float64})
    plan.layout==:single || error("regular matching requires a single query-bound grid")
    lock(HOPPET_NATIVE_LOCK) do
        hoppet_valid(plan,initial)
        result=zeros(size(initial));status=Ref{Cint}(0)
        ccall(plan.state.functions[:cf_hoppet_collins],Cvoid,
            (Int64,Cint,Cint,Ptr{Cdouble},Ptr{Cdouble},Ref{Cint}),
            plan.state.context,size(initial,1),size(initial,2),Matrix(initial),result,status)
        status[]==0 || error("HOPPET regular matching rejected (status $(status[]))")
        result
    end
end
function hoppet_collins_regular(plan::HoppetPlan,initial::AbstractArray{Float64,3})
    reshape(hoppet_collins_regular(plan,reshape(initial,size(initial,1),:)),size(initial))
end

struct HoppetPointBlock
    indices::Vector{Int}
    rows::FixedRowMatrix
    HoppetPointBlock(indices::Vector{Int},rows::AbstractMatrix{Float64})=new(copy(indices),fixed_rows(rows))
end

"Detached query-bound rows from retained HOPPET; no native context in the fit."
struct HoppetPointPlan
    grid::Vector{Float64}
    queries::Vector{Float64}
    Q2::Vector{Float64}
    Q02::Float64
    tau::Vector{Float64}
    quantity::Symbol
    nf::Int
    intervals::Int
    construction::Symbol
    blocks::Vector{HoppetPointBlock}
    coefficient_payload_bytes::Int64
    row_native_max_scaled::Float64
    signature::String
end

"Distinct process-matched rows; physical alpha_s, never the auxiliary coupling."
struct HoppetMatchingPlan
    grid::Vector{Float64}
    queries::Vector{Float64}
    Q2::Vector{Float64}
    Q02::Float64
    tau::Vector{Float64}
    quantity::Symbol
    nf::Int
    intervals::Int
    construction::Symbol
    blocks::Vector{HoppetPointBlock}
    coefficient_payload_bytes::Int64
    row_native_max_scaled::Float64
    process::Symbol
    alpha::Vector{Float64}
    signature::String
end

function prepare_hoppet_points(buildfile,queries,Q2,Q02,coupling;intervals=1024,nf=5,
                               quantity=:transversity,order=-6,eps=1e-10,du=.0125,matching=nothing,
                               construction=:batched_probes)
    intervals isa Integer && !(intervals isa Bool) && 11<=intervals<=20000 || error("fixed interval count 11..20000 required")
    qs=[fixed_scalar(x,"query fraction") for x in queries]
    !isempty(qs) && all(x -> 0<x<=1,qs) && any(<(1),qs) || error("interior queries in (0,1), with optional endpoint, required")
    contract=JSON.parsefile(normpath(joinpath(@__DIR__,"../../config/hoppet-point-contract.json")))
    verification=HoppetBuildVerification(verified_hoppet_build(buildfile))
    prepared=Dict{Float64,Any}(); coordinates=Float64[1.]; scales=times=nothing
    row_errors=Float64[]
    for query in unique(qs)
        query==1 && continue
        plan=HoppetPlan(verification,Q2,Q02,coupling;nf=nf,quantity=quantity,xmin=query,
            dy=-log(query)/intervals,order=order,layout=:single,eps=eps,du=du,construction=construction)
        try
            n,nt=length(plan.grid),length(plan.tau)
            rows=hoppet_rows(plan,1)
            x=plan.grid
            test=hcat(x.*(1 .-x).^3,-2x.^2 .* (1 .-x).^2)
            evolved=apply_operator(plan,test)
            reference=permutedims(evolved[1,:,:])
            if matching!==nothing
                a=matching.alpha./pi
                prefix=1 .+ a.*matching_delta(matching.process,matching.channel)
                if matching.channel==:collins
                    regular_rows=hoppet_rows(plan,1;collins=true)
                    regular=permutedims(hoppet_collins_regular(plan,evolved)[1,:,:])
                    rows=prefix.*rows .+ a.*regular_rows
                    reference=prefix.*reference .+ a.*regular
                else
                    rows=prefix.*rows
                    reference=prefix.*reference
                end
            end
            actual=rows*test
            row_error=maximum(abs.(actual-reference)./max.(abs.(reference),.01maximum(abs,reference),1e-15))
            isfinite(row_error) && row_error<=contract["validation_policy"]["row_native_max_scaled"] || error("query-row/native mismatch")
            push!(row_errors,row_error)
            prepared[query]=(grid=copy(x),rows=fixed_rows(rows))
            append!(coordinates,x)
            scales,times=copy(plan.Q2),copy(plan.tau)
        finally
            close(plan)
        end
    end
    grid=sort!(unique(coordinates)); blocks=HoppetPointBlock[]
    for query in qs
        if query==1
            push!(blocks,HoppetPointBlock(Int[],zeros(length(times),0)))
        else
            native=prepared[query]
            indices=searchsortedfirst.(Ref(grid),native.grid)
            push!(blocks,HoppetPointBlock(indices,native.rows))
        end
    end
    bytes=sum(sizeof(native.rows) for native in values(prepared))
    verified_hoppet_build(buildfile)==verification.receipt || error("HOPPET build changed during point preparation")
    (grid=grid,queries=qs,Q2=scales,Q02=fixed_scalar(Q02,"Q02"),tau=times,quantity=quantity,nf=nf,
        intervals=intervals,construction=construction,blocks=blocks,coefficient_payload_bytes=bytes,row_native_max_scaled=maximum(row_errors))
end

function HoppetPointPlan(buildfile,queries,Q2,Q02,coupling;intervals=1024,nf=5,
                         quantity=:transversity,order=-6,eps=1e-10,du=.0125,construction=:batched_probes)
    data=prepare_hoppet_points(buildfile,queries,Q2,Q02,coupling;intervals=intervals,nf=nf,
        quantity=quantity,order=order,eps=eps,du=du,construction=construction)
    HoppetPointPlan(data...,hoppet_point_snapshot(data...))
end

function HoppetMatchingPlan(buildfile,queries,Q2,Q02,coupling;process,
                            quantity=:KPSY_Hhat3_homogeneous,intervals=1024,nf=5,order=-6,eps=1e-10,du=.0125,
                            construction=:batched_probes)
    quantity in (:transversity,:KPSY_Hhat3_homogeneous) || error("unsupported polarized matching quantity")
    channel=quantity==:transversity ? :transversity : :collins
    matching_channel(process,channel)
    q2=[fixed_scalar(q,"matching Q2") for q in Q2];q0=fixed_scalar(Q02,"matching Q02")
    alpha=coupling isa AsymptoticTwoLoopCoupling ? alpha_s.(Ref(coupling),q2) :
        [alpha_lo(q,q0,fixed_scalar(coupling,"matching alpha0");nf=nf) for q in q2]
    all(a -> isfinite(a) && a>0,alpha) || error("finite positive physical matching coupling required")
    data=prepare_hoppet_points(buildfile,queries,q2,q0,coupling;quantity=quantity,intervals=intervals,
        nf=nf,order=order,eps=eps,du=du,construction=construction,
        matching=(process=process,channel=channel,alpha=alpha))
    HoppetMatchingPlan(data...,process,alpha,hoppet_point_snapshot(data...;process=process,alpha=alpha))
end

function hoppet_point_snapshot(grid,queries,Q2,Q02,tau,quantity,nf,intervals,construction,blocks,payload,row_error;
                               process=nothing,alpha=nothing)
    numerical_snapshot((grid,queries,Q2,Q02,tau,quantity,nf,intervals,construction,
        [(block.indices,block.rows) for block in blocks],payload,row_error,process,alpha))
end
function hoppet_point_signature(plan::HoppetPointPlan)
    hoppet_point_snapshot(plan.grid,plan.queries,plan.Q2,plan.Q02,plan.tau,plan.quantity,plan.nf,
        plan.intervals,plan.construction,plan.blocks,plan.coefficient_payload_bytes,plan.row_native_max_scaled)
end
function hoppet_point_signature(plan::HoppetMatchingPlan)
    hoppet_point_snapshot(plan.grid,plan.queries,plan.Q2,plan.Q02,plan.tau,plan.quantity,plan.nf,
        plan.intervals,plan.construction,plan.blocks,plan.coefficient_payload_bytes,plan.row_native_max_scaled;
        process=plan.process,alpha=plan.alpha)
end

function hoppet_point_valid(plan,initial;endpoint=true)
    hoppet_point_signature(plan)==plan.signature || error("mutated HOPPET point settings/grid")
    size(initial,1)==length(plan.grid) && size(initial,2)>0 && all(isfinite,initial) || error("invalid point-plan initial shape/values")
    endpoint && !all(iszero,initial[end,:]) && error("h(1)=0 is required")
end

function apply_point_rows(plan,initial::AbstractMatrix{T}) where {T<:Real}
    hoppet_point_valid(plan,initial)
    result=zeros(promote_type(T,Float64),length(plan.queries),size(initial,2),length(plan.tau))
    for (i,block) in enumerate(plan.blocks)
        isempty(block.indices) && continue
        # A contiguous Float64 gather lets BLAS share all tangent columns. The
        # high-precision reference retains the independent generic arithmetic.
        selected=T===Float64 ? initial[block.indices,:] : view(initial,block.indices,:)
        result[i,:,:]=permutedims(block.rows*selected)
    end
    result
end
apply_operator(plan::Union{HoppetPointPlan,HoppetMatchingPlan},initial::AbstractMatrix{<:Real})=apply_point_rows(plan,initial)
function apply_operator(plan::Union{HoppetPointPlan,HoppetMatchingPlan},initial::AbstractMatrix{<:ForwardDiff.Dual})
    if initial_primal(zero(eltype(initial))) isa Float64
        return hoppet_dual_batch(x -> apply_operator(plan,x),initial)
    end
    apply_point_rows(plan,initial)
end
apply_operator(plan::Union{HoppetPointPlan,HoppetMatchingPlan},initial::AbstractVector)=dropdims(apply_operator(plan,reshape(initial,:,1));dims=2)

function operator_adjoint(plan::Union{HoppetPointPlan,HoppetMatchingPlan},g::AbstractArray{T,3}) where {T<:Real}
    hoppet_point_signature(plan)==plan.signature || error("mutated HOPPET point settings/grid")
    size(g,1)==length(plan.queries) && size(g,3)==length(plan.tau) && all(isfinite,g) || error("invalid point-plan cotangent")
    result=zeros(promote_type(T,Float64),length(plan.grid),size(g,2))
    for (i,block) in enumerate(plan.blocks)
        isempty(block.indices) && continue
        result[block.indices,:].+=transpose(block.rows)*permutedims(g[i,:,:])
    end
    result[end,:].=zero(T)
    result
end
operator_adjoint(plan::Union{HoppetPointPlan,HoppetMatchingPlan},g::AbstractMatrix)=vec(operator_adjoint(plan,reshape(g,size(g,1),1,size(g,2))))
