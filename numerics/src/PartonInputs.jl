# Reconstructed input interfaces. These do not admit production physics jobs.
const PARTON_MANIFEST = normpath(joinpath(@__DIR__, "../../config/parton-sources.json"))
const PARTON_NATIVE_LOCK = ReentrantLock()
parton_error(message) = String(message[1:something(findfirst(iszero,message),length(message)+1)-1])

function verified_parton_asset(directory, id)
    sources = JSON.parsefile(PARTON_MANIFEST)["sources"]
    matches = filter(s -> s["id"] == id, sources)
    length(matches) == 1 || error("unknown or duplicate parton source")
    source = only(matches)
    path = joinpath(directory, source["group"], source["filename"])
    isfile(path) && bytes2hex(sha256(read(path))) == source["sha256"] ||
        error("parton source checksum mismatch: $id")
    path
end

function verify_parton_build(path)
    receipt = JSON.parsefile(path)
    receipt["scope"] == "reconstructed_theory" && !receipt["production_ready"] &&
        receipt["status"] == "BUILT_NOT_VALIDATED" && receipt["LHAPDF_version"] == "6.5.6" ||
        error("incorrect parton build scope/version")
    for key in ("source_hashes", "output_hashes")
        !isempty(receipt[key]) || error("missing parton build hashes")
        for (file, hash) in receipt[key]
            isfile(file) && bytes2hex(sha256(read(file))) == hash || error("stale parton build: $file")
        end
    end
    receipt
end

mutable struct CT10Input
    handle::Ptr{Nothing}
    context::Ptr{Nothing}
    batch::Ptr{Nothing}
    destroy::Ptr{Nothing}
    limits::NTuple{4,Float64}
    policy::Symbol
end

function CT10Input(buildfile; policy=:strict)
    policy in (:strict, :reference) || error("choose strict or explicit reference extrapolation")
    receipt = verify_parton_build(buildfile)
    root = realpath(receipt["input_root"])
    for id in ("CT10_info", "CT10_central")
        verified_parton_asset(root,id)
    end
    handle = Libdl.dlopen(receipt["ct10_library"])
    context = C_NULL
    lock(PARTON_NATIVE_LOCK) do
        try
            version = unsafe_string(ccall(Libdl.dlsym(handle,:collins_lha_version),Cstring,()))
            version == "6.5.6" || error("wrong LHAPDF library version")
            message = zeros(UInt8,2048)
            context = ccall(Libdl.dlsym(handle,:collins_lha_open),Ptr{Cvoid},
                (Cstring,Cint,Ptr{UInt8},Csize_t),root,policy == :strict ? 0 : 1,message,length(message))
            context != C_NULL || error(parton_error(message))
            limits = zeros(4)
            ccall(Libdl.dlsym(handle,:collins_lha_limits),Cint,(Ptr{Cvoid},Ptr{Cdouble}),context,limits) == 0 ||
                error("missing CT10 grid limits")
            all(isapprox.(limits,[1e-8,1.0,1.3^2,1e10];rtol=1e-12)) || error("unexpected CT10 grid limits")
            input = CT10Input(handle,context,Libdl.dlsym(handle,:collins_lha_batch),
                Libdl.dlsym(handle,:collins_lha_close),Tuple(limits),policy)
            finalizer(close,input)
            input
        catch
            context != C_NULL && ccall(Libdl.dlsym(handle,:collins_lha_close),Cvoid,(Ptr{Cvoid},),context)
            Libdl.dlclose(handle)
            rethrow()
        end
    end
end

function Base.close(input::CT10Input)
    lock(PARTON_NATIVE_LOCK) do
        if input.context != C_NULL
            ccall(input.destroy,Cvoid,(Ptr{Cvoid},),input.context)
            input.context = C_NULL
            Libdl.dlclose(input.handle)
            input.handle = C_NULL
        end
    end
    nothing
end
Base.isopen(input::CT10Input) = input.context != C_NULL

"CT10 x*f, rows=PDG flavors, columns=queries; fixed coordinates, no AD through C++."
function ct10_batch(input::CT10Input,pids,x,Q2)
    length(x) == length(Q2) || error("CT10 query size mismatch")
    !isempty(pids) && all(p -> p isa Integer && !(p isa Bool) && (p == 21 || 1 <= abs(p) <= 5),pids) ||
        error("unsupported CT10 flavor")
    flavors = Cint.(pids)
    xs = [fixed_scalar(v,"x") for v in x]
    qs = [fixed_scalar(v,"Q2") for v in Q2]
    all(v -> 0 < v <= 1,xs) && all(>(0),qs) || error("invalid parton coordinates")
    result = Matrix{Float64}(undef,length(flavors),length(xs))
    message = zeros(UInt8,2048)
    lock(PARTON_NATIVE_LOCK) do
        isopen(input) || error("closed CT10 input")
        status = GC.@preserve input ccall(input.batch,Cint,
            (Ptr{Cvoid},Ptr{Cint},Cint,Ptr{Cdouble},Ptr{Cdouble},Csize_t,Ptr{Cdouble},Ptr{UInt8},Csize_t),
            input.context,flavors,length(flavors),xs,qs,length(xs),result,message,length(message))
        status == 0 || error(parton_error(message))
    end
    result
end
parton_pdf(input::CT10Input,pid,x,Q2) = only(ct10_batch(input,[pid],[x],[Q2]))/x

const DSSV_POWERS = [(1.,3.),(1.,4.),(.5,8.),(.5,8.),(.5,8.),(.5,5.)]
# Source **0.3 is default REAL, even though its base/result are DOUBLE PRECISION.
const DSS_POWER_03 = Float64(Float32(.3))
const DSS_POWERS = [(.5,4.),(.5,4.),(.5,4.),(DSS_POWER_03,7.),(DSS_POWER_03,7.),
    (DSS_POWER_03,4.),(.5,4.),(.5,4.),(.5,4.)]

function fortran_real(token)
    occursin(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?$",token) ||
        error("unsupported Fortran DATA literal")
    occursin(r"[dD]",token) ? parse(Float64,replace(lowercase(token),'d'=>'e')) : Float64(parse(Float32,token))
end

function legacy_axis(code,name)
    # Fixed-form columns 1..6 contain comments/labels/continuation, not data.
    body = join([line[7:end] for line in readlines(code) if length(line)>=7 &&
        !(first(line) in ('c','C','*','!'))],"\n")
    matched = match(Regex("(?is)\\bDATA\\s+"*name*"\\s*/([^/]+)/"),body)
    isnothing(matched) && error("missing legacy DATA axis")
    tokens = split(replace(matched.captures[1],r"\s+"=>""),',';keepempty=false)
    values = fortran_real.(tokens)
    all(isfinite,values) && all(diff(values).>0) || error("invalid legacy DATA axis")
    values
end

struct LegacyPartonInput
    kind::Symbol
    x::Vector{Float64}
    q2::Vector{Float64}
    lx::Vector{Float64}
    lq::Vector{Float64}
    table::Array{Float64,3} # component, x, Q²
    powers::Vector{Tuple{Float64,Float64}}
    xmin::Float64
end

function legacy_table(kind,x,q2,raw,powers,xmin)
    size(raw) == (length(powers),length(q2),length(x)-1) || error("legacy grid dimensions")
    all(isfinite,raw) || error("nonfinite legacy grid")
    table = zeros(length(powers),length(x),length(q2))
    for ix in 1:length(x)-1, iq in eachindex(q2), f in eachindex(powers)
        a,b = powers[f]
        table[f,ix,iq] = raw[f,iq,ix]/(x[ix]^a*(1-x[ix])^b)
    end
    LegacyPartonInput(kind,copy(x),copy(q2),log.(x),log.(q2),table,copy(powers),xmin)
end

function load_legacy(directory,kind)
    isdssv = kind == :DSSV08
    kind in (:DSSV08,:DSS14) || error("unknown legacy parton set")
    id = isdssv ? "DSSV08" : "DSS14"
    code = verified_parton_asset(directory,id*"_code")
    grid = verified_parton_asset(directory,id*"_grid")
    x,q2 = legacy_axis(code,"XB"),legacy_axis(code,"QS")
    length(x)==47 && length(q2)==(isdssv ? 30 : 24) || error("wrong legacy DATA axis size")
    powers = isdssv ? DSSV_POWERS : DSS_POWERS
    lines = filter(!isempty,strip.(readlines(grid)))
    length(lines) == (length(x)-1)*length(q2) || error("wrong legacy grid row count")
    values = Float64[]
    for line in readlines(grid)
        isempty(strip(line)) && continue
        fields = if isdssv
            split(line)
        else
            ncodeunits(line) == 90 || error("invalid DSS14 fixed-width row")
            [line[i:i+9] for i in 1:10:90]
        end
        length(fields) == length(powers) || error("invalid legacy row width")
        append!(values,parse.(Float64,fields))
    end
    raw = reshape(values,length(powers),length(q2),length(x)-1)
    legacy_table(kind,x,q2,raw,powers,isdssv ? first(x) : .05)
end
load_dssv08(directory) = load_legacy(directory,:DSSV08)
load_dss14(directory) = load_legacy(directory,:DSS14)

function legacy_coordinates(input,x,Q2,policy)
    policy in (:strict,:reference) || error("unknown parton extrapolation policy")
    x,q2 = fixed_scalar(x,"x or z"),fixed_scalar(Q2,"Q2")
    input.xmin <= x <= 1 && q2 > 0 || error("invalid parton coordinate domain")
    if policy == :strict || input.kind == :DSSV08
        first(input.q2) <= q2 <= last(input.q2) || error("parton scale outside strict range")
    end
    x,q2
end

"Legacy x*Δq or z*D components; DSS14 returns its raw total/valence basis here."
function weighted_components(input::LegacyPartonInput,x,Q2;policy=:strict)
    x,q2 = legacy_coordinates(input,x,Q2,policy)
    u,v = log(x),log(q2)
    i = clamp(searchsortedlast(input.lx,u),1,length(input.lx)-1)
    j = clamp(searchsortedlast(input.lq,v),1,length(input.lq)-1)
    a = (u-input.lx[i])/(input.lx[i+1]-input.lx[i])
    b = (v-input.lq[j])/(input.lq[j+1]-input.lq[j])
    t = input.table
    interpolated = (1-a)*(1-b)*t[:,i,j] + a*(1-b)*t[:,i+1,j] +
        (1-a)*b*t[:,i,j+1] + a*b*t[:,i+1,j+1]
    [interpolated[f]*x^p*(1-x)^r for (f,(p,r)) in enumerate(input.powers)]
end

"Unweighted DSSV08 u,d,ubar,dbar,g. Strange/heavy mapping is deliberately unadmitted."
function parton_helicity(input::LegacyPartonInput,pid,x,Q2)
    input.kind == :DSSV08 || error("helicity requires DSSV08")
    pid isa Integer && !(pid isa Bool) && pid in (2,1,-2,-1,21) || error("unsupported DSSV08 helicity flavor")
    uv,dv,ub,db,_,g = weighted_components(input,x,Q2)
    (pid == 2 ? uv+ub : pid == 1 ? dv+db : pid == -2 ? ub : pid == -1 ? db : g)/x
end

function pion_components(raw,charge)
    charge isa Integer && !(charge isa Bool) && charge in (-1,0,1) || error("invalid pion charge")
    ut,dt,st,ct,bt,g,uv,dv,sv = raw
    s = charge
    [(ut+s*uv)/2,(ut-s*uv)/2,(dt+s*dv)/2,(dt-s*dv)/2,
     (st+s*sv)/2,(st-s*sv)/2,ct/2,bt/2,g]
end

"Unweighted D_i^π using PDG flavors; c=cbar and b=bbar are each half the grid total."
function pion_ff(input::LegacyPartonInput,pid,charge,z,Q2;policy=:strict)
    input.kind == :DSS14 || error("pion FF requires DSS14")
    pid isa Integer && !(pid isa Bool) && (pid == 21 || 1 <= abs(pid) <= 5) || error("unsupported DSS14 flavor")
    components = pion_components(weighted_components(input,z,Q2;policy=policy),charge)
    index = pid == 21 ? 9 : abs(pid) == 5 ? 8 : abs(pid) == 4 ? 7 :
        pid == 2 ? 1 : pid == -2 ? 2 : pid == 1 ? 3 : pid == -1 ? 4 : pid == 3 ? 5 : 6
    components[index]/z
end

function parton_domain(input::CT10Input,x,Q2)
    x,q2 = fixed_scalar(x,"x"),fixed_scalar(Q2,"Q2")
    lo,hi,qlo,qhi = input.limits
    Dict("within_grid"=>lo <= x <= hi && qlo <= q2 <= qhi,
         "policy"=>string(input.policy),"production_ready"=>false)
end
function parton_domain(input::LegacyPartonInput,x,Q2;policy=:strict)
    policy in (:strict,:reference) || error("unknown parton extrapolation policy")
    x,q2 = fixed_scalar(x,"x or z"),fixed_scalar(Q2,"Q2")
    Dict("within_grid"=>input.xmin <= x <= 1 && first(input.q2) <= q2 <= last(input.q2),
         "policy"=>string(policy),"production_ready"=>false)
end
