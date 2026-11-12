const INITIAL_CONTRACT = normpath(joinpath(@__DIR__,"../../config/initial-contract.json"))
const INITIAL_PARAMETER_NAMES = (:Nh_u,:Nh_d,:a_u,:a_d,:b_u,:b_d,:Nc_u,:Nc_d,
    :alpha_u,:alpha_d,:beta_u,:beta_d,:gc)
const INITIAL_FLAVORS = [-5,-4,-3,-2,-1,1,2,3,4,5]

struct KPSYParameters{T<:Real}
    values::NTuple{13,T}
    gh::Float64
end

"Direct physical coordinates; closed Nh/alpha/beta boundaries are retained exactly."
function KPSYParameters(theta::AbstractVector;gh=.042)
    length(theta)==13 && all(v -> v isa Real && !(v isa Bool) && isfinite(v),theta) ||
        error("13 finite real KPSY parameters required")
    width = fixed_scalar(gh,"fixed gh")
    width>0 || error("positive gh required")
    all(i -> abs(theta[i])<=1,1:2) || error("Nh must lie in the closed unit interval")
    all(i -> theta[i]>=0,(3,4,5,6,9,10,11,12)) || error("nonnegative shape parameters required")
    theta[13]<width || error("gc must be below gh")
    values = Base.tail(promote(1.,theta...))
    KPSYParameters{typeof(first(values))}(values,width)
end

"Rounded Table 1 values for reconstruction checks, never a new optimizer result."
function kpsy_published_parameters()
    KPSYParameters(Float64.(JSON.parsefile(INITIAL_CONTRACT)["published_parameters"]))
end

function initial_fraction(x)
    value = fixed_scalar(x,"initial fraction")
    0<value<=1 || error("initial fraction must lie in (0,1]")
    value
end

initial_primal(t::Real) = t
initial_primal(t::ForwardDiff.Dual) = initial_primal(ForwardDiff.value(t))

function entropy_term(t)
    if iszero(initial_primal(t))
        t isa ForwardDiff.Dual && error("normalization derivative at a=0 or b=0 is singular")
        return zero(t)
    end
    t*log(t)
end

function peak_shape(x,a,b)
    exp(a*log(x)+b*log1p(-x)+entropy_term(a+b)-entropy_term(a)-entropy_term(b))
end

function initial_flavor(pid)
    pid isa Integer && !(pid isa Bool) && 1<=abs(pid)<=5 || error("quark PDG flavor required")
    pid
end

"Unweighted h1 at Q0; the historical model sets sea and heavy transversity to zero."
function initial_transversity(p::KPSYParameters,pid,x,soffer)
    initial_flavor(pid)
    x = initial_fraction(x)
    bound = fixed_scalar(soffer,"fixed Soffer input")
    bound>=0 || error("negative Soffer input; no clipping is permitted")
    (x==1 || !(pid in (2,1))) && return zero(p.values[1])
    n,a,b = pid==2 ? (p.values[1],p.values[3],p.values[5]) : (p.values[2],p.values[4],p.values[6])
    n*peak_shape(x,a,b)*bound
end

function initial_shape(z,n,a,b,D)
    z==1 && return zero(n+a+b)
    n*exp(a*log(z)+b*log1p(-z))*D
end

"Hhat3 from the four independent pi+ FF bases (u,d,s,sbar), in that order."
function initial_collins(p::KPSYParameters,pid,charge,z,ff)
    initial_flavor(pid)
    charge isa Integer && !(charge isa Bool) && charge in (-1,0,1) || error("pion charge must be -1/0/+1")
    z = initial_fraction(z)
    length(ff)==4 && all(v -> v isa Union{AbstractFloat,Integer} && !(v isa Bool) && isfinite(v),ff) ||
        error("four fixed finite pi+ FF bases required")
    (z==1 || abs(pid)>=4) && return zero(p.values[1])
    favored = initial_shape(z,p.values[7],p.values[9],p.values[11],ff[1])
    unfavored = initial_shape(z,p.values[8],p.values[10],p.values[12],ff[2])
    if abs(pid)==3
        return initial_shape(z,p.values[8],p.values[10],p.values[12],ff[pid==3 ? 3 : 4])
    end
    charge==0 && return (favored+unfavored)/2
    isfavored = charge==1 ? pid in (2,-1) : pid in (-2,1)
    isfavored ? favored : unfavored
end

function trento_first_moment(hhat,z,mass)
    z = initial_fraction(z)
    mass = fixed_scalar(mass,"hadron mass")
    mass>0 && hhat isa Real && isfinite(hhat) || error("finite real Hhat3 and positive hadron mass required")
    -hhat/(2z*mass)
end

struct InitialConditions
    x::Vector{Float64}
    z::Vector{Float64}
    soffer::Matrix{Float64} # rows u,d
    ff::Matrix{Float64} # rows pi+/u, pi+/d, pi+/s, pi+/sbar
    Q02::Float64
    provenance::Dict{String,Any}
end

function initial_background_grid(n,minimum,transition)
    low_count=n÷4+1;high_count=n-low_count+1
    low=exp.(range(log(minimum),log(transition);length=low_count))
    low[1]=minimum;low[end]=transition
    high=1 .- range(sqrt(1-transition),0.;length=high_count).^2
    vcat(low,high[2:end])
end

"Align the reference grid with locked input knots and declared prediction probes."
function initial_source_grids(pdf::CT10Input,helicity::LegacyPartonInput,ff::LegacyPartonInput;
                              background_nodes=3201,x_probes=Float64[],z_probes=Float64[])
    n=background_nodes
    n isa Integer && !(n isa Bool) && n>=21 || error("at least 21 background nodes required")
    helicity.kind==:DSSV08 && ff.kind==:DSS14 || error("DSSV08/DSS14 source grids required")
    xp=[fixed_scalar(v,"initial x probe") for v in x_probes]
    zp=[fixed_scalar(v,"initial z probe") for v in z_probes]
    all(v -> 1e-5<=v<=1,xp) && all(v -> .05<=v<=1,zp) || error("probe outside strict initial support")
    raw=readlines(verified_parton_asset(pdf.input_root,"CT10_central"))
    pdfx=parse.(Float64,split(raw[4]))
    xs=sort(unique(vcat(initial_background_grid(n,1e-5,.02),filter(v -> 1e-5<=v<=1,pdfx),helicity.x,xp)))
    zs=sort(unique(vcat(initial_background_grid(n,.05,.2),filter(v -> .05<=v<=1,ff.x),zp)))
    (x=xs,z=zs,background_nodes=Int(n),definition="background plus locked input knots and declared probes")
end

function initial_conditions(x,z,soffer,ff;provenance=Dict{String,Any}("input_kind"=>"analytic"))
    xs,zs = collinear_grid(x),collinear_grid(z)
    size(soffer)==(2,length(xs)) && all(isfinite,soffer) && all(>=(0),soffer) ||
        error("finite nonnegative u,d Soffer table required")
    size(ff)==(4,length(zs)) && all(isfinite,ff) || error("four finite pion FF bases required")
    all(iszero,soffer[:,end]) && all(iszero,ff[:,end]) || error("initial input endpoint must be exactly zero")
    InitialConditions(xs,zs,Float64.(soffer),Float64.(ff),2.4,copy(provenance))
end

"Strict u/d Soffer inputs at the contracted initial scale."
function prepare_initial_soffer(pdf::CT10Input,helicity::LegacyPartonInput,x)
    helicity.kind==:DSSV08 || error("DSSV08 input required")
    xs=collinear_grid(x)
    first(xs)>=1e-5 || error("initial grid outside strict input support")
    all(x -> parton_domain(pdf,x,2.4)["within_grid"],xs) || error("initial PDF extrapolation is not admitted")
    values = ct10_batch(pdf,[2,1],xs,fill(2.4,length(xs)))
    soffer = [(values[f,j]/x+parton_helicity(helicity,pid,x,2.4))/2 for (f,pid) in enumerate((2,1)),(j,x) in enumerate(xs)]
    xs,soffer
end

"Strict four pi+ fragmentation bases, without an unnecessary PDF dependency."
function prepare_initial_ff(ff::LegacyPartonInput,z)
    ff.kind==:DSS14 || error("DSS14 input required")
    zs=collinear_grid(z)
    first(zs)>=.05 || error("initial grid outside strict input support")
    bases = hcat([pion_ff(ff,[2,1,3,-3],1,y,2.4) for y in zs]...)
    zs,bases
end

"Load strict Q0 support once; no provider callbacks enter initial parameter evaluation."
function prepare_initial_conditions(pdf::CT10Input,helicity::LegacyPartonInput,ff::LegacyPartonInput,x,z)
    xs,soffer=prepare_initial_soffer(pdf,helicity,x)
    zs,bases=prepare_initial_ff(ff,z)
    initial_conditions(xs,zs,soffer,bases;provenance=Dict{String,Any}("input_kind"=>"pinned_CT10_DSSV08_DSS14",
        "scope"=>"reconstructed_theory","production_ready"=>false,"initial_extrapolation"=>false,
        "input_manifest_sha256"=>bytes2hex(sha256(read(PARTON_MANIFEST)))))
end

function initial_distributions(inputs::InitialConditions,p::KPSYParameters;charge=1)
    h = [initial_transversity(p,pid,x,pid==2 ? inputs.soffer[1,j] : pid==1 ? inputs.soffer[2,j] : 0.)
        for pid in INITIAL_FLAVORS,(j,x) in enumerate(inputs.x)]
    H = [initial_collins(p,pid,charge,z,view(inputs.ff,:,j)) for pid in INITIAL_FLAVORS,(j,z) in enumerate(inputs.z)]
    (h1=h,Hhat3=H,flavors=copy(INITIAL_FLAVORS),gc=p.values[13])
end

function initial_operator(op,grid,Q02,quantity)
    op isa CollinearLOOperator && op.quantity==quantity && op.grid==grid && op.Q02==Q02 &&
        op.nf==5 && op.coupling_model==:asymptotic_two_loop && op.lambda_GeV==.225 ||
        error("initial family requires a matching fixed-grid KPSY two-loop LO operator")
end

"Connect the initial family to separately tagged fixed-grid evolution references."
function evolve_initial_conditions(inputs::InitialConditions,p::KPSYParameters,h_operator,H_operator;charge=1)
    initial_operator(h_operator,inputs.x,inputs.Q02,:transversity)
    initial_operator(H_operator,inputs.z,inputs.Q02,:KPSY_Hhat3_homogeneous)
    h_operator.Q2==H_operator.Q2 || error("h1/Hhat3 output scales must agree")
    initial = initial_distributions(inputs,p;charge=charge)
    (h1=permutedims(apply_operator(h_operator,permutedims(initial.h1))),
     Hhat3=permutedims(apply_operator(H_operator,permutedims(initial.Hhat3))),
     flavors=initial.flavors,gc=initial.gc,Q2=h_operator.Q2)
end
