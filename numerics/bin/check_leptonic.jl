# Component reconstruction only: no experimental data, historical fit or production admission.
using CollinsFit, ForwardDiff, LinearAlgebra, QuadGK, SpecialFunctions, JSON, SHA, Dates
include("LeptonicReference.jl")
using .LeptonicReference

const ROOT = normpath(joinpath(@__DIR__,"../.."))
const CONTRACT = joinpath(ROOT,"config/leptonic-contract.json")
const SOURCES = vcat([@__FILE__,CONTRACT,joinpath(@__DIR__,"LeptonicReference.jl"),
    joinpath(ROOT,"config/css-contract.json"),joinpath(ROOT,"numerics/Project.toml"),
    joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes() = Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function radial_edges(B)
    vcat([0.,1e-7,1e-6,1e-5,1e-4,.001,.005,.02,.05,.1,.2,.4,.7,1.,1.5,2.,3.,4.,6.,8.],
         filter(<=(B),[12.,16.,20.,24.,32.]))
end

function check_leptonic()
    sources = source_hashes()
    contract = JSON.parsefile(CONTRACT)
    policy = contract["validation_policy"]
    floor = policy["relative_floor"]
    scaled(a,b) = maximum(abs.(a.-b)./max.(abs.(b),floor))
    checks = Dict{String,Any}[]
    check(name,value,tolerance) = push!(checks,Dict("name"=>name,"error"=>value,
        "tolerance"=>tolerance,"passed"=>isfinite(value) && value<=tolerance))
    # Coefficient Mellin moments have distinct harmonic/trigamma closed forms.
    moment_errors = Float64[]
    CF = 4/3
    for N in 3:7, channel in (:pdf,:ff,:collins,:transversity)
        expected = if channel == :pdf
            [CF/(2N*(N+1)),.5*(1/(N+1)-1/(N+2))]
        elseif channel == :ff
            [CF*(1/(2N*(N+1))-trigamma(N)-trigamma(N+2)),
             CF*(1/(2(N+1))-2/(N-1)^2+2/N^2-1/(N+1)^2)]
        elseif channel == :collins
            [-2CF*trigamma(N+1),0.]
        else
            [0.,0.]
        end
        actual = [quadgk(v -> v^(N-1)*CollinsFit.matching_regular(channel,v)[f],0.,1.;
            rtol=1e-12,atol=1e-13)[1] for f in 1:2]
        push!(moment_errors,maximum(abs.(actual.-expected)))
    end
    check("regular Wilson coefficients versus closed Mellin moments",maximum(moment_errors),policy["kernel_moment_max_abs"])
    endpoint = [abs(first(CollinsFit.matching_regular(c,1.))+8/3) for c in (:ff,:collins)]
    push!(endpoint,abs(last(CollinsFit.matching_regular(:ff,1.))-2/3))
    push!(endpoint,abs(CollinsFit.matching_delta(:SIDIS,:transversity)+8/3))
    push!(endpoint,abs(CollinsFit.matching_delta(:SIA,:collins)-(pi^2-8)/3))
    check("regular endpoint and process-specific delta coefficients",maximum(endpoint),policy["endpoint_max_abs"])
    q(y) = y^.4*(1-y)^3
    g(y) = y^.7*(1-y)^2
    cases = [(p,c,x,a) for (p,c) in ((:SIDIS,:pdf),(:SIDIS,:transversity),(:SIDIS,:ff),
        (:SIDIS,:collins),(:SIA,:ff),(:SIA,:collins)) for x in (1e-5,.001,.03,.2,.7,.99,1.) for a in (0.,.3)]
    reference = [reference_match(p,c,x,a,q,g) for (p,c,x,a) in cases]
    matching = []
    previous = nothing
    for order in policy["matching_orders"]
        values = [begin
            rule = matching_rule(p,c,x,a;order=order)
            match_leg(rule,q.(rule.coordinates);gluon=c in (:pdf,:ff) ? g.(rule.coordinates) : nothing)
        end for (p,c,x,a) in cases]
        push!(matching,Dict("order"=>order,"max_scaled_error"=>scaled(values,reference)))
        if order==last(policy["matching_orders"])
            check("matching convolution versus 256-bit different-variable quadrature",scaled(values,reference),policy["matching_max_scaled"])
            check("matching node refinement",scaled(values,previous),policy["node_refinement_max_scaled"])
        end
        previous = values
    end
    polynomial = []
    for x in (.001,.2,.8), p in 0:3
        J(m) = m == -1 ? -log(x) : -expm1((m+1)*log(x))/(m+1)
        rule = matching_rule(:SIDIS,:pdf,x,.3)
        expected = (1-(8/3)*.3/pi)*x^p + .3/pi*((2/3)*(J(p-1)-x*J(p-2))+
            .5*(x*J(p-1)-x^2*J(p-2))) # gluon y^(p+1)
        actual = match_leg(rule,rule.coordinates.^p;gluon=rule.coordinates.^(p+1))
        push!(polynomial,abs(actual-expected))
    end
    check("PDF convolution Jacobian and limits versus polynomial antiderivatives",maximum(polynomial),policy["kernel_moment_max_abs"])
    gaussian_refinements = []
    for setting in policy["Gaussian_refinements"]
        radial = RadialRule(collect(0.:setting["panel_width"]:16.);order=setting["order"])
        errors = []
        for A in (.4,.9,2.1), P in (0.,.02,.3,.9,2.), z in (.2,.5,.8)
            amp = exp.(-A.*radial.b.^2)
            k = P/z
            sidis,sia = sidis_fourier(radial,P,z,amp,amp),sia_fourier(radial,P,z,.6,amp,amp)
            expected = [gaussian_hankel(0,k,A)/z^2,-gaussian_hankel(1,k,A)/(2z^3),
                gaussian_hankel(0,k,A)/z^2,gaussian_hankel(2,k,A)/(4z^3*.6)]
            push!(errors,scaled([sidis.UU,sidis.Collins,sia.UU,sia.Collins],expected))
        end
        push!(gaussian_refinements,Dict("order"=>setting["order"],"panel_width"=>setting["panel_width"],
            "max_scaled_error"=>maximum(errors),"cases"=>length(errors),
            "passes_target"=>maximum(errors)<=policy["gaussian_transform_max_scaled"]))
    end
    check("all structure-function Fourier prefactors versus exact Gaussians",last(gaussian_refinements)["max_scaled_error"],policy["gaussian_transform_max_scaled"])
    angular = [abs(besselj(n,t)-angular_bessel(n,t)) for n in 0:2 for t in (0.,.03,.7,3.,12.)]
    check("Bessel harmonics versus independent angular Fourier projections",maximum(angular),policy["angular_fourier_max_abs"])

    model,gc = KPSYCSS(),.012
    ids = [2,-2,1,-1,3,-3]
    A(q,b) = (1+.1q)*(1+.03b^2)*exp(-.02b^2)
    D(q,b) = (.8-.05q)*(1+.02b)*exp(-.01b^2)
    H(q,b) = (.2+.03q)*exp(-.03b^2)
    K(q,b) = (-.1+.04q)*(1+.01b^2)*exp(-.02b^2)
    locations = [(:SIDIS,2.4,.1,.3,nothing),(:SIDIS,10.,.4,.6,nothing),
        (:SIDIS,100.,.7,.8,nothing),(:SIA,2.4,.1,.3,.5),(:SIA,110.,.4,.6,.7),(:SIA,100.,.7,.8,.8)]
    css_reference = [reference_structure(p,model,Q2,P,z,z2,ids,A,D,H,K,gc,32.) for (p,Q2,P,z,z2) in locations]
    assembly = []
    node_values = []
    for order in policy["radial_orders"]
        radial = RadialRule(radial_edges(32.);order=order)
        values = [begin
            plan = leptonic_plan(p,model,radial,Q2,P,z;z2=z2,flavors=ids)
            matrices = [[f(q,b) for q in ids,b in radial.b] for f in (A,D,H,K)]
            result = leptonic_structure(plan,matrices...;gc=gc)
            [result.UU,result.Collins]
        end for (p,Q2,P,z,z2) in locations]
        push!(node_values,reduce(vcat,values))
        push!(assembly,Dict("order"=>order,"max_scaled_error"=>scaled(last(node_values),
            reduce(vcat,[[r.UU,r.Collins] for r in css_reference]))))
    end
    check("matched flavor/CSS assembly versus independent adaptive radial reference",last(assembly)["max_scaled_error"],policy["CSS_assembly_max_scaled"])
    check("radial node refinement",scaled(node_values[end],node_values[end-1]),policy["node_refinement_max_scaled"])
    tails = []
    for B in policy["radial_bounds"]
        radial = RadialRule(radial_edges(B);order=last(policy["radial_orders"]))
        values = [begin
            plan = leptonic_plan(p,model,radial,Q2,P,z;z2=z2,flavors=ids)
            result = leptonic_structure(plan,[[f(q,b) for q in ids,b in radial.b] for f in (A,D,H,K)]...;gc=gc)
            [result.UU,result.Collins]
        end for (p,Q2,P,z,z2) in locations]
        push!(tails,Dict("B_inverse_GeV"=>B,"values"=>reduce(vcat,values)))
    end
    check("finite radial bound extension 24 to 32",scaled(tails[end]["values"],tails[end-1]["values"]),policy["tail_extension_max_scaled"])

    radial = RadialRule(radial_edges(24.);order=32)
    derivatives = []
    for process in (:SIDIS,:SIA)
        plan = leptonic_plan(process,model,radial,10.,.35,.4;z2=process==:SIA ? .6 : nothing,flavors=ids)
        r = matching_rule(:SIDIS,:collins,.3,.3)
        function objective(t)
            h(y) = exp(t[1])*y^.4*(1-y)^(2+exp(t[2]))
            normalization = match_leg(r,h.(r.coordinates))
            first = [normalization*H(q,b) for q in ids,b in radial.b]
            second = [(1+t[2]*b/(1+b))*K(q,b) for q in ids,b in radial.b]
            width_gc = model.gh-.03exp(t[3])
            leptonic_collins(plan,first,second,width_gc)
        end
        t,direction = [.2,-.3,.1],normalize([.3,-.5,.4])
        gradient = ForwardDiff.gradient(objective,t)
        Hmat = ForwardDiff.hessian(objective,t)
        fderrors,hvperrors = Float64[],Float64[]
        for step in (1e-3,1e-4,1e-5)
            fd = (objective(t+step*direction)-objective(t-step*direction))/(2step)
            push!(fderrors,scaled([dot(gradient,direction)],[fd]))
            fd_hvp = (ForwardDiff.gradient(objective,t+step*direction)-
                ForwardDiff.gradient(objective,t-step*direction))/(2step)
            push!(hvperrors,scaled(Hmat*direction,fd_hvp))
        end
        check("$process full matching/width gradient plateau",maximum(fderrors[2:3]),policy["derivative_max_scaled"])
        check("$process full matching/width HVP plateau",maximum(hvperrors[2:3]),policy["derivative_max_scaled"])
        check("$process Hessian symmetry",maximum(abs.(Hmat.-transpose(Hmat))),policy["endpoint_max_abs"])
        push!(derivatives,Dict("process"=>string(process),"gradient_errors"=>fderrors,"hvp_errors"=>hvperrors))
    end
    check("source hashes stable during job",source_hashes()==sources ? 0. : 1.,0.)
    Dict("scope"=>"reconstructed_theory","production_ready"=>false,"original_KPSY_reproduction"=>false,
        "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
        "contract_id"=>contract["id"],"julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),
        "matching_cases"=>length(cases),"matching_refinements"=>matching,
        "Gaussian_refinements"=>gaussian_refinements,"CSS_locations"=>locations,"CSS_refinements"=>assembly,
        "finite_bound_refinements"=>tails,"derivatives"=>derivatives,
        "input_scope"=>"analytic component inputs only; no experiment or fitted initial distributions",
        "limitations"=>contract["limitations"])
end

output = isempty(ARGS) ? joinpath(ROOT,"runs/leptonic/validation.json") : ARGS[1]
try
    report = check_leptonic()
    writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
