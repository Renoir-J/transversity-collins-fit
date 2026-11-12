# Independent coupling/CSS reconstruction checks, separate from production admission.
using CollinsFit, ForwardDiff, LinearAlgebra, QuadGK, JSON, SHA, Dates
include("CollinearReference.jl")
using .CollinearReference

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const CONTRACT = joinpath(ROOT, "config/css-contract.json")
const SOURCES = vcat([@__FILE__, CONTRACT, joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(@__DIR__,"CollinearReference.jl"), joinpath(ROOT,"numerics/Project.toml"),
    joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"), readdir(joinpath(ROOT,"numerics/src"); join=true)))
source_hashes() = Dict(relpath(p, ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function reference_alpha(Q2, nf, lambda)
    setprecision(256) do
        L = log(BigFloat(Q2)/BigFloat(lambda)^2)
        b0 = 11 - 2BigFloat(nf)/3
        b1 = 102 - 38BigFloat(nf)/3
        4BigFloat(pi)/(b0*L) * (1 - b1/b0^2*log(L)/L)
    end
end

function reference_sudakov(model, b, Q2)
    # Use a distinct algebraic mu_b expression and adaptive integration.
    lower2 = 4exp(-2Float64(Base.MathConstants.eulergamma)) * (1/b^2+1/model.bmax^2)
    start, stop = log(lower2), log(Q2)
    function integrand(t)
        a = Float64(reference_alpha(exp(t),model.coupling.nf,model.coupling.lambda_GeV))/pi
        A2 = (2/3) * (3*(67/18-pi^2/6) - 5model.coupling.nf/9)
        ((4/3)*a+A2*a^2)*(log(Q2)-t) - 2a
    end
    direction = stop >= start ? 1 : -1
    direction * quadgk(integrand,min(start,stop),max(start,stop); rtol=1e-12,atol=1e-13)[1]
end

function exact_two_loop_rg(Q2, Q02, alpha0, nf)
    # Informational comparator with a different definition, not the chosen oracle.
    beta0, beta1 = 11-2nf/3, 102-38nf/3
    c = beta1/beta0
    F(a) = inv(a) + c*log(a/(1+c*a))
    target = F(alpha0/(4pi)) + beta0*log(Q2/Q02)
    low, high = 1e-5, 0.3
    for _ in 1:100
        middle = (low+high)/2
        if F(middle) > target
            low = middle
        else
            high = middle
        end
    end
    4pi*(low+high)/2
end

function check_css()
    sources = source_hashes()
    contract = JSON.parsefile(CONTRACT)
    policy = contract["validation_policy"]
    checks = Dict{String,Any}[]
    check(name,error,tolerance) = push!(checks,Dict("name"=>name,"error"=>error,
        "tolerance"=>tolerance,"passed"=>isfinite(error) && error <= tolerance))
    coupling_cases = []
    for nf in 3:6, q2 in (0.6,2.4,100.0,91.1876^2,1e6)
        coupling = AsymptoticTwoLoopCoupling(;nf=nf)
        alpha = alpha_s(coupling,q2)
        expected = Float64(reference_alpha(q2,nf,coupling.lambda_GeV))
        lo, hi = minmax(log(2.4),log(q2))
        numerical_tau = sign(log(q2/2.4)) * quadgk(t ->
            Float64(reference_alpha(exp(t),nf,coupling.lambda_GeV))/(2pi),lo,hi; atol=1e-14)[1]
        push!(coupling_cases,Dict("nf"=>nf,"Q2"=>q2,"alpha_s"=>alpha,
            "alpha_relative_error"=>abs(alpha/expected-1),
            "time_integral_error"=>abs(evolution_time(coupling,q2,2.4)-numerical_tau)))
    end
    check("coupling versus 256-bit formula",maximum(c["alpha_relative_error"] for c in coupling_cases),
        policy["coupling_max_relative"])
    check("closed evolution time versus independent adaptive integral",maximum(c["time_integral_error"] for c in coupling_cases),
        policy["integrated_coupling_max_abs"])
    coupling = AsymptoticTwoLoopCoupling()
    defaults = KPSYCSS()
    configured = contract["defaults"]
    match = defaults.coupling.nf == contract["coupling"]["nf"] &&
        defaults.coupling.lambda_GeV == contract["coupling"]["Lambda_GeV"] &&
        all(getproperty(defaults,property)==configured[key] for (property,key) in
            ((:bmax,"bmax_inverse_GeV"),(:Q02,"Q02_GeV2"),(:g2,"g2"),(:gq,"gq_GeV2"),(:gh,"gh_GeV2"))) &&
        length(defaults.nodes)==last(policy["orders"])
    check("implemented defaults match the scientific contract",match ? 0.0 : 1.0,0.0)
    alpha0 = alpha_s(coupling,2.4)
    check("paper alpha0 at its printed precision",abs(alpha0-0.327),0.0005)
    check("paper alpha_MZ at its printed precision",abs(alpha_s(coupling,91.1876^2)-0.118),0.0005)
    check("evolution time composition",abs(evolution_time(coupling,100.0,2.4)-
        evolution_time(coupling,16.0,2.4)-evolution_time(coupling,100.0,16.0)),
        policy["integrated_coupling_max_abs"])
    check("nearby-scale evolution time",abs(evolution_time(coupling,2.4*(1+1e-10),2.4)-
        alpha0/(2pi)*log1p(1e-10)),policy["integrated_coupling_max_abs"])

    models = [KPSYCSS(;order=order) for order in policy["orders"]]
    model = last(models)
    c02 = 4exp(-2Float64(Base.MathConstants.eulergamma))
    negative_b = inv(sqrt(100exp(-0.1)/c02-inv(model.bmax^2)))
    locations = [(b,q2) for b in (0.01,0.1,0.7,1.5,5.0,20.0) for q2 in (2.4,10.0,100.0,91.1876^2)]
    push!(locations,(negative_b,100.0))
    reference = [reference_sudakov(model,b,q2) for (b,q2) in locations]
    quadrature = []
    predictions = []
    for setting in models
        values = [sudakov_perturbative(setting,b,q2) for (b,q2) in locations]
        push!(predictions,values)
        push!(quadrature,Dict("nodes"=>length(setting.nodes),"max_abs_error"=>maximum(abs.(values-reference)),
            "predictions"=>values))
    end
    check("fixed-node Sudakov versus adaptive reference",last(quadrature)["max_abs_error"],policy["sudakov_max_abs"])
    check("Sudakov node refinement",maximum(abs.(predictions[end]-predictions[end-1])),
        policy["quadrature_refinement_max_abs"])
    check("Gauss weights",abs(sum(model.weights)-2),policy["quadrature_refinement_max_abs"])
    check("no silent nonnegative Sudakov clamp",last(last(predictions)) < 0 ? 0.0 : 1.0,0.0)
    constant_errors = []
    A, B = CollinsFit.css_coefficients(0.2,5)
    for L in (-3.0,0.1,2.0)
        actual = CollinsFit.css_integral(1.0,exp(L),_ -> (A,B),model.nodes,model.weights)
        push!(constant_errors,abs(actual-(A*L^2/2+B*L)))
    end
    check("constant-coupling signed integral",maximum(constant_errors),policy["sudakov_max_abs"])
    b, q2 = 0.7, 100.0
    step = 1e-5
    numerical = (sudakov_perturbative(model,b,q2*exp(step))-sudakov_perturbative(model,b,q2*exp(-step)))/(2step)
    integral_A = quadgk(t -> first(CollinsFit.css_coefficients(alpha_s(coupling,exp(t)),5)),
        log(css_mu_b2(model,b)),log(q2); atol=1e-13)[1]
    expected = integral_A + last(CollinsFit.css_coefficients(alpha_s(coupling,q2),5))
    check("Sudakov derivative with respect to ln Q2",abs(numerical-expected),policy["derivative_max_scaled"])

    np_errors = []
    gc = 0.02 # analytic test parameter, not a reported fit result
    for b in (0.01,0.7,5.0,20.0), z in (0.2,0.4,0.8), q2 in (2.4,100.0)
        push!(np_errors,abs(sidis_collins_snp(model,b,q2,z,gc)-sidis_snp(model,b,q2,z)+gc*b^2/z^2))
        push!(np_errors,abs(sia_collins_snp(model,b,q2,z,0.5,gc)-sia_snp(model,b,q2,z,0.5)+gc*b^2*(1/z^2+4)))
        push!(np_errors,abs(sia_snp(model,b,q2,z,0.5)-sia_snp(model,b,q2,0.5,z)))
    end
    check("nonperturbative z factors, channels and exchange",maximum(np_errors),policy["nonperturbative_identity_max_abs"])
    check("no universal logarithm at Q0",abs(sidis_snp(model,0.7,2.4,0.4)-
        (model.gq+model.gh/0.4^2)*0.7^2),policy["nonperturbative_identity_max_abs"])
    b, q2, z1, z2 = 1.2, 100.0, 0.4, 0.55
    perturbative = sudakov_perturbative(model,b,q2)
    f(g) = exp(-perturbative-sidis_collins_snp(model,b,q2,z1,g)) +
           0.7exp(-perturbative-sia_collins_snp(model,b,q2,z1,z2,g))
    s1, s2 = b^2/z1^2, b^2*(1/z1^2+1/z2^2)
    v1 = exp(-perturbative-sidis_collins_snp(model,b,q2,z1,gc))
    v2 = 0.7exp(-perturbative-sia_collins_snp(model,b,q2,z1,z2,gc))
    gradient = ForwardDiff.derivative(f,gc)
    hessian = ForwardDiff.derivative(g -> ForwardDiff.derivative(f,g),gc)
    check("Collins-width gradient",abs(gradient/(s1*v1+s2*v2)-1),policy["derivative_max_scaled"])
    check("Collins-width second derivative",abs(hessian/(s1^2*v1+s2^2*v2)-1),policy["derivative_max_scaled"])
    finite_difference = [Dict("step"=>h,"error"=>abs(((f(gc+h)-f(gc-h))/(2h))/gradient-1)) for h in (1e-3,1e-4,1e-5,1e-6)]
    plateau = any(finite_difference[i]["error"] <= policy["derivative_max_scaled"] &&
        finite_difference[i+1]["error"] <= policy["derivative_max_scaled"] for i in 1:length(finite_difference)-1)
    check("two adjacent Collins-width finite-difference steps",plateau ? 0.0 : 1.0,0.0)

    # The Mellin reference uses a separately integrated time, not the closed core expression.
    tau_reference = quadgk(t -> Float64(reference_alpha(exp(t),5,0.225))/(2pi),log(2.4),log(100.0); atol=1e-14)[1]
    xs = [0.005,0.02,0.08,0.2,0.5,0.8,0.95]
    expected = [mellin_polynomial_tau(x,tau_reference) for x in xs]
    contour = [mellin_polynomial_tau(x,tau_reference;angle=2pi/3) for x in xs]
    check("two-loop Mellin contour refinement",maximum(abs.(expected-contour)),policy["sudakov_max_abs"])
    refinements = []
    for n in policy["collinear_grids"]
        grid = validation_grid(n,1e-6)
        operator = collinear_lo_operator(grid,100.0,2.4,coupling)
        values = apply_operator(operator,grid.*(1 .-grid).^3)
        actual = [interpolate_grid(grid,values,x) for x in xs]
        moments = [numerical_moment(grid,values,k) for k in 1:3]
        target = [6/prod(k+j for j in 1:4)*exp(tau_reference*transversity_moment(k)) for k in 1:3]
        entry = Dict("nodes"=>n,"mellin_scaled_error"=>max_scaled(actual,expected),
            "moment_relative_error"=>maximum(abs.(moments./target.-1)))
        if !isempty(refinements)
            check("two-loop Mellin refinement to $n",entry["mellin_scaled_error"]/last(refinements)["mellin_scaled_error"],policy["refinement_error_ratio_max"])
            check("two-loop moment refinement to $n",entry["moment_relative_error"]/last(refinements)["moment_relative_error"],policy["refinement_error_ratio_max"])
        end
        push!(refinements,entry)
    end
    check("two-loop collinear versus independent Mellin",last(refinements)["mellin_scaled_error"],policy["collinear_max_scaled"])
    check("two-loop collinear moments",last(refinements)["moment_relative_error"],policy["collinear_max_moment_relative"])
    check("source hashes stable during job",source_hashes()==sources ? 0.0 : 1.0,0.0)
    Dict("scope"=>"reconstructed_theory","production_ready"=>false,"original_KPSY_reproduction"=>false,
        "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL", "checks"=>checks,"sources"=>sources,
        "contract_id"=>contract["id"],"contract_sha256"=>sources[relpath(CONTRACT,ROOT)],
        "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),"coupling_cases"=>coupling_cases,
        "settings"=>Dict("coupling_model"=>"asymptotic_two_loop","nf"=>model.coupling.nf,
            "Lambda_GeV"=>model.coupling.lambda_GeV,"bmax_inverse_GeV"=>model.bmax,"Q02_GeV2"=>model.Q02,
            "g2"=>model.g2,"gq_GeV2"=>model.gq,"gh_GeV2"=>model.gh,"test_gc_GeV2"=>gc,
            "quadrature_nodes"=>length(model.nodes),"mu_b_cap_at_Q"=>false),
        "defined_model_comparison_informational"=>Dict("asymptotic_alpha0"=>alpha0,
            "asymptotic_alpha_MZ"=>alpha_s(coupling,91.1876^2),
            "exact_RG_alpha_MZ_same_alpha0"=>exact_two_loop_rg(91.1876^2,2.4,alpha0,5),
            "same_contract"=>false),"sudakov_locations"=>locations,"sudakov_reference"=>reference,
        "quadrature"=>quadrature,"collinear_refinements"=>refinements,"finite_difference"=>finite_difference,
        "input_scope"=>"analytic inputs only; no experimental fit or production admission")
end

output = isempty(ARGS) ? joinpath(ROOT,"runs/css/validation.json") : ARGS[1]
try
    report = check_css()
    writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "failed"=>[c for c in report["checks"] if !c["passed"]])))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
