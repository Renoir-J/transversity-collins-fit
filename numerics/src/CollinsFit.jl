module CollinsFit

using LinearAlgebra, Random, Statistics
using ForwardDiff, Optim, QuadGK, JSON
using SHA, Libdl

export Context, load_context, coefficients, predict, predict!, objective, value_and_gradient!,
       fit, gradient_report, benchmark_report, closure_report, performance_report,
       transversity_moment, alpha_lo, evolve_moment, apply_operator, operator_adjoint,
       CollinearLOOperator, collinear_lo_generator, collinear_lo_operator,
       AsymptoticTwoLoopCoupling, alpha_s, evolution_time, KPSYCSS,
       css_bstar, css_mu_b2, sudakov_perturbative,
       sidis_snp, sidis_collins_snp, sia_snp, sia_collins_snp,
       CT10Input, load_dssv08, load_dss14, weighted_components,
       ct10_batch, parton_pdf, parton_helicity, pion_ff, parton_domain,
       MatchingRule, matching_rule, match_leg, RadialRule, LeptonicPlan,
       sidis_fourier, sia_fourier, leptonic_plan, leptonic_structure,
       leptonic_unpolarized, leptonic_collins,
       prepare_unpolarized_sidis, prepare_unpolarized_sia,
       KPSYParameters, InitialConditions, kpsy_published_parameters,
       initial_transversity, initial_collins, trento_first_moment,
       prepare_initial_conditions, initial_distributions, evolve_initial_conditions,
       initial_source_grids,
       HoppetPlan, HoppetPointPlan, hoppet_evaluate, hoppet_generator,
       profile_nuisance, writejson, run_job

struct Context
    data::Dict{String,Any}
    kernel::Dict{String,Any}
    covariance::Matrix{Float64}
    nuisance::Matrix{Float64}
end

matrix(rows) = isempty(rows) ? zeros(0, 0) : reduce(vcat, permutedims.(Vector{Float64}.(rows)))

function load_context(datafile, kernelfile)
    d, k = JSON.parsefile(datafile), JSON.parsefile(kernelfile)
    d["scope"] in ("synthetic_closure", "production") || error("experimental candidates are not admitted fit data")
    length(d["points"]) == length(k["bins"]) || error("kernel/data size mismatch")
    [p["id"] for p in d["points"]] == [b["point_id"] for b in k["bins"]] || error("kernel point order mismatch")
    d["scope"] == k["scope"] || error("kernel/data scope mismatch")
    C = matrix(d["covariance"])
    isposdef(Symmetric(C)) || error("nonpositive covariance")
    B = matrix(get(d, "nuisance_design", [Float64[] for _ in d["points"]]))
    Context(d, k, C, B)
end

function coefficients(theta, context)
    [parameter["transform"] == "exp" ? exp(theta[i]) :
     parameter["transform"] == "tanh" ? tanh(theta[i]) :
     parameter["transform"] == "identity" ? theta[i] : error("unknown parameter transform")
     for (i, parameter) in enumerate(context.kernel["parameters"])]
end

function predict!(out, theta, context)
    c = coefficients(theta, context)
    for (i, bin) in enumerate(context.kernel["bins"])
        numerator = zero(eltype(c)) + get(bin, "constant", 0.0)
        for term in get(bin, "linear", [])
            numerator += term["weight"] * c[term["index"]]
        end
        for term in get(bin, "bilinear", [])
            numerator += term["weight"] * c[term["left"]] * c[term["right"]]
        end
        denominator = bin["denominator"]
        denominator > 0 || error("nonpositive bin-integrated unpolarized denominator")
        out[i] = bin["sign"] * numerator / denominator
    end
    out
end

predict(theta, context) = predict!(similar(theta, length(context.kernel["bins"])), theta, context)

function profile_nuisance(residual, covariance, B)
    isempty(B) && return (sum(residual .* (covariance \ residual)), eltype(residual)[])
    weighted_B = covariance \ B
    a = (I + transpose(B) * weighted_B) \ (transpose(B) * (covariance \ residual))
    shifted = residual - B * a
    sum(shifted .* (covariance \ shifted)) + sum(abs2, a), a
end

function objective(theta, context; observations=nothing, indices=nothing)
    observed = isnothing(observations) ? [p["value"] for p in context.data["points"]] : observations
    residual = observed - predict(theta, context)
    chosen = isnothing(indices) ? eachindex(residual) : indices
    chi2, _ = profile_nuisance(residual[chosen], context.covariance[chosen, chosen], context.nuisance[chosen, :])
    chi2
end

function value_and_gradient!(gradient, theta, context)
    ForwardDiff.gradient!(gradient, t -> objective(t, context), theta)
    objective(theta, context)
end

apply_operator(U, initial) = U * initial
operator_adjoint(U, cotangent) = transpose(U) * cotangent

"LO spacelike transversity Mellin kernel, d/dlnQ² = αs/(2π) P_N. Fixed nf only."
function transversity_moment(n::Integer)
    n >= 1 || error("positive integer Mellin moment required")
    (4 / 3) * (1.5 - 2 * sum(1 / k for k in 1:n))
end

function alpha_lo(Q2, Q02, alpha0; nf=3)
    Q2 > 0 && Q02 > 0 && alpha0 > 0 && 3 <= nf <= 6 || error("invalid fixed-nf coupling inputs")
    denominator = 1 + alpha0 * (11 - 2nf / 3) * log(Q2 / Q02) / (4pi)
    denominator > 0 || error("Landau-pole domain")
    alpha0 / denominator
end

function evolve_moment(moment, n, Q2, Q02, alpha0; nf=3)
    moment * (alpha_lo(Q2, Q02, alpha0; nf=nf) / alpha0)^(-2transversity_moment(n) / (11 - 2nf / 3))
end

include("RunningCoupling.jl")
include("CollinearLO.jl")
include("CSS.jl")
include("PartonInputs.jl")
include("Matching.jl")
include("Leptonic.jl")
include("LeptonicInputs.jl")
include("InitialConditions.jl")
include("Hoppet.jl")

function fit(context; seed=20261111, starts=4, replicas=0, observations=nothing, indices=nothing)
    rng = MersenneTwister(seed)
    initial = Float64.(context.kernel["initial"])
    f = t -> objective(t, context; observations=observations, indices=indices)
    candidates = []
    for j in 1:starts
        start = j == 1 ? initial : initial + 0.5randn(rng, length(initial))
        gradient! = (g, t) -> ForwardDiff.gradient!(g, f, t)
        result = optimize(f, gradient!, start, LBFGS(), Optim.Options(iterations=2000, g_tol=1e-8))
        push!(candidates, Dict("theta" => Optim.minimizer(result), "chi2" => Optim.minimum(result),
                              "converged" => Optim.converged(result), "iterations" => Optim.iterations(result)))
    end
    sort!(candidates; by=c -> c["chi2"])
    best = candidates[1]
    theta = best["theta"]
    gradient = ForwardDiff.gradient(f, theta)
    accepted = [c for c in candidates if c["converged"] ||
                norm(ForwardDiff.gradient(f, c["theta"]), Inf) < 1e-5]
    agreement = length(accepted) >= min(2, starts) && all(abs(c["chi2"] - best["chi2"]) < 1e-5 for c in accepted)
    converged = (best["converged"] || norm(gradient, Inf) < 1e-5) && agreement
    replica_values = []
    if replicas > 0
        # Both uncorrelated noise and shared Gaussian nuisance fluctuations enter replicas.
        total_covariance = context.covariance + context.nuisance * transpose(context.nuisance)
        L = cholesky(Symmetric(total_covariance)).L
        center = predict(theta, context)
        for j in 1:replicas
            y = center + L * randn(rng, length(center))
            r = fit(context; seed=seed + j, starts=2, observations=y, indices=indices)
            r["converged"] || error("replica $j did not converge")
            push!(replica_values, r["theta"])
        end
    end
    Dict("scope" => context.data["scope"], "sample_kind" => "frequentist_parametric_replicas",
         "theta" => theta, "coefficients" => coefficients(theta, context), "chi2" => best["chi2"],
         "predictions" => predict(theta, context), "converged" => converged,
         "gradient_inf" => norm(gradient, Inf), "starts" => candidates, "multistart_agreement" => agreement,
         "replicas" => replica_values, "seed" => seed,
         "limitations" => ["Replicas are not posterior samples", "Fixed covariance and additive Gaussian nuisances only"])
end

function gradient_report(context; seed=20261111)
    rng = MersenneTwister(seed)
    atol, rtol = 1e-6, 1e-4
    steps = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    comparisons = []
    n = length(context.kernel["initial"])
    for point in 1:12
        theta = Float64.(context.kernel["initial"]) + 0.3randn(rng, n)
        funcs = Function[t -> objective(t, context)]
        append!(funcs, [t -> predict(t, context)[i] for i in eachindex(context.kernel["bins"])])
        directions = [Matrix{Float64}(I, n, n)[:, i] for i in 1:n]
        append!(directions, [normalize(randn(rng, n)) for _ in 1:8])
        for f in funcs, direction in directions
            ad = dot(ForwardDiff.gradient(f, theta), direction)
            refs = [(f(theta + h * direction) - f(theta - h * direction)) / (2h) for h in steps]
            tolerance(x, y) = atol + rtol * max(abs(x), abs(y))
            stable = [j for j in 1:length(refs)-1 if abs(refs[j] - refs[j+1]) <= tolerance(refs[j], refs[j+1])]
            passed = any(abs(ad - refs[j+1]) <= tolerance(ad, refs[j+1]) for j in stable)
            push!(comparisons, passed)
        end
    end
    U, v, w = randn(rng, 5, 5), randn(rng, 5), randn(rng, 5)
    adjoint_error = abs(dot(w, apply_operator(U, v)) - dot(operator_adjoint(U, w), v))
    linear_error = norm(apply_operator(U, 2v) - 2apply_operator(U, v))
    theta = Float64.(context.kernel["initial"])
    H = ForwardDiff.hessian(t -> objective(t, context), theta)
    direction = normalize(randn(rng, n))
    h = 1e-5
    hvp_fd = (ForwardDiff.gradient(t -> objective(t, context), theta + h * direction) -
              ForwardDiff.gradient(t -> objective(t, context), theta - h * direction)) / (2h)
    hvp_error = norm(H * direction - hvp_fd) / max(1.0, norm(hvp_fd))
    passed = all(comparisons) && adjoint_error < 1e-10 && linear_error < 1e-12 &&
             isapprox(H, transpose(H); atol=1e-8, rtol=1e-8) && hvp_error < 1e-4
    Dict("passed" => passed, "scope" => context.data["scope"], "checks" => length(comparisons),
         "failed_derivatives" => count(!, comparisons), "parameter_points" => 12,
         "random_directions_per_point" => 8, "steps" => steps,
         "adjoint_error" => adjoint_error, "linearity_error" => linear_error, "hvp_relative_error" => hvp_error,
         "atol" => atol, "rtol" => rtol, "integration_refinement" => "N/A: exact finite synthetic kernel; production requires separate tests")
end

function benchmark_report(context)
    theta = Float64.(context.kernel["truth_theta"])
    predicted = predict(theta, context)
    golden = Float64.(context.kernel["independent_predictions"])
    absolute_error = maximum(abs.(predicted - golden))
    sigma_error = maximum(abs.(predicted - golden) ./ sqrt.(diag(context.covariance)))
    # Independent adaptive integration of the distribution-valued LO kernel moment.
    moments = [abs(transversity_moment(n) - (4/3) * (1.5 +
                quadgk(z -> 2(z^n - 1) / (1-z), 0.0, 1.0; rtol=1e-11)[1])) for n in 1:8]
    # The first moment evolves with exponent 4/(33-2nf), and is not conserved.
    ratio = evolve_moment(1.0, 1, 100.0, 4.0, 0.3)
    expected = (alpha_lo(100.0, 4.0, 0.3) / 0.3)^(4/27)
    Dict("passed" => absolute_error < 1e-12 && sigma_error < 1e-9 && maximum(moments) < 1e-9 &&
                     abs(ratio - expected) < 1e-14 && ratio < 1,
         "scope" => context.data["scope"], "max_abs_error" => absolute_error,
         "max_sigma_error" => sigma_error, "LO_mellin_integral_error" => maximum(moments),
         "tensor_charge_scale_ratio" => ratio,
         "physics_qualification" => "Synthetic process kernels plus fixed-nf LO analytic checks; no KPSY or STAR reproduction")
end

function performance_report(context; samples=20)
    theta = Float64.(context.kernel["initial"])
    f = t -> objective(t, context)
    cold = @elapsed ForwardDiff.gradient(f, theta)
    f(theta); ForwardDiff.gradient(f, theta)
    times, gradients = Float64[], Float64[]
    for _ in 1:samples
        push!(times, @elapsed f(theta))
        push!(gradients, @elapsed ForwardDiff.gradient(f, theta))
    end
    allocations = @allocated ForwardDiff.gradient(f, theta)
    median_time = median(gradients)
    Dict("passed" => median_time < 10.0 && allocations < 1_000_000_000,
         "scope" => context.data["scope"], "samples" => samples, "first_gradient_seconds" => cold,
         "objective_median_seconds" => median(times), "gradient_median_seconds" => median_time,
         "gradient_p95_seconds" => quantile(gradients, 0.95), "gradient_allocated_bytes" => allocations,
         "process_maxrss_bytes" => Sys.iswindows() ? nothing : Sys.maxrss(),
         "baseline_speedup" => nothing, "speedup_reason" => "No accepted same-contract timing baseline",
         "budget" => Dict("kind" => "synthetic_pilot_only", "gradient_median_seconds_max" => 10,
                          "gradient_allocated_bytes_max" => 1_000_000_000),
         "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads())
end

function closure_report(context; seed=20261111)
    truth = Float64.(context.kernel["truth_theta"])
    exact = predict(truth, context)
    noisefree = fit(context; observations=exact, starts=6, seed=seed)
    error = maximum(abs.(predict(noisefree["theta"], context) - exact))
    rng = MersenneTwister(seed)
    L = cholesky(Symmetric(context.covariance + context.nuisance * transpose(context.nuisance))).L
    statistic = Float64[]
    failures = 0
    for j in 1:12
        result = fit(context; observations=exact + L * randn(rng, length(exact)), starts=2, seed=seed+j)
        failures += !result["converged"]
        push!(statistic, result["chi2"])
    end
    # A deliberately missing process offset must be detected, not absorbed as a validation success.
    stress = copy(exact)
    stress[1:8] .+= 20sqrt.(diag(context.covariance))[1:8]
    mismatch = fit(context; observations=stress, starts=3, seed=seed)["chi2"]
    dof = length(exact) - length(truth)
    Dict("passed" => noisefree["converged"] && error < 1e-7 && failures == 0 &&
                     0.4dof < mean(statistic) < 1.7dof && mismatch > 3dof,
         "scope" => context.data["scope"], "noisefree_max_prediction_error" => error,
         "pseudoexperiments" => 12, "convergence_failures" => failures,
         "mean_chi2" => mean(statistic), "dof" => dof, "misspecified_chi2" => mismatch,
         "coverage_qualification" => "Pilot goodness-of-fit check; parameter-interval coverage not yet certified for production")
end

function writejson(path, result)
    mkpath(dirname(abspath(path)))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        JSON.print(io, result, 2)
        write(io, '\n')
    end
    mv(temporary, path; force=true)
end

function run_job(job)
    BLAS.set_num_threads(1)
    context = load_context(job["data"], job["kernel"])
    action = job["action"]
    result = if action == "benchmark"
        benchmark_report(context)
    elseif action == "gradient"
        gradient_report(context; seed=job["seed"])
    elseif action == "performance"
        performance_report(context)
    elseif action == "closure"
        closure_report(context; seed=job["seed"])
    elseif action == "fit"
        fit(context; seed=job["seed"], starts=6, replicas=get(job, "replicas", 32))
    elseif action == "ablate"
        processes = unique(p["process"] for p in context.data["points"])
        Dict("scope" => context.data["scope"], "holdouts" => [Dict("held_out" => process,
             "result" => fit(context; seed=job["seed"], starts=3,
                 indices=findall(p -> p["process"] != process, context.data["points"]))) for process in processes])
    else
        error("unknown numerical job $action")
    end
    writejson(job["output"], result)
    get(result, "passed", get(result, "converged", true)) || error("numerical validation failed")
    result
end

end

