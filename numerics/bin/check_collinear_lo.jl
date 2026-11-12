# Independent reconstructed-LO checks; never a production admission receipt.
using CollinsFit, ForwardDiff, LinearAlgebra, JSON, SHA, Dates
include("CollinearReference.jl")
using .CollinearReference

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const CONTRACT = joinpath(ROOT, "config/evolution-contract.json")
const SOURCES = [@__FILE__, CONTRACT, joinpath(ROOT,"numerics/src/CollinearLO.jl"),
    joinpath(ROOT,"numerics/src/CollinsFit.jl"), joinpath(ROOT,"numerics/Project.toml"),
    joinpath(@__DIR__,"CollinearReference.jl"), joinpath(ROOT,"numerics/src/RunningCoupling.jl"),
    joinpath(ROOT,"numerics/src/CSS.jl"),
    joinpath(ROOT,"numerics/Manifest.toml")]

source_hashes() = Dict(relpath(path, ROOT)=>bytes2hex(sha256(read(path))) for path in SOURCES)

function check_collinear_lo()
    sources = source_hashes()
    contract = JSON.parsefile(CONTRACT)
    policy = contract["validation_policy"]
    checks = Dict{String,Any}[]
    check(name, error, tolerance) = push!(checks, Dict("name"=>name, "error"=>error,
        "tolerance"=>tolerance, "passed"=>isfinite(error) && error <= tolerance))
    xmin = policy["xmin"]
    z = -1000.0 + 1000im
    asymptotic_psi = log(z) - 1/(2z) - 1/(12z^2) + 1/(120z^4)
    check("complex-contour digamma tail", abs(contour_digamma(z)-asymptotic_psi),
          policy["mellin_contour_refinement_max_abs"])
    Q02, Q2, alpha0, nf = 4.0, 100.0, 0.3, 3
    locations = [0.005, 0.02, 0.08, 0.2, 0.5, 0.8, 0.95]
    reference = [mellin_polynomial(x, Q2, Q02, alpha0; nf=nf) for x in locations]
    other_contour = [mellin_polynomial(x, Q2, Q02, alpha0; nf=nf, angle=2pi/3) for x in locations]
    check("Mellin contour refinement", maximum(abs.(reference-other_contour)),
          policy["mellin_contour_refinement_max_abs"])
    check("Mellin inversion at initial scale", maximum(abs.([
        mellin_polynomial(x, Q02, Q02, alpha0; nf=nf) - x*(1-x)^3 for x in locations])),
        policy["mellin_contour_refinement_max_abs"])
    grids = []
    previous = nothing
    for n in policy["grids"]
        grid = validation_grid(n, xmin)
        operator = nothing
        elapsed = @elapsed operator = collinear_lo_operator(grid, Q2, Q02, alpha0; nf=nf)
        initial = grid .* (1 .- grid).^3
        evolved = apply_operator(operator, initial)
        moments = [numerical_moment(grid, evolved, k) for k in 1:3]
        expected = [evolve_moment(6/prod(k+j for j in 1:4), k, Q2, Q02, alpha0; nf=nf) for k in 1:3]
        moment_error = maximum(abs.(moments ./ expected .- 1))
        values = [interpolate_grid(grid, evolved, x) for x in locations]
        mellin_error = max_scaled(values, reference)
        affine_exact = [(4/3) * (2x*log(x) + (1.5+2log1p(-x))*(1-x)) for x in grid[1:end-1]]
        affine_actual = collinear_lo_generator(grid) * (1 .- grid)
        check("affine generator n=$n", maximum(abs.(affine_actual[1:end-1] - affine_exact)),
              policy["affine_generator_max_abs"])
        entry = Dict("nodes"=>n, "moment_relative_error"=>moment_error,
            "mellin_scaled_error"=>mellin_error, "moments"=>moments,
            "reference_moments"=>expected, "point_predictions"=>values,
            "operator_build_seconds_informational"=>elapsed)
        if !isnothing(previous)
            check("moment refinement to n=$n", moment_error/previous["moment_relative_error"],
                  policy["refinement_error_ratio_max"])
            check("Mellin refinement to n=$n", mellin_error/previous["mellin_scaled_error"],
                  policy["refinement_error_ratio_max"])
        end
        push!(grids, entry)
        previous = entry
    end
    check("fine moments N=1,2,3", previous["moment_relative_error"], policy["fine_moment_max_relative"])
    check("fine x-space versus independent Mellin inversion", previous["mellin_scaled_error"],
          policy["fine_mellin_max_scaled"])

    grid = validation_grid(81, xmin)
    op = collinear_lo_operator(grid, Q2, Q02, alpha0; nf=nf)
    identity = collinear_lo_operator(grid, Q02, Q02, alpha0; nf=nf)
    middle = collinear_lo_operator(grid, 16.0, Q02, alpha0; nf=nf)
    final = collinear_lo_operator(grid, Q2, 16.0, alpha_lo(16.0, Q02, alpha0; nf=nf); nf=nf)
    check("identity evolution", maximum(abs.(identity.matrix-I)), policy["identity_semigroup_max_abs"])
    check("scale composition", maximum(abs.(final.matrix*middle.matrix-op.matrix)),
          policy["identity_semigroup_max_abs"])
    reverse = collinear_lo_operator(grid, Q02, Q2, alpha_lo(Q2, Q02, alpha0; nf=nf); nf=nf)
    check("forward/backward composition", maximum(abs.(reverse.matrix*op.matrix-I)),
          policy["identity_semigroup_max_abs"])
    hhat = collinear_lo_operator(grid, Q2, Q02, alpha0; nf=nf, quantity=:KPSY_Hhat3_homogeneous)
    check("explicit KPSY Hhat3 homogeneous kernel", maximum(abs.(hhat.matrix-op.matrix)), 0.0)

    basis = hcat(grid .* (1 .- grid).^3, grid.^2 .* (1 .- grid).^4, grid .* (1 .- grid).^5)
    initial(t) = basis * [t[1], t[2]^2, exp(t[3])]
    prediction(t) = apply_operator(op, initial(t))
    theta = [0.8, -0.3, 0.1]
    exact = op.matrix * basis * Diagonal([1.0, 2theta[2], exp(theta[3])])
    J = ForwardDiff.jacobian(prediction, theta)
    check("initial-parameter Jacobian", max_scaled(J, exact), policy["derivative_max_scaled"])
    errors = []
    for step in (1e-3, 1e-4, 1e-5, 1e-6)
        numerical = hcat([(prediction(theta + step*I(3)[:,j]) - prediction(theta - step*I(3)[:,j])) /
                          (2step) for j in 1:3]...)
        push!(errors, Dict("step"=>step, "error"=>max_scaled(numerical, J)))
    end
    plateau = any(errors[i]["error"] <= policy["derivative_max_scaled"] &&
                  errors[i+1]["error"] <= policy["derivative_max_scaled"] for i in 1:length(errors)-1)
    check("two adjacent finite-difference steps", plateau ? 0.0 : 1.0, 0.0)
    weights, direction = sin.(grid .* 3), [0.3, -0.5, 0.2]
    tangent = basis * Diagonal([1.0, 2theta[2], exp(theta[3])]) * direction
    check("operator adjoint on endpoint-zero tangent", abs(dot(weights, op.matrix*tangent) -
          dot(operator_adjoint(op, weights), tangent)), policy["identity_semigroup_max_abs"])
    scalar(t) = dot(weights, prediction(t))
    exact_hvp = Diagonal([0.0, 2dot(weights, op.matrix*basis[:,2]),
                         exp(theta[3])*dot(weights, op.matrix*basis[:,3])]) * direction
    check("parameter HVP", max_scaled(ForwardDiff.hessian(scalar, theta)*direction, exact_hvp),
          policy["derivative_max_scaled"])
    check("source hashes stable during job", sources == source_hashes() ? 0.0 : 1.0, 0.0)

    Dict("scope"=>"reconstructed_theory", "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL",
        "production_ready"=>false, "original_KPSY_reproduction"=>false,
        "contract_id"=>contract["id"], "contract_sha256"=>sources[relpath(CONTRACT, ROOT)],
        "julia_version"=>string(VERSION), "completed_at"=>string(now(UTC)),
        "sources"=>sources,
        "input"=>"analytic polynomial h0=x*(1-x)^3; no experimental input",
        "coupling"=>Dict("order"=>"one_loop", "nf"=>nf, "Q02"=>Q02, "Q2"=>Q2, "alpha0"=>alpha0),
        "locations"=>locations, "mellin_reference"=>reference, "grids"=>grids,
        "finite_difference_steps"=>errors, "checks"=>checks)
end

output = isempty(ARGS) ? joinpath(ROOT, "runs/collinear-lo/validation.json") : ARGS[1]
try
    report = check_collinear_lo()
    writejson(output, report)
    println(JSON.json(Dict("scope"=>report["scope"], "status"=>report["status"],
        "checks"=>length(report["checks"]), "failed"=>[c for c in report["checks"] if !c["passed"]])))
    report["status"] == "PASS" || exit(2)
catch exception
    writejson(output, Dict("scope"=>"reconstructed_theory", "status"=>"FAIL", "production_ready"=>false,
                          "error"=>sprint(showerror, exception)))
    rethrow()
end
