using Test, CollinsFit, ForwardDiff, LinearAlgebra, JSON
include("css_tests.jl")
include("parton_tests.jl")
include("leptonic_tests.jl")
include("initial_tests.jl")
include("polarized_tests.jl")
include("observable_tests.jl")
include("sidis_convention_tests.jl")
include("observable_reference_tests.jl")
include("leptonic_bin_tests.jl")
include("cache_snapshot_tests.jl")
include("job_scope_tests.jl")

@testset "LO transversity moments and fixed-operator AD" begin
    @test transversity_moment(1) ≈ -2/3
    @test evolve_moment(0.9, 1, 4.0, 4.0, 0.3) ≈ 0.9
    @test evolve_moment(0.9, 1, 100.0, 4.0, 0.3) < 0.9
    @test_throws ErrorException transversity_moment(0)
    U = [0.8 0.1; 0.2 0.7]
    f = t -> apply_operator(U, [t[1]^2, exp(t[2])])
    theta = [0.4, -0.1]
    @test ForwardDiff.jacobian(f, theta) ≈ U * Diagonal([0.8, exp(-0.1)])
    @test dot([0.3, 0.7], f(theta)) ≈ dot(operator_adjoint(U, [0.3, 0.7]), [theta[1]^2, exp(theta[2])])
end

@testset "reconstructed LO collinear evolution domains and derivatives" begin
    grid = vcat(exp.(range(log(1e-5), log(0.05); length=21)), collect(range(0.1, 1; length=31)))
    op = collinear_lo_operator(grid, 100.0, 4.0, 0.3)
    initial = grid .* (1 .- grid).^3
    @test all(isfinite, apply_operator(op, initial))
    @test iszero(last(apply_operator(op, initial)))
    @test collinear_lo_operator(grid, 100.0, 4.0, 0.3; quantity=:KPSY_Hhat3_homogeneous).matrix == op.matrix
    affine = (4/3) .* (2 .* grid[1:end-1] .* log.(grid[1:end-1]) .+
                      (1.5 .+ 2 .* log1p.(-grid[1:end-1])) .* (1 .- grid[1:end-1]))
    @test maximum(abs.((collinear_lo_generator(grid) * (1 .- grid))[1:end-1] - affine)) < 1e-10
    basis = grid .* (1 .- grid).^2
    f(t) = apply_operator(op, t[1]^2 .* basis)
    @test ForwardDiff.jacobian(f, [0.4])[:,1] ≈ op.matrix * (0.8 .* basis)
    weights = cos.(grid)
    @test dot(weights, op.matrix*basis) ≈ dot(operator_adjoint(op, weights), basis)
    @test iszero(last(operator_adjoint(op, weights)))
    @test evolve_moment(1.0, 1, 100.0, 4.0, 0.3; nf=5) ≈ (alpha_lo(100.0, 4.0, 0.3; nf=5)/0.3)^(4/23)
    for invalid in ([0.0, 0.2, 1.0], [0.1, 0.3, 0.2, 1.0], [0.1, 0.2, 0.9], [0.1, NaN, 1.0])
        @test_throws ErrorException collinear_lo_generator(invalid)
    end
    @test_throws ErrorException collinear_lo_operator(grid, 100.0, 4.0, 0.3; nf=3.5)
    @test_throws ErrorException collinear_lo_operator(grid, Inf, 4.0, 0.3)
    @test_throws ErrorException ForwardDiff.derivative(a -> sum(collinear_lo_operator(grid, 100.0, 4.0, a).matrix), 0.3)
    @test_throws ErrorException ForwardDiff.derivative(a -> sum(collinear_lo_generator(grid .* a)), 1.0)
    @test_throws ErrorException collinear_lo_operator(grid, 1e-30, 4.0, 0.3)
    for quantity in (:Collins, :Trento_H1perp1, :x_times_h1)
        @test_throws ErrorException collinear_lo_operator(grid, 100.0, 4.0, 0.3; quantity=quantity)
    end
    @test_throws ErrorException apply_operator(op, ones(length(grid)))
    @test_throws ErrorException apply_operator(op, initial[1:end-1])
    @test_throws ErrorException operator_adjoint(op, [NaN for _ in grid])
end

@testset "experimental intake cannot enter a numerical fit" begin
    root = normpath(joinpath(@__DIR__, "../.."))
    data = JSON.parsefile(joinpath(root, "fixtures/closure/data.json"))
    kernel = JSON.parsefile(joinpath(root, "fixtures/closure/kernel.json"))
    for scope in ("experimental_candidate", "historical_benchmark_candidate")
        data["scope"] = scope
        kernel["scope"] = scope
        mktempdir() do directory
            datafile, kernelfile = joinpath(directory, "data.json"), joinpath(directory, "kernel.json")
            CollinsFit.writejson(datafile, data)
            CollinsFit.writejson(kernelfile, kernel)
            @test_throws ErrorException load_context(datafile, kernelfile)
        end
    end
end

@testset "synthetic optimizer and prediction recovery" begin
    root = normpath(joinpath(@__DIR__, "../.."))
    context = load_context(joinpath(root, "fixtures/closure/data.json"), joinpath(root, "fixtures/closure/kernel.json"))
    truth = Float64.(context.kernel["truth_theta"])
    predictions = predict(truth, context)
    result = fit(context; observations=predictions, starts=3)
    @test result["converged"]
    @test maximum(abs.(predict(result["theta"], context) - predictions)) < 1e-7
end

@testset "profile likelihood agrees with covariance marginal quadratic" begin
    C = [0.04 0.01; 0.01 0.09]
    B = reshape([0.1, -0.2], 2, 1)
    residual = [0.3, -0.1]
    chi2, nuisance = profile_nuisance(residual, C, B)
    @test chi2 ≈ dot(residual, (C + B * B') \ residual)
    @test length(nuisance) == 1
end

