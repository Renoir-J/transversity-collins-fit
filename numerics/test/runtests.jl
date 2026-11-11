using Test, CollinsFit, ForwardDiff, LinearAlgebra

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

