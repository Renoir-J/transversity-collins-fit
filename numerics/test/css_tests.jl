@testset "two-loop coupling and CSS domains/conventions" begin
    coupling = AsymptoticTwoLoopCoupling()
    model = KPSYCSS(;order=32)
    @test alpha_s(coupling,2.4) ≈ 0.327 atol=0.0005 rtol=0
    @test alpha_s(coupling,91.1876^2) ≈ 0.118 atol=0.0005 rtol=0
    @test iszero(evolution_time(coupling,2.4,2.4))
    @test evolution_time(coupling,100.0,2.4) ≈ -evolution_time(coupling,2.4,100.0)
    @test css_bstar(model,1e-7) ≈ 1e-7 rtol=1e-12
    @test 0 < css_bstar(model,1000.0) < model.bmax
    @test iszero(sudakov_perturbative(model,1.0,css_mu_b2(model,1.0)))
    @test sidis_snp(model,1.0,2.4,0.5) ≈ 0.274
    @test sidis_collins_snp(model,1.0,2.4,0.5,0.02) ≈ 0.194
    @test sia_snp(model,1.0,2.4,0.5,0.5) ≈ 0.336
    @test sia_collins_snp(model,1.0,2.4,0.5,0.5,0.02) ≈ 0.176
    @test ForwardDiff.derivative(g -> sidis_collins_snp(model,1.0,2.4,0.5,g),0.02) ≈ -4
    @test ForwardDiff.derivative(g -> sia_collins_snp(model,1.0,2.4,0.5,0.5,g),0.02) ≈ -8
    grid = collect(range(0.01,1;length=21))
    op = collinear_lo_operator(grid,100.0,2.4,coupling)
    @test op.coupling_model == :asymptotic_two_loop && op.lambda_GeV == 0.225
    @test op.alpha0 == alpha_s(coupling,2.4)
    hhat = collinear_lo_operator(grid,100.0,2.4,coupling;quantity=:KPSY_Hhat3_homogeneous)
    @test hhat.matrix == op.matrix
    @test_throws ErrorException AsymptoticTwoLoopCoupling(;nf=3.5)
    @test_throws ErrorException AsymptoticTwoLoopCoupling(;lambda_GeV=0)
    @test_throws ErrorException alpha_s(coupling,0.04)
    @test_throws ErrorException evolution_time(coupling,0.04,2.4)
    @test_throws ErrorException KPSYCSS(;bmax=5)
    @test_throws ErrorException KPSYCSS(;order=1)
    @test_throws ErrorException KPSYCSS(;g2=-0.1)
    @test_throws ErrorException KPSYCSS(;gh=0)
    @test_throws ErrorException sudakov_perturbative(model,0.0,100.0)
    @test_throws ErrorException sidis_snp(model,1.0,2.4,0.0)
    @test_throws ErrorException sidis_snp(model,1.0,2.4,1.1)
    @test_throws ErrorException sidis_collins_snp(model,1.0,2.4,0.5,model.gh)
    @test_throws ErrorException sia_collins_snp(model,1.0,2.4,0.5,0.5,NaN)
    @test_throws ErrorException ForwardDiff.derivative(l -> alpha_s(AsymptoticTwoLoopCoupling(;lambda_GeV=l),2.4),0.225)
    @test_throws ErrorException ForwardDiff.derivative(q -> sum(collinear_lo_operator(grid,q,2.4,coupling).matrix),100.0)
    @test_throws ErrorException ForwardDiff.derivative(g -> KPSYCSS(;g2=g).g2,0.84)
end
