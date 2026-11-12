@testset "synthetic log-bilinear grid, weights, flavor/charge conversions" begin
    # Analytic synthetic functions, not copies of downloaded numerical arrays.
    xs,qs = [.05,.2,.6,1.],[1.,4.,16.]
    f(k,x,q) = k+2log(x)+3log(q)+.7log(x)*log(q)
    powers = CollinsFit.DSS_POWERS
    raw = [f(k,x,q)*x^a*(1-x)^b for (k,(a,b)) in enumerate(powers), q in qs, x in xs[1:end-1]]
    table = CollinsFit.legacy_table(:DSS14,xs,qs,raw,powers,.05)
    for (x,q) in ((.09,2.8),(.31,8.1))
        expected = [f(k,x,q)*x^a*(1-x)^b for (k,(a,b)) in enumerate(powers)]
        @test weighted_components(table,x,q) ≈ expected rtol=1e-13
    end
    @test iszero(sum(abs,weighted_components(table,1.,4.)))
    plus = CollinsFit.pion_components([10.,20.,30.,40.,50.,60.,2.,4.,6.],1)
    @test plus == [6.,4.,12.,8.,18.,12.,20.,25.,60.]
    minus = CollinsFit.pion_components([10.,20.,30.,40.,50.,60.,2.,4.,6.],-1)
    neutral = CollinsFit.pion_components([10.,20.,30.,40.,50.,60.,2.,4.,6.],0)
    @test minus == plus[[2,1,4,3,6,5,7,8,9]]
    @test neutral == (plus+minus)/2
    @test pion_ff(table,4,1,.2,4.) ≈ pion_ff(table,-4,1,.2,4.)
    @test pion_ff(table,5,1,.2,4.) ≈ pion_ff(table,-5,1,.2,4.)
    @test pion_ff(table,2,1,.2,4.)*.2 ≈ (raw[1,2,2]+raw[7,2,2])/2
    @test pion_ff(table,[2,-2,4,-4,21],1,.2,4.) ≈ [pion_ff(table,p,1,.2,4.) for p in (2,-2,4,-4,21)]
    @test_throws ErrorException pion_ff(table,[2,6],1,.2,4.)
    @test_throws ErrorException weighted_components(table,.049,4.)
    @test_throws ErrorException weighted_components(table,.2,.5)
    @test_throws ErrorException pion_ff(table,6,1,.2,4.)
    @test_throws ErrorException pion_ff(table,2,2,.2,4.)
    @test_throws ErrorException parton_helicity(table,2,.2,4.)
    @test_throws ErrorException ForwardDiff.derivative(x -> sum(weighted_components(table,x,4.)),.2)
    @test_throws ErrorException weighted_components(table,NaN,4.)
    @test_throws ErrorException weighted_components(table,.2,Inf)
    @test !parton_domain(table,.2,.5;policy=:reference)["within_grid"]
    @test !parton_domain(table,.2,.5;policy=:reference)["production_ready"]
    @test all(isfinite,weighted_components(table,.2,.5;policy=:reference))
    @test CollinsFit.fortran_real("0.3") == Float64(Float32(.3))
    @test CollinsFit.fortran_real("0.3D0") == .3
    @test CollinsFit.fortran_real("1.0d-5") == 1e-5
    @test_throws ErrorException CollinsFit.fortran_real("1.0/2.0")
    mktempdir() do directory
        mkpath(joinpath(directory,"CT10nlo"))
        write(joinpath(directory,"CT10nlo/CT10nlo.info"),"explicitly synthetic wrong input")
        @test_throws ErrorException CollinsFit.verified_parton_asset(directory,"CT10_info")
    end
end
