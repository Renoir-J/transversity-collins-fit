@testset "physical KPSY coordinates, exact boundaries and pion flavor relations" begin
    p = kpsy_published_parameters()
    @test p.values[2]==-1.
    @test p.values[11]==p.values[12]==0.
    x = p.values[4]/(p.values[4]+p.values[6])
    @test initial_transversity(p,1,x,2.) ≈ -2. rtol=1e-13
    @test initial_transversity(p,2,.1,2.)>0
    @test all(iszero(initial_transversity(p,q,.1,2.)) for q in (-2,-1,3,-3,4,-4,5,-5))
    z,ff = .3,[3.,1.,.5,.25]
    favored = p.values[7]*z^p.values[9]*ff[1]
    unfavored = p.values[8]*z^p.values[10]*ff[2]
    @test initial_collins(p,2,1,z,ff) ≈ favored
    @test initial_collins(p,-1,1,z,ff) ≈ favored
    @test initial_collins(p,-2,-1,z,ff) ≈ favored
    @test initial_collins(p,1,-1,z,ff) ≈ favored
    @test initial_collins(p,-2,1,z,ff) ≈ unfavored
    @test initial_collins(p,2,0,z,ff) ≈ (favored+unfavored)/2
    @test initial_collins(p,3,-1,z,ff) ≈ unfavored*.5
    @test initial_collins(p,-3,1,z,ff) ≈ unfavored*.25
    @test iszero(initial_collins(p,4,1,z,ff))
    @test trento_first_moment(favored,z,.14)>0
    theta = collect(p.values)
    @test all(iszero,ForwardDiff.gradient(t -> initial_collins(KPSYParameters(t),2,1,1.,zeros(4)),theta))
    @test all(iszero,ForwardDiff.hessian(t -> initial_collins(KPSYParameters(t),2,1,1.,zeros(4)),theta))
    @test all(iszero,ForwardDiff.gradient(t -> initial_transversity(KPSYParameters(t),1,1.,0.),theta))
    beta_derivative = ForwardDiff.gradient(t -> initial_collins(KPSYParameters(t),2,1,z,ff),theta)[11]
    @test beta_derivative ≈ favored*log1p(-z)
    # These are shape value limits; their normalization gradients are singular.
    theta[3]=0.;theta[5]=2.
    @test initial_transversity(KPSYParameters(theta),2,.3,1.) ≈ .85*.7^2
    @test_throws ErrorException ForwardDiff.gradient(t -> initial_transversity(KPSYParameters(t),2,.3,1.),theta)
    theta[5]=0.
    @test initial_transversity(KPSYParameters(theta),2,.3,1.) ≈ .85
    for (i,value) in ((1,1.01),(3,-.1),(11,-.1),(13,.042),(7,Inf))
        invalid=collect(p.values);invalid[i]=value
        @test_throws ErrorException KPSYParameters(invalid)
    end
    @test_throws ErrorException initial_transversity(p,2,.3,-1.)
    @test_throws ErrorException initial_collins(p,21,1,z,ff)
    @test_throws ErrorException initial_collins(p,2,2,z,ff)
    @test_throws ErrorException trento_first_moment(favored,z,0.)
    @test_throws ErrorException trento_first_moment(1+im,z,.14)
end

@testset "generic initial matrices and matched evolution settings" begin
    p = kpsy_published_parameters()
    x,z = collect(range(.01,1.;length=31)),collect(range(.05,1.;length=27))
    bound = [(f==1 ? 2. : 1.)*u*(1-u)^3 for f in 1:2,u in x]
    ff = [(1+.2f)*u*(1-u)^2 for f in 1:4,u in z]
    inputs = CollinsFit.initial_conditions(x,z,bound,ff)
    initial = initial_distributions(inputs,p)
    flipped=collect(p.values);flipped[[1,2,7,8]].*=-1
    mirrored=initial_distributions(inputs,KPSYParameters(flipped))
    @test mirrored.h1 == -initial.h1
    @test mirrored.Hhat3 == -initial.Hhat3
    @test size(initial.h1)==(10,31)
    @test size(initial.Hhat3)==(10,27)
    coupling = AsymptoticTwoLoopCoupling()
    hop = collinear_lo_operator(x,10.,2.4,coupling)
    Cop = collinear_lo_operator(z,10.,2.4,coupling;quantity=:KPSY_Hhat3_homogeneous)
    evolved = evolve_initial_conditions(inputs,p,hop,Cop)
    @test all(iszero,evolved.h1[:,end])
    @test all(iszero,evolved.Hhat3[:,end])
    @test evolved.h1[7,:] ≈ apply_operator(hop,initial.h1[7,:])
    jacobian = ForwardDiff.jacobian(t -> vec(evolve_initial_conditions(inputs,KPSYParameters(t),hop,Cop).Hhat3),collect(p.values))
    @test all(isfinite,jacobian)
    @test_throws ErrorException evolve_initial_conditions(inputs,p,hop,hop)
    @test_throws ErrorException evolve_initial_conditions(inputs,p,collinear_lo_operator(x,10.,2.4,.3;nf=5),Cop)
    @test_throws ErrorException apply_operator(hop,ones(length(x),2))
    @test_throws ErrorException apply_operator(hop,zeros(length(x)-1,2))
end
