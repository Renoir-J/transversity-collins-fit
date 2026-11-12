@testset "matching conventions, generic values and invalid inputs" begin
    for n in (4,8,16)
        x,w = CollinsFit.gauss_rule(n)
        @test maximum(abs(sum(w.*x.^k)-(iseven(k) ? 2/(k+1) : 0.)) for k in 0:2n-1) < 2e-15
        @test x == -reverse(x)
        @test w == reverse(w)
    end
    nodes,weights = CollinsFit.gauss_rule(8)
    nodes[1]=weights[1]=NaN
    fresh = CollinsFit.gauss_rule(8)
    @test all(isfinite,fresh[1]) && all(isfinite,fresh[2])
    rule = matching_rule(:SIDIS,:transversity,.2,.3)
    @test length(rule.coordinates)==1
    @test match_leg(rule,[2.]) ≈ 2*(1-8*.3/(3pi))
    @test_throws ErrorException match_leg(rule,[2.];gluon=[0.])
    tree = matching_rule(:SIA,:ff,.2,0.)
    @test match_leg(tree,[3.])==3.
    rule = matching_rule(:SIDIS,:pdf,.2,.3)
    @test_throws ErrorException match_leg(rule,ones(length(rule.coordinates)))
    @test_throws ErrorException match_leg(rule,ones(3);gluon=zeros(3))
    @test_throws ErrorException match_leg(rule,fill(NaN,length(rule.coordinates));gluon=zeros(length(rule.coordinates)))
    # Analytic constant-input PDF convolution, independent of the rule weights.
    x,a = .2,.3
    base = 1-(8/3)*a/pi + a/pi*(2/3)*(-log(x)-1+x)
    value(t) = match_leg(rule,fill(t,length(rule.coordinates));gluon=zeros(length(rule.coordinates)))
    @test value(2.) ≈ 2base rtol=1e-12
    @test ForwardDiff.derivative(value,2.) ≈ base rtol=1e-12
    for (process,channel,x,alpha) in ((:SIA,:pdf,.2,.3),(:SIDIS,:bad,.2,.3),
        (:SIDIS,:ff,0.,.3),(:SIDIS,:ff,.2,-.1))
        @test_throws ErrorException matching_rule(process,channel,x,alpha)
    end
    @test_throws ErrorException matching_rule(:SIDIS,:ff,.2,.3;knots=[.4,.3])
    @test_throws ErrorException matching_rule(:SIDIS,:ff,.2,.3;order=true)
    @test_throws ErrorException ForwardDiff.derivative(a -> sum(matching_rule(:SIDIS,:ff,.2,a).quark_weights),.3)
end

@testset "Fourier signs, charge conjugation, shapes and width stability" begin
    radial = RadialRule(collect(0.:.5:16.);order=16)
    a,P,z1,z2 = .7,.3,.4,.6
    amplitude = exp.(-a.*radial.b.^2)
    k = P/z1
    gaussian(n) = k^n*exp(-k^2/(4a))/(2^(n+1)*a^(n+1)*2pi)
    sidis = sidis_fourier(radial,P,z1,amplitude,amplitude)
    sia = sia_fourier(radial,P,z1,z2,amplitude,amplitude)
    @test sidis.UU ≈ gaussian(0)/z1^2 rtol=1e-12
    @test sidis.Collins ≈ -gaussian(1)/(2z1^3) rtol=1e-12
    @test sia.UU ≈ gaussian(0)/z1^2 rtol=1e-12
    @test sia.Collins ≈ gaussian(2)/(4z1^3*z2) rtol=1e-12
    @test iszero(sidis_fourier(radial,0.,z1,amplitude,amplitude).Collins)
    @test iszero(sia_fourier(radial,0.,z1,z2,amplitude,amplitude).Collins)
    @test_throws ErrorException RadialRule([.01,1.])
    @test_throws ErrorException RadialRule([0.,1.,.5])
    @test_throws ErrorException sidis_fourier(radial,-.1,z1,amplitude,amplitude)
    @test_throws ErrorException sia_fourier(radial,.1,0.,z2,amplitude,amplitude)
    model = KPSYCSS()
    s = leptonic_plan(:SIDIS,model,radial,2.4,P,z1;flavors=[2,-2])
    e = leptonic_plan(:SIA,model,radial,110.,P,z1;z2=z2,flavors=[2,-2])
    one = ones(2,length(radial.b))
    first = vcat(one[1:1,:],3one[1:1,:])
    second = vcat(5one[1:1,:],7one[1:1,:])
    @test leptonic_unpolarized(s,first,second)/leptonic_unpolarized(s,one,one) ≈ 13.
    @test leptonic_unpolarized(e,first,second)/leptonic_unpolarized(e,one,one) ≈ 11.
    @test_throws ErrorException leptonic_plan(:SIA,model,radial,110.,P,z1)
    @test_throws ErrorException leptonic_plan(:SIDIS,model,radial,2.4,P,z1;flavors=[2])
    @test_throws ErrorException leptonic_plan(:SIDIS,model,radial,2.4,P,z1;flavors=[2,-2,2])
    @test_throws ErrorException leptonic_plan(:SIDIS,model,radial,2.4,P,z1;flavors=[21])
    @test_throws ErrorException leptonic_structure(s,ones(1,2),one,one,one;gc=.01)
    @test_throws ErrorException leptonic_collins(s,one,one,model.gh)
    narrow = leptonic_plan(:SIDIS,model,radial,2.4,.1,.05;flavors=[2,-2])
    @test isfinite(leptonic_collins(narrow,one,one,model.gh-1e-10))
    @test isfinite(ForwardDiff.derivative(g -> leptonic_collins(s,one,one,g),.01))
end
