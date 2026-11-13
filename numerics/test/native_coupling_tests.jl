include("../bin/NativeCouplingReference.jl")
using .NativeCouplingReference

@testset "synthetic native-coupling protocol and independent beta references" begin
    response="CF-COUPLING-1\n2 5 2\n0.327\n0.25\n"
    @test parse_native_response(response,2,5,2)==[.327,.25]
    for malformed in ("",replace(response,"CF-COUPLING-1"=>"OTHER"),replace(response,"2 5 2"=>"2 4 2"),
        "CF-COUPLING-1\n2 5 2\n0.327\n",replace(response,"0.25"=>"NaN"),response*"0.1\n")
        @test_throws Exception parse_native_response(malformed,2,5,2)
    end
    @test startswith(native_request(.327,2.4,[2.4,10.],5,2,.05),"CF-COUPLING-1\n2 5 2\n")
    @test_throws ErrorException native_request(.327,2.4,[1.],5,2,.05)
    @test_throws ErrorException native_request(.327,2.4,[10.],true,2,.05)
    @test_throws ErrorException native_request(.327,2.4,[10.],5,3,.05)
    @test_throws ErrorException native_request(.327,2.4,[10.],5,2,0.)
    @test_throws ErrorException native_request(NaN,2.4,[10.],5,2,.05)

    for nf in 3:6
        initial=rg_reference(.327,2.4,2.4,nf,2)
        @test initial.alpha==.327 && initial.tau==0.
        one=rg_reference(.327,1000.,2.4,nf,1)
        analytic=.327/(1+(11-2nf/3)*.327*log(1000/2.4)/(4pi))
        @test one.alpha ≈ analytic rtol=3e-16
        two=rg_reference(.327,1000.,2.4,nf,2)
        @test 0<two.alpha<one.alpha<.327
        @test two.implicit_residual<1e-60
        middle=rg_reference(.327,10.,2.4,nf,2)
        chain=rg_reference(middle.alpha,1000.,10.,nf,2)
        @test chain.alpha ≈ two.alpha rtol=4e-16
        @test middle.tau+chain.tau ≈ two.tau rtol=5e-16
    end
    @test_throws ErrorException ForwardDiff.derivative(a->rg_reference(a,10.,2.4,5,2).alpha,.327)
    rule=log_rule(2.4,1000.,32)
    # Constant alpha has an exact logarithmic evolution time; this catches Q/Q2 and 2pi factors.
    @test sum(rule.weights)*.3 ≈ .3*log(1000/2.4)/(2pi) rtol=3e-15
    @test isempty(log_rule(2.4,2.4,16).nodes)
    @test_throws ErrorException log_rule(10.,2.4,16)
    @test !comparison("changed result",[.2],[.3],1e-10)["passed"]
    @test !comparison("nonfinite result",[NaN],[.3],1e-10)["passed"]
end
