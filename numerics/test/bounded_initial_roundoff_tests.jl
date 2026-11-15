using Test, CollinsFit
if !isdefined(@__MODULE__,:BoundedJointAudit)
    include("../bin/BoundedJointAudit.jl")
end
using .BoundedJointAudit: expected_first_evaluation
@testset "pinned native initial mapping uses exact basic binary64 operations" begin
    truth=[.5,-.7,.8,.9,2.8,3.4,.25,-.4,.7,.9,1.5,2.1,.02]
    first=[.42,-.62,.9,1.,2.6,3.1,.22,-.36,.85,1.05,1.35,1.9,.018]
    second=[.58,-.78,.7,.8,3.1,3.7,.29,-.44,.6,.8,1.7,2.3,.023]
    boundary=copy(truth);boundary[1]=1.;boundary[2]=-1.;boundary[3]=0.;boundary[5]=0.
    near=copy(truth);near[1]=.999;near[2]=-.999;near[3]=.01;near[5]=.001;near[13]=.041
    signed_zero=copy(truth);signed_zero[3]=-0.;signed_zero[5]=0.
    cases=[truth,first,second,boundary,near,signed_zero]
    for steps in ([.1,.1,.2,.2,.3,.3,.1,.1,.2,.2,.3,.3,.003],fill(.1,13)),start in cases
        settings=BoundedJointSettings(algorithm=:LN_BOBYQA,initial_step=steps,maxeval=1,source="exact native initial map control")
        result=CollinsFit.bounded_joint_start(t->sum(abs2,t),start,settings,.042)
        @test result["status"]=="INCOMPLETE" && result["solver_stop"]=="MAXEVAL_REACHED"
        @test length(result["trace"])==result["callback_attempts"]==result["solver_evaluations"]==1
        actual=only(result["trace"])["theta"]
        expected=expected_first_evaluation(:LN_BOBYQA,start,settings.initial_step,.042)
        @test reinterpret(UInt64,actual)==reinterpret(UInt64,expected)
        @test all(isfinite,KPSYParameters(actual).values)
    end
    settings=BoundedJointSettings(algorithm=:LN_BOBYQA,source="roundtrip proof")
    mapped=expected_first_evaluation(:LN_BOBYQA,first,settings.initial_step,.042)
    @test mapped[6]==prevfloat(first[6]) && mapped[6]!=first[6]
    @test expected_first_evaluation(:LN_SBPLX,first,settings.initial_step,.042)==first
    changed=copy(mapped);changed[6]=nextfloat(changed[6])
    @test reinterpret(UInt64,changed)!=reinterpret(UInt64,mapped)
end
