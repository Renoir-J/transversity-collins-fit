using Test, CollinsFit, ForwardDiff
if !isdefined(@__MODULE__,:JointLikelihoodReference)
    include("../bin/JointLikelihoodReference.jl")
end

@testset "direct physical optimizer settings retain complete boundaries" begin
    settings=BoundedJointSettings(source="synthetic unit settings")
    bounds=kpsy_parameter_bounds(.042)
    @test bounds.lower[1:2]==(-1.,-1.) && bounds.upper[1:2]==(1.,1.)
    @test all(i->bounds.lower[i]==0 && bounds.upper[i]==Inf,(3,4,5,6,9,10,11,12))
    @test bounds.lower[7:8]==(-Inf,-Inf) && bounds.upper[7:8]==(Inf,Inf)
    @test bounds.lower[13]==-Inf && bounds.upper[13]==prevfloat(.042)
    metadata=CollinsFit.bounded_settings_metadata(settings,.042)
    @test metadata["upper_bounds"][7]=="Inf" && metadata["lower_bounds"][13]=="-Inf"
    @test settings.maxeval==20000 && settings.trace
    for budget in (0,-1,true,1.5,typemax(Int))
        @test_throws ErrorException BoundedJointSettings(maxeval=budget,source="bad budget")
    end
    @test_throws ErrorException BoundedJointSettings(initial_step=zeros(13),source="zero steps")
    @test_throws ErrorException BoundedJointSettings(initial_step=ones(12),source="short steps")
    @test_throws ErrorException BoundedJointSettings(ftol_abs=NaN,source="bad tolerance")
    @test_throws ErrorException BoundedJointSettings(maxtime=-1.,source="bad time")
    @test_throws ErrorException BoundedJointSettings(trace=1,source="bad trace")
    @test_throws ErrorException BoundedJointSettings(source=" ")
    @test_throws ErrorException ForwardDiff.derivative(x->BoundedJointSettings(ftol_abs=x,source="no fitted setting").ftol_abs,1.)
    external=fill(.1,13);owned=BoundedJointSettings(initial_step=external,source="owned settings")
    external[1]=99.
    @test owned.initial_step[1]==.1
end

@testset "independent thirteen-coordinate interior and closed-boundary quadratics" begin
    settings=BoundedJointSettings(source="synthetic independent quadratic",trace=true)
    scales=[1.,1.,1.,1.,3.,3.,1.,1.,1.,1.,2.,2.,.02]
    start=[.2,-.3,.6,.8,2.,2.5,.1,-.2,.7,.8,1.1,1.2,.01]
    truth=[.4,-.6,.8,.9,2.8,3.4,.25,-.4,.7,.9,1.5,2.1,.02]
    for (target,expected) in ((truth,truth),
            ([1.4,-1.3,-.4,.6,1.2,-.5,8.,-6.,0.,1.1,1.5,-.2,.018],
             [1.,-1.,0.,.6,1.2,0.,8.,-6.,0.,1.1,1.5,0.,.018]))
        objective=t->sum(abs2,(t-target)./scales)
        result=CollinsFit.bounded_joint_start(objective,start,settings,.042)
        @test result["status"]=="PASS"
        @test result["solver_stop"] in ("FTOL_REACHED","XTOL_REACHED")
        @test maximum(abs.(result["theta"]-expected)./scales)<=1e-5
        @test abs(result["chi2"]-objective(expected))<=1e-5
        @test result["evaluation_accounting_matches"] && result["postcheck_passed"]
        @test length(result["trace"])==result["successful_evaluations"]==result["callback_attempts"]
        @test all(row->all(isfinite,KPSYParameters(row["theta"]).values),result["trace"])
        @test result["postcheck_evaluations"]==1
    end
end

@testset "resource exhaustion and callback failures never masquerade as convergence" begin
    start=[.2,-.3,.6,.8,2.,2.5,.1,-.2,.7,.8,1.1,1.2,.01]
    short=BoundedJointSettings(maxeval=1,source="deliberate resource stop")
    result=CollinsFit.bounded_joint_start(t->sum(abs2,t),start,short,.042)
    @test result["status"]=="INCOMPLETE" && result["solver_stop"]=="MAXEVAL_REACHED"
    @test result["optimizer_tolerance_stop"]===false
    @test result["solver_evaluations"]==1 && result["callback_attempts"]==1
    @test result["theta"]==start && result["postcheck_passed"]
    timed=BoundedJointSettings(maxtime=1e-9,source="deliberate time stop")
    stopped=CollinsFit.bounded_joint_start(t->sum(abs2,t),start,timed,.042)
    @test stopped["status"]=="INCOMPLETE" && stopped["solver_stop"]=="MAXTIME_REACHED"
    @test stopped["optimizer_tolerance_stop"]===false
    failed=CollinsFit.bounded_joint_start(t->error("explicit synthetic callback failure"),start,short,.042)
    @test failed["status"]=="FAIL" && failed["successful_evaluations"]==0
    @test occursin("explicit synthetic callback failure",failed["error"])
    nonfinite=CollinsFit.bounded_joint_start(t->Inf,start,short,.042)
    @test nonfinite["status"]=="FAIL" && nonfinite["retained_best_evaluation"]===nothing
    mutation=CollinsFit.bounded_joint_start(t->(t[1]=.3;0.),start,short,.042)
    @test mutation["status"]=="FAIL" && occursin("mutated parameter",mutation["error"])
end

@testset "complete joint job retains every start and its synthetic scope" begin
    fixture=JointLikelihoodReference.fixture_likelihood();l=fixture.likelihood;theta=fixture.theta
    short=BoundedJointSettings(maxeval=1,trace=false,source="synthetic full-job resource test")
    starts=[copy(theta),copy(theta)];starts[2][1]+=.01
    result=fit_joint(l,starts;settings=short)
    @test result["status"]=="INCOMPLETE" && length(result["starts"])==2
    @test all(r->r["status"]=="INCOMPLETE" && isempty(r["trace"]),result["starts"])
    @test result["scope"]=="synthetic_closure" && result["theory_scope"]=="reconstructed_theory"
    @test result["parameter_count"]==13 && result["point_count"]==8
    @test result["production_ready"]===false && result["physical_extraction"]===false && result["G2PERF_pass"]===false
    @test result["global_optimum_certified"]===false && result["interval_coverage_verified"]===false
    @test result["real_data_fit_implemented"]===false && result["sources_and_runtime_stable"]
    @test result["chi2"]==joint_objective(l,result["theta"])
    first_start=result["starts"][1]["start"][1];starts[1][1]=.99
    @test result["starts"][1]["start"][1]==first_start
    @test_throws ErrorException fit_joint(l,Vector{Float64}[];settings=short)
    bad=copy(theta);bad[1]=1.01
    @test_throws ErrorException fit_joint(l,[theta,bad];settings=short)
    damaged=deepcopy(l);damaged.observations[1]+=.001
    @test_throws ErrorException fit_joint(damaged,[theta];settings=short)
    boundary=copy(theta);boundary[1]=1.;boundary[2]=-1.;boundary[3]=0.;boundary[5]=0.
    bounded=fit_joint(l,[boundary];settings=short)
    @test bounded["status"]=="INCOMPLETE"
    @test bounded["starts"][1]["start"]==boundary && isfinite(bounded["chi2"])
end
