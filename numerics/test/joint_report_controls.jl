using Test, JSON
include("../bin/JointLikelihoodAudit.jl")
length(ARGS)==1 || error("usage: joint_report_controls.jl REPORT")
report=JSON.parsefile(only(ARGS))
@testset "joint evidence rejects altered reports" begin
    @test JointLikelihoodAudit.audit_joint_report(report)["status"]=="PASS"
    alterations=(
        r->(r["scope"]="validated_production"),
        r->(r["production_ready"]=true),
        r->(r["likelihood"]["physical_extraction"]=true),
        r->(r["status"]="FAIL"),
        r->(r["sources"]["numerics/src/JointLikelihood.jl"]="stale"),
        r->(r["input_contract_sha256"]="stale"),
        r->(r["policy"]["derivative_max_scaled"]=1.),
        r->reverse!(r["data"]["point_ids"]),
        r->(r["data"]["covariance"][1][1]+=.001),
        r->(r["data"]["nuisance_prior"][1][1]+=.1),
        r->(r["data"]["nuisance_ids"][1]="invented"),
        r->(r["profile"]["chi2"]+=.01),
        r->(r["derivatives"]["finite_difference_trials"][1]["gradient_function_values_decimal"][1][1]="99"),
        r->(r["derivatives"]["hvp"][1]+=.01),
        r->pop!(r["checks"]),
        r->push!(r["checks"],deepcopy(first(r["checks"]))),
        r->begin c=first(r["checks"]);c["actual"].+=1.;c["reference"].+=1.;c["error"]=0.;c["passed"]=true;end,
        r->(first(r["negative_controls"])["rejected"]=false))
    for alter in alterations
        changed=deepcopy(report);alter(changed)
        @test_throws ErrorException JointLikelihoodAudit.audit_joint_report(changed)
    end
end
