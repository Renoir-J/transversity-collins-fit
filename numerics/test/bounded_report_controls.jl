using Test, CollinsFit, JSON
include("../bin/BoundedJointAudit.jl")
length(ARGS)==1 || error("one complete bounded fitting report required")
report=JSON.parsefile(only(ARGS))

@testset "bounded fitting reader preserves independent recovery and scope" begin
    result=BoundedJointAudit.audit_bounded_report(report)
    @test result["status"]=="PASS" && result["accepted_fitting_gate"]===true
    @test result["scope"]=="synthetic_closure_retained_report_review" && result["production_ready"]===false
    @test result["physical_extraction"]===false && result["G2PERF_pass"]===false
    @test result["independent_local_rank"]==13 && length(result["recomputed_checks"])==11
    @test result["all_trace_parameter_domains_verified"] && !result["all_trace_objectives_replayed"]
    @test result["global_optimum_certified"]===false && result["interval_coverage_verified"]===false
    mutators=[
        ("false producer PASS",r->(r["status"]="INCOMPLETE"),"false producer status"),
        ("production promotion",r->(r["production_ready"]=true),"invalid admission flag"),
        ("child data promotion",r->(r["fit"]["experimental_data_admitted"]=true),"invalid child admission"),
        ("boundary promotion",r->(r["boundary_resource_control"]["production_ready"]=true),"invalid boundary admission"),
        ("scope promotion",r->(r["scope"]="production"),"changed scope"),
        ("changed precision",r->(r["fit"]["settings"]["ftol_abs"]*=100),"changed solver settings"),
        ("changed algorithm",r->(r["fit"]["settings"]["algorithm"]="unreviewed"),"changed solver settings"),
        ("missing source",r->pop!(r["sources"],"numerics/src/BoundedJointFit.jl"),"incomplete producer inventory"),
        ("rehashed source claim",r->(r["sources"]["numerics/src/BoundedJointFit.jl"]=repeat("0",64)),"stale source"),
        ("changed runtime",r->(r["fit"]["runtime"]["native_version"]="0"),"stale solver wrapper/library"),
        ("changed observations",r->(r["data"]["observations"][1]+=.0001),"changed synthetic data"),
        ("changed covariance",r->(r["data"]["covariance"][1][1]*=2),"changed synthetic data"),
        ("changed relative nuisance",r->(r["data"]["relative"][1][1]+=.1),"changed synthetic data"),
        ("changed point order",r->reverse!(r["data"]["point_ids"]),"changed data order"),
        ("changed start",r->(r["starts_requested"][1][1]+=.01),"changed truth/starts"),
        ("one-ULP changed initial callback",r->(r["fit"]["starts"][1]["trace"][1]["theta"][6]=nextfloat(r["fit"]["starts"][1]["trace"][1]["theta"][6])),"changed native initial evaluation"),
        ("out of bounds iterate",r->(r["fit"]["starts"][1]["trace"][2]["theta"][1]=1.01),"Nh must lie"),
        ("trace order",r->(r["fit"]["starts"][1]["trace"][2]["evaluation"]=1),"changed trace order"),
        ("missing trace",r->pop!(r["fit"]["starts"][1]["trace"]),"incomplete evaluation trace"),
        ("retained best changed",r->(r["fit"]["starts"][1]["retained_best_evaluation"]["chi2"]+=.25),"changed retained best"),
        ("counter mismatch",r->(r["fit"]["starts"][1]["successful_evaluations"]-=1),"incomplete evaluation trace"),
        ("budget stop hidden",r->(r["fit"]["starts"][1]["solver_stop"]="MAXEVAL_REACHED"),"false convergence"),
        ("selected theta changed",r->(r["fit"]["theta"][1]+=.01),"changed selected start"),
        ("selected chi2 changed",r->(r["fit"]["chi2"]+=.25),"changed numerical operand"),
        ("prediction changed",r->(r["fit"]["predictions"][1]+=.0001),"changed numerical operand"),
        ("pull changed",r->(r["fit"]["nuisance"][1]+=.1),"changed numerical operand"),
        ("parameter distance forged",r->(r["closure"]["scaled_parameter_distances_to_both_sign_branches"][1]=0.),"changed numerical operand"),
        ("sign branch relabelled",r->(r["closure"]["selected_sign_branch"]=3),"changed closure claim"),
        ("local rank forged",r->(r["local_information"]["rank"]=12),"changed local information"),
        ("start dropped",r->pop!(r["fit"]["starts"]),"missing requested start"),
        ("check bit forged",r->(r["checks"][1]["passed"]=false),"false check bit"),
        ("missing check",r->pop!(r["checks"]),"incomplete check ledger"),
        ("duplicate check",r->(r["checks"][2]=deepcopy(r["checks"][1])),"incomplete check ledger"),
        ("boundary stop promoted",r->(r["boundary_resource_control"]["starts"][1]["solver_stop"]="FTOL_REACHED"),"boundary resource stop promoted"),
        ("recovery threshold changed",r->(r["contract"]["validation"]["closure_max_scaled_parameter_error"]*=10),"changed contract"),
        ("prospective policy changed",r->(r["solver_comparison_contract"]["scientific_acceptance_thresholds_changed"]=true),"changed prospective solver comparison")
    ]
    for (name,mutate,expected_error) in mutators
        altered=deepcopy(report);mutate(altered)
        exception=try BoundedJointAudit.audit_bounded_report(altered);nothing catch error;error end
        @test exception isa ErrorException
        @test exception!==nothing && occursin(expected_error,sprint(showerror,exception))
        println(JSON.json(Dict("name"=>name,"rejected"=>exception!==nothing,"expected_error"=>expected_error)))
    end
end
