using CollinsFit, JSON, SHA
include("BoundedJointAudit.jl")

function check_bounded_cache(report_path,audit_path,accepted_path,algorithm,budget)
    report=JSON.parsefile(report_path)
    report["algorithm"]==algorithm && report["maxeval_per_start"]==parse(Int,budget) || error("stale requested bounded configuration")
    result=BoundedJointAudit.audit_bounded_report(report)
    result["status"]=="PASS" && result["accepted_fitting_gate"]===true || error("failed or incomplete fitting report")
    read(audit_path)==read(accepted_path) || error("changed accepted bounded audit")
    accepted=JSON.parsefile(accepted_path)
    accepted["input_sha256"]==bytes2hex(sha256(read(report_path))) || error("changed bounded report input")
    accepted["auditor_sources"]==BoundedJointAudit.audit_source_hashes() || error("stale bounded auditor/workflow sources")
    for (key,value) in result
        haskey(accepted,key) && accepted[key]==value || error("stale or changed bounded audit field: $key")
    end
    println(JSON.json(Dict("status"=>"PASS","scope"=>"synthetic_closure_cached_fitting_revalidation",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false)))
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)==5 || error("usage: check_bounded_joint_cache.jl REPORT AUDIT ACCEPTED ALGORITHM MAXEVAL")
    check_bounded_cache(ARGS...)
end
