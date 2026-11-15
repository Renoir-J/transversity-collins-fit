using JSON, SHA
include("JointLikelihoodAudit.jl")

function check_joint_cache(report_path,audit_path,accepted_path)
    result=JointLikelihoodAudit.audit_joint_report(JSON.parsefile(report_path))
    read(audit_path)==read(accepted_path) || error("changed accepted joint audit")
    accepted=JSON.parsefile(accepted_path)
    accepted["input_sha256"]==bytes2hex(sha256(read(report_path))) || error("changed joint report input")
    expected=Dict(name=>bytes2hex(sha256(read(joinpath(@__DIR__,name)))) for name in
        ("audit_joint_likelihood.jl","JointLikelihoodAudit.jl"))
    accepted["auditor_sources"]==expected || error("stale joint report auditor")
    for (key,value) in result
        haskey(accepted,key) && accepted[key]==value || error("stale or changed joint audit field: $key")
    end
    println(JSON.json(Dict("status"=>"PASS","scope"=>"synthetic_closure_cached_report_revalidation",
        "kind"=>result["kind"],"production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false)))
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)==3 || error("usage: check_joint_cache.jl REPORT AUDIT ACCEPTED")
    check_joint_cache(ARGS...)
end
