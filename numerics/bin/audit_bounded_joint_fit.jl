using CollinsFit, JSON, SHA, Dates
include("BoundedJointAudit.jl")

function main()
    length(ARGS)==2 || error("usage: audit_bounded_joint_fit.jl REPORT FRESH_OUTPUT")
    input,output=ARGS;ispath(output) && error("fresh audit path required")
    result=Dict{String,Any}("status"=>"FAIL","review_status"=>"FAIL",
        "scope"=>"synthetic_closure_retained_report_review","accepted_fitting_gate"=>false,
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false)
    try
        merge!(result,BoundedJointAudit.audit_bounded_report(JSON.parsefile(input)))
    catch exception
        result["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        result["input_sha256"]=bytes2hex(sha256(read(input)))
        result["auditor_sources"]=BoundedJointAudit.audit_source_hashes()
        result["completed_at"]=string(now(UTC));mkpath(dirname(output));write(output,JSON.json(result,2)*"\n")
    end
    println(JSON.json(result));result["status"]=="PASS" && result["accepted_fitting_gate"]===true || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
