using JSON, SHA, Dates
include("JointLikelihoodAudit.jl")

function main()
    length(ARGS)==2 || error("usage: audit_joint_likelihood.jl REPORT OUTPUT")
    input,output=ARGS;ispath(output) && error("use a fresh audit path")
    result=Dict{String,Any}("status"=>"FAIL","scope"=>"synthetic_closure_report_audit",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false)
    try
        merge!(result,JointLikelihoodAudit.audit_joint_report(JSON.parsefile(input)))
    catch exception
        result["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        result["input_sha256"]=bytes2hex(sha256(read(input)))
        result["auditor_sources"]=Dict(basename(path)=>bytes2hex(sha256(read(path))) for path in
            (@__FILE__,joinpath(@__DIR__,"JointLikelihoodAudit.jl")))
        result["completed_at"]=string(now(UTC));mkpath(dirname(output));write(output,JSON.json(result,2)*"\n")
    end
    println(JSON.json(result));result["status"]=="PASS" || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
