using CollinsFit, JSON, SHA, Dates
include("SIAJobBindings.jl")
include("SIASourceControls.jl")

function main(args=ARGS)
    length(args)==6 || error("usage: audit_sia_source_controls.jl FROZEN SOURCE_PLAN PLAN_SHA REPORT REPORT_SHA FRESH_AUDIT")
    frozen,plan,sha,input,inputsha,output=args
    context=SIAJobBindings.load(frozen,plan,sha)
    output=SIAJobBindings.fresh_output(output,plan)
    report=merge(SIASourceControls.flags(),Dict{String,Any}("status"=>"FAIL","scope"=>SIASourceControls.SCOPE,
        "theory_scope"=>"synthetic_closure","source_controls_success"=>false,"started_at"=>string(now(UTC))))
    SIAJobBindings.atomic_report(output,report)
    try
        SIAJobBindings.digest(input)==inputsha || error("changed externally bound source-control report")
        report=SIASourceControls.audit(context.request.data,context.parent,context.binding,JSON.parsefile(input))
        SIAJobBindings.recheck(context)
        SIAJobBindings.digest(input)==inputsha || error("source-control report bytes changed during audit")
        report["source_control_report_sha256"]=inputsha
    catch exception
        report["status"]="FAIL";report["source_controls_success"]=false
        report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));SIAJobBindings.atomic_report(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"fit_performed"=>false)))
    report["status"]==SIASourceControls.REVIEWED ? 0 : 2
end
if abspath(PROGRAM_FILE)==(@__FILE__);exit(main());end
