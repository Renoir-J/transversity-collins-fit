using CollinsFit, JSON, Dates
include("HistoricalSIDISJobRecords.jl")
const SA=HistoricalSIDISJobRecords

function main()
    length(ARGS)==8 || error("REQUEST PHASE INDEX_OR_ZERO NATIVE_BUILD PARTON_BUILD RUN REPORT FRESH_AUDIT")
    request,phase,index_text,native,inputs,run,reportfile,output=ARGS;index=parse(Int,index_text)
    !ispath(output) || error("preserve earlier independent audits")
    record=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "phase"=>phase,SA.false_flags()...)
    try
        ctx=SA.bindings(request,native,inputs)
        if phase in SA.POINT_PHASES
            reportfile==SA.point_report(abspath(run),phase,index) || error("wrong point audit path")
            packet=SA.audit_point(ctx,abspath(run),phase,index)
            record["point_index"]=index;record["point_id"]=packet["point"]["id"]
            record["saved_reference_steps_reviewed"]=length(get(packet,"reference_steps",[]))
        else
            index==0 && phase in SA.GLOBAL_PHASES && reportfile==SA.global_report(abspath(run),phase) || error("wrong global audit identity")
            packet=JSON.parsefile(reportfile)
            record["review"]=SA.audit_global(ctx,abspath(run),phase,packet)
            if phase=="diagnostic"
                record["scope"]=ctx.contract["scope"]
                record["accepted_prediction_diagnostic"]=record["all132_coverage_verified"]=true
            end
        end
        SA.stable(ctx)
        record["source_fingerprint"]=ctx.fp;record["request_sha256"]=ctx.request.sha256
        record["policy_sha256"]=SA.digest(CollinsFit.HISTORICAL_SIDIS_JOB_POLICY)
        record["report_sha256"]=SA.digest(reportfile);record["status"]=SA.AUDITED
        record["unreexecuted_producer_measurements"]=["native evolution and rows","physical reference assembly","AD","phase timings","RSS"]
        record["primary_raw_tables_reparsed"]=false
    catch exception
        record["status"]="FAIL";record["error"]=sprint(showerror,exception,catch_backtrace())
        for key in SA.FLAGS;record[key]=false;end
    finally
        record["completed_at"]=string(now(UTC));writejson(output,record)
    end
    println(JSON.json(Dict("status"=>record["status"],"phase"=>phase,"error"=>get(record,"error",nothing))));flush(stdout)
    record["status"]==SA.AUDITED || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
