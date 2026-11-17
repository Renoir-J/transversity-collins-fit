include("SIDISContinuationPaths.jl")
const CP=SIDISContinuationPaths
length(ARGS)==7 || error("FROZEN_CHECKOUT MANIFEST MANIFEST_SHA PHASE INDEX REPORT FRESH_AUDIT")
const CM=CP.load_manifest(ARGS[1],ARGS[2],ARGS[3])
include(joinpath(CM.frozen,"numerics/bin/audit_historical_sidis_job.jl"))
include("SIDISContinuationDispatch.jl")

function continuation_main()
    phase=ARGS[4];index=parse(Int,ARGS[5]);reportfile=ARGS[6];output=ARGS[7]
    CP.require_new_output(CM,output)
    record=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "phase"=>phase,"continuation"=>CP.provenance(CM),SA.false_flags()...)
    try
        ctx=SA.bindings(CM.data["request"],CM.data["native"],CM.data["parton"])
        if phase in SA.POINT_PHASES
            reportfile==SA.point_report(CM.run,phase,index) || error("wrong point audit path")
            packet=SA.audit_point(ctx,CM.run,phase,index)
            record["point_index"]=index;record["point_id"]=packet["point"]["id"]
            record["saved_reference_steps_reviewed"]=length(get(packet,"reference_steps",[]))
        else
            index==0 && phase in SA.GLOBAL_PHASES && reportfile==SA.global_report(CM.run,phase) || error("wrong global audit path")
            packet=JSON.parsefile(reportfile)
            record["review"]=SA.audit_global(ctx,CM.run,phase,packet)
            if phase=="diagnostic"
                record["scope"]=ctx.contract["scope"]
                record["accepted_prediction_diagnostic"]=record["all132_coverage_verified"]=true
            end
        end
        if !(phase=="native123" && index<=13)
            packet["continuation"]==CP.provenance(CM) || error("unbound new producer")
        end
        SA.stable(ctx);CP.stable(CM)
        record["source_fingerprint"]=ctx.fp;record["request_sha256"]=ctx.request.sha256
        record["policy_sha256"]=SA.digest(CollinsFit.HISTORICAL_SIDIS_JOB_POLICY)
        record["report_sha256"]=SA.digest(reportfile);record["status"]=SA.AUDITED
        record["unreexecuted_producer_measurements"]=["native evolution and rows","physical reference assembly","AD","phase timings","RSS"]
        record["primary_raw_tables_reparsed"]=false
    catch exception
        record["status"]="FAIL";record["error"]=sprint(showerror,exception,catch_backtrace())
        for flag in SA.FLAGS;record[flag]=false;end
    finally
        record["completed_at"]=string(now(UTC));CP.write_record(output,record)
    end
    println(JSON.json(Dict("status"=>record["status"],"phase"=>phase,"error"=>get(record,"error",nothing))));flush(stdout)
    record["status"]==SA.AUDITED || exit(2)
end
continuation_main()
