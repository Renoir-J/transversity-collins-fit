include("SIDISRecoveryBindings.jl")
const RB=SIDISRecoveryBindings
length(ARGS)==8 || error("FROZEN PREDECESSOR MANIFEST MANIFEST_SHA PHASE INDEX REPORT FRESH_AUDIT")
const RM=RB.load_manifest(ARGS[1:4]...)
include(joinpath(RM.frozen,"numerics/bin/audit_historical_sidis_job.jl"))
include("SIDISRecoveryDispatch.jl")

function recovery_main()
    phase=ARGS[5];index=parse(Int,ARGS[6]);reportfile=ARGS[7];output=ARGS[8]
    RB.require_new_output(RM,output)
    record=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "phase"=>phase,"recovery"=>RB.provenance(RM),SA.false_flags()...)
    try
        ctx=SA.bindings(RM.data["request"],RM.data["native"],RM.data["parton"])
        if phase in SA.POINT_PHASES
            reportfile==SA.point_report(RM.run,phase,index) || error("wrong actual point audit path")
            packet=SA.audit_point(ctx,RM.run,phase,index)
            record["point_index"]=index;record["point_id"]=packet["point"]["id"]
            record["saved_reference_steps_reviewed"]=length(get(packet,"reference_steps",[]))
        else
            index==0 && phase in SA.GLOBAL_PHASES && reportfile==SA.global_report(RM.run,phase) || error("wrong fresh global audit path")
            packet=JSON.parsefile(reportfile)
            get(packet,"recovery",nothing)==RB.provenance(RM) || error("unbound current global producer")
            record["review"]=SA.audit_global(ctx,RM.run,phase,packet)
            if phase=="diagnostic"
                record["scope"]=ctx.contract["scope"]
                record["accepted_prediction_diagnostic"]=record["all132_coverage_verified"]=true
            end
        end
        SA.stable(ctx);RB.stable(RM)
        record["source_fingerprint"]=ctx.fp;record["request_sha256"]=ctx.request.sha256
        record["policy_sha256"]=SA.digest(CollinsFit.HISTORICAL_SIDIS_JOB_POLICY)
        record["report_sha256"]=SA.digest(reportfile);record["status"]=SA.AUDITED
        record["unreexecuted_producer_measurements"]=["native evolution and rows","physical reference assembly","AD","phase timings","RSS"]
        record["primary_raw_tables_reparsed"]=false
    catch exception
        record["status"]="FAIL";record["error"]=sprint(showerror,exception,catch_backtrace())
        for flag in SA.FLAGS;record[flag]=false;end
    finally
        record["completed_at"]=string(now(UTC));RB.write_record(output,record)
    end
    println(JSON.json(Dict("status"=>record["status"],"phase"=>phase,"error"=>get(record,"error",nothing))));flush(stdout)
    record["status"]==SA.AUDITED || exit(2)
end
recovery_main()
