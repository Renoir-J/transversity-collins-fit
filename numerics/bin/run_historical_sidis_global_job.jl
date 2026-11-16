using CollinsFit, JSON, Dates
include("HistoricalSIDISJobRecords.jl")
const SG=HistoricalSIDISJobRecords

function populate_global!(report,ctx,run,phase)
    report["dependency_sha256"]=Dict(p=>SG.digest(p) for p in SG.global_dependencies(ctx,run,phase))
    report["source_ids"]=SG.all_ids(ctx)
    if phase=="source-math"
        SG.HistoricalSIDISSourceMath.populate!(report,ctx)
        return report
    end
    initial=SG.point_packets(ctx,run,"native123")
    first=sidis_job_native_decision(SG.all_ids(ctx),initial)
    fourth=phase!="native-select3" && first["selected_level"]==0 ? SG.point_packets(ctx,run,"native4") : nothing
    native=fourth===nothing ? first : sidis_job_native_decision(SG.all_ids(ctx),initial;fourth)
    if phase in ("native-select3","native-selection")
        merge!(report,phase=="native-select3" ? first : native)
        return report
    end
    native["selected_level"] in (3,4) || error("no accepted common native level")
    selected=SG.global_packet(ctx,run,"native-selection")
    selected["selected_level"]==native["selected_level"] || error("stale native selection")
    refs=SG.point_packets(ctx,run,"tailref123")
    first_derivative=sidis_job_derivative_decision(SG.all_ids(ctx),refs,native["selected_level"])
    if phase=="derivative-select3"
        merge!(report,first_derivative)
        return report
    end
    fourth_refs=first_derivative["selected_reference_level"]==0 ? SG.point_packets(ctx,run,"reference4") : nothing
    derivative=fourth_refs===nothing ? first_derivative : sidis_job_derivative_decision(SG.all_ids(ctx),refs,native["selected_level"];fourth=fourth_refs)
    derivative["selected_reference_level"] in (3,4) || error("all132/all13 inward reference policy exhausted without acceptance")
    math=SG.global_packet(ctx,run,"source-math");SG.HistoricalSIDISSourceMath.audit(math,ctx)
    report["native_level"]=native["selected_level"];report["reference_level"]=derivative["selected_reference_level"]
    report["coordinates"]=13;report["full132_native_tail_reference_derivative_policy_completed"]=true
    report["points"]=SG.expected_rows(ctx,refs,native["selected_level"],derivative["selected_reference_level"])
    report["job_ledger"]=SG.ledger(ctx,run,phase)
    report["scope"]=ctx.contract["scope"]
    report["accepted_prediction_diagnostic"]=report["all132_coverage_verified"]=true
    report["status"]="CONDITIONAL_FULL132_POINT_DIAGNOSTIC_RECORDED_NOT_FITTING"
    report
end

function main()
    length(ARGS)==6 || error("REQUEST PHASE NATIVE_BUILD PARTON_BUILD RUN FRESH_REPORT")
    request,phase,native,inputs,run,output=ARGS
    phase in SG.GLOBAL_PHASES || error("unknown global SIDIS phase")
    abspath(output)==SG.global_report(abspath(run),phase) && !ispath(output) || error("canonical fresh global report required")
    report=Dict{String,Any}("status"=>"FAIL","phase"=>phase,"scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,SG.false_flags()...)
    try
        ctx=SG.bindings(request,native,inputs);merge!(report,SG.header(ctx,phase))
        populate_global!(report,ctx,abspath(run),phase);SG.stable(ctx)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
        for key in SG.FLAGS;report[key]=false;end
    finally
        report["completed_at"]=string(now(UTC));writejson(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"phase"=>phase,"error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]=="FAIL" && exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
