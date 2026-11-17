include("SIDISContinuationPaths.jl")
const CP=SIDISContinuationPaths
length(ARGS)==6 || error("FROZEN_CHECKOUT MANIFEST MANIFEST_SHA PHASE ZERO FRESH_REPORT")
const CM=CP.load_manifest(ARGS[1],ARGS[2],ARGS[3])
include(joinpath(CM.frozen,"numerics/bin/run_historical_sidis_global_job.jl"))
include("SIDISContinuationDispatch.jl")

function continuation_main()
    phase=ARGS[4];index=parse(Int,ARGS[5]);output=ARGS[6]
    phase in SG.GLOBAL_PHASES && index==0 || error("declared global phase/index0 required")
    CP.require_new_output(CM,output;canonical=SG.global_report(CM.run,phase))
    report=Dict{String,Any}("status"=>"FAIL","phase"=>phase,"scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "continuation"=>CP.provenance(CM),SG.false_flags()...)
    try
        ctx=SG.bindings(CM.data["request"],CM.data["native"],CM.data["parton"])
        merge!(report,SG.header(ctx,phase))
        populate_global!(report,ctx,CM.run,phase)
        SG.stable(ctx);CP.stable(CM)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
        for flag in SG.FLAGS;report[flag]=false;end
    finally
        report["completed_at"]=string(now(UTC));CP.write_record(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"phase"=>phase,"error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]=="FAIL" && exit(2)
end
continuation_main()
