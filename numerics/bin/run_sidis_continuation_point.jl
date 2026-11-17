include("SIDISContinuationPaths.jl")
const CP=SIDISContinuationPaths
length(ARGS)==6 || error("FROZEN_CHECKOUT MANIFEST MANIFEST_SHA PHASE INDEX FRESH_REPORT")
const CM=CP.load_manifest(ARGS[1],ARGS[2],ARGS[3])
include(joinpath(CM.frozen,"numerics/bin/run_historical_sidis_point_job.jl"))
include("SIDISContinuationDispatch.jl")

function continuation_main()
    phase=ARGS[4];index=parse(Int,ARGS[5]);output=ARGS[6]
    phase in SJ.POINT_PHASES || error("unknown point phase")
    phase=="native123" && index<=13 && error("never dispatch an original completed point")
    CP.require_new_output(CM,output;canonical=SJ.point_report(CM.run,phase,index))
    BLAS.set_num_threads(1)
    report=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "phase"=>phase,"point_index"=>index,"continuation"=>CP.provenance(CM),SJ.false_flags()...)
    try
        ctx=SJ.bindings(CM.data["request"],CM.data["native"],CM.data["parton"])
        merge!(report,SJ.header(ctx,phase;index))
        run_sidis_point!(report,ctx,phase,index,CM.run,output*".progress.json")
        SJ.stable(ctx);CP.stable(CM)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));CP.write_record(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"phase"=>phase,"point_index"=>index,"error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]==SJ.RECORDED || exit(2)
end
continuation_main()
