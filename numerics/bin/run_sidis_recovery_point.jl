include("SIDISRecoveryBindings.jl")
const RB=SIDISRecoveryBindings
length(ARGS)==7 || error("FROZEN PREDECESSOR MANIFEST MANIFEST_SHA PHASE INDEX FRESH_REPORT")
const RM=RB.load_manifest(ARGS[1:4]...)
include(joinpath(RM.frozen,"numerics/bin/run_historical_sidis_point_job.jl"))
include("SIDISRecoveryDispatch.jl")

function recovery_main()
    phase=ARGS[5];index=parse(Int,ARGS[6]);output=ARGS[7]
    RB.RP.point_index(phase,index)
    RB.RP.retained(RM.run.paths,phase,index)===nothing || error("never dispatch a retained complete pair")
    RB.require_new_output(RM,output;canonical=SJ.point_report(RM.run,phase,index))
    BLAS.set_num_threads(1)
    report=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "phase"=>phase,"point_index"=>index,"recovery"=>RB.provenance(RM),SJ.false_flags()...)
    try
        ctx=SJ.bindings(RM.data["request"],RM.data["native"],RM.data["parton"])
        merge!(report,SJ.header(ctx,phase;index))
        run_sidis_point!(report,ctx,phase,index,RM.run,output*".progress.json")
        SJ.stable(ctx);RB.stable(RM)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));RB.write_record(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"phase"=>phase,"point_index"=>index,"error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]==SJ.RECORDED || exit(2)
end
recovery_main()
