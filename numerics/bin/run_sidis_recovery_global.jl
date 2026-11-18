include("SIDISRecoveryBindings.jl")
const RB=SIDISRecoveryBindings
length(ARGS)==7 || error("FROZEN PREDECESSOR MANIFEST MANIFEST_SHA PHASE ZERO FRESH_REPORT")
const RM=RB.load_manifest(ARGS[1:4]...)
include(joinpath(RM.frozen,"numerics/bin/run_historical_sidis_global_job.jl"))
include("SIDISRecoveryDispatch.jl")

function recovery_main()
    phase=ARGS[5];index=parse(Int,ARGS[6]);output=ARGS[7]
    phase in SG.GLOBAL_PHASES && index==0 || error("declared fresh global phase/index0 required")
    RB.require_new_output(RM,output;canonical=SG.global_report(RM.run,phase))
    report=Dict{String,Any}("status"=>"FAIL","phase"=>phase,"scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "recovery"=>RB.provenance(RM),SG.false_flags()...)
    try
        ctx=SG.bindings(RM.data["request"],RM.data["native"],RM.data["parton"])
        merge!(report,SG.header(ctx,phase))
        populate_global!(report,ctx,RM.run,phase)
        SG.stable(ctx);RB.stable(RM)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
        for flag in SG.FLAGS;report[flag]=false;end
    finally
        report["completed_at"]=string(now(UTC));RB.write_record(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"phase"=>phase,"error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]=="FAIL" && exit(2)
end
recovery_main()
