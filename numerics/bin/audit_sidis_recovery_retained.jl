include("SIDISRecoveryBindings.jl")
const RB=SIDISRecoveryBindings
length(ARGS)==5 || error("FROZEN PREDECESSOR MANIFEST MANIFEST_SHA FRESH_AUDIT")
const RM=RB.load_manifest(ARGS[1:4]...)
include(joinpath(RM.frozen,"numerics/bin/HistoricalSIDISJobRecords.jl"))
include("SIDISRecoveryDispatch.jl")
const SR=HistoricalSIDISJobRecords

function recovery_main()
    output=ARGS[5];RB.require_new_output(RM,output)
    record=Dict{String,Any}("status"=>"FAIL","scope"=>"saved_retained_arithmetic_only_not_cold_numerics_or_current_global_acceptance",
        "recovery"=>RB.provenance(RM),"current_global_prerequisites_accepted"=>false,"jobs_reviewed"=>Any[],
        "retained_invocations"=>Dict{String,Any}(),SR.false_flags()...)
    try
        ctx=SR.bindings(RM.data["request"],RM.data["native"],RM.data["parton"])
        for row in RM.data["retained_jobs"]
            phase=row["phase"];index=row["index"];producer=RB.origin(RM.run,phase,index)
            SR.point_folder(producer,phase,index)==row["folder"] || error("wrong actual retained producing path")
            packet=SR.audit_point(ctx,producer,phase,index;complete_invocation=true)
            check_recovery_transport(RM.run,phase,index,packet;complete_invocation=true)
            path=joinpath(row["folder"],"invocation.json")
            record["retained_invocations"][path]=SR.digest(path)
            push!(record["jobs_reviewed"],[phase,index])
        end
        SR.stable(ctx);RB.stable(RM)
        record["source_fingerprint"]=ctx.fp;record["request_sha256"]=ctx.request.sha256
        record["status"]="RETAINED_POINT_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE"
        record["unreexecuted_producer_measurements"]=["native evolution and rows","physical reference assembly","AD","old timings and RSS"]
    catch exception
        record["status"]="FAIL";record["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        record["completed_at"]=string(now(UTC));RB.write_record(output,record)
    end
    println(JSON.json(Dict("status"=>record["status"],"pairs_reviewed"=>length(record["jobs_reviewed"]),"error"=>get(record,"error",nothing))));flush(stdout)
    record["status"]=="RETAINED_POINT_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE" || exit(2)
end
recovery_main()
