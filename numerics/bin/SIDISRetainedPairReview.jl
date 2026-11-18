"Whole saved-pair reader with self-contained report finalization and failure retention."
module SIDISRetainedPairReview
using JSON, Dates

function run(manifest,output,source_api,bindings,transport_check)
    bindings.require_new_output(manifest,output)
    record=Dict{String,Any}("status"=>"FAIL","scope"=>"saved_retained_arithmetic_only_not_cold_numerics_or_current_global_acceptance",
        "theory_scope"=>get(manifest.data,"theory_scope","reconstructed_theory"),
        "recovery"=>bindings.provenance(manifest),"current_global_prerequisites_accepted"=>false,"jobs_reviewed"=>Any[],
        "retained_invocations"=>Dict{String,Any}(),source_api.false_flags()...)
    try
        ctx=source_api.bindings(manifest.data["request"],manifest.data["native"],manifest.data["parton"])
        for row in manifest.data["retained_jobs"]
            phase=row["phase"];index=row["index"];producer=bindings.origin(manifest.run,phase,index)
            source_api.point_folder(producer,phase,index)==row["folder"] || error("wrong actual retained producing path")
            packet=source_api.audit_point(ctx,producer,phase,index;complete_invocation=true)
            transport_check(manifest.run,phase,index,packet;complete_invocation=true)
            path=joinpath(row["folder"],"invocation.json")
            record["retained_invocations"][path]=source_api.digest(path)
            push!(record["jobs_reviewed"],[phase,index])
        end
        source_api.stable(ctx);bindings.stable(manifest)
        record["source_fingerprint"]=ctx.fp;record["request_sha256"]=ctx.request.sha256
        record["status"]="RETAINED_POINT_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE"
        record["unreexecuted_producer_measurements"]=["native evolution and rows","physical reference assembly","AD","old timings and RSS"]
    catch exception
        record["status"]="FAIL";record["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        record["completed_at"]=string(now(UTC));bindings.write_record(output,record)
    end
    println(JSON.json(Dict("status"=>record["status"],"pairs_reviewed"=>length(record["jobs_reviewed"]),"error"=>get(record,"error",nothing))));flush(stdout)
    record["status"]=="RETAINED_POINT_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE" ? 0 : 2
end
end
