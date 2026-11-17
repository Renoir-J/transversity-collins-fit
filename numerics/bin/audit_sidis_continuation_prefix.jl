include("SIDISContinuationPaths.jl")
const CP=SIDISContinuationPaths
length(ARGS)==4 || error("FROZEN_CHECKOUT MANIFEST MANIFEST_SHA FRESH_AUDIT")
const CM=CP.load_manifest(ARGS[1],ARGS[2],ARGS[3])
using CollinsFit, JSON, Dates
include(joinpath(CM.frozen,"numerics/bin/HistoricalSIDISJobRecords.jl"))
include("SIDISContinuationDispatch.jl")
const SP=HistoricalSIDISJobRecords

function continuation_main()
    output=ARGS[4];CP.require_new_output(CM,output)
    record=Dict{String,Any}("status"=>"FAIL","scope"=>CP.SCOPE,"continuation"=>CP.provenance(CM),
        "points_reviewed"=>Int[],"original_invocations"=>Dict{String,String}(),CP.false_flags()...)
    try
        ctx=SP.bindings(CM.data["request"],CM.data["native"],CM.data["parton"])
        for i in 1:13
            folder=SP.point_folder(CM.run,"native123",i)
            folder==SP.point_folder(CM.run.original,"native123",i) || error("wrong original prefix resolution")
            SP.audit_point(ctx,CM.run,"native123",i;complete_invocation=true)
            push!(record["points_reviewed"],i)
            record["original_invocations"][joinpath(folder,"invocation.json")]=SP.digest(joinpath(folder,"invocation.json"))
        end
        SP.point_folder(CM.run,"native123",14)==joinpath(CM.run.root,"native123/p0014") || error("partial point14 must not resolve to old data")
        SP.stable(ctx);CP.stable(CM)
        record["source_fingerprint"]=ctx.fp;record["request_sha256"]=ctx.request.sha256
        record["native_evolution_AD_or_original_resources_reexecuted"]=false
        record["original_files_verified_unchanged"]=141
        record["status"]="ORIGINAL13_TYPED_PATH_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE"
    catch exception
        record["status"]="FAIL";record["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        record["completed_at"]=string(now(UTC));CP.write_record(output,record)
    end
    println(JSON.json(Dict("status"=>record["status"],"error"=>get(record,"error",nothing))));flush(stdout)
    record["status"]=="FAIL" && exit(2)
end
continuation_main()
