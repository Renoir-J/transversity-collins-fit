using CollinsFit, JSON, Dates
include("FixedContextStorageRecords.jl")
const FS=FixedContextStorageRecords

function main(args=ARGS)
    length(args)==6 || error("usage: run_fixed_context_storage_case.jl SIDIS_REQUEST SIA_REQUEST NATIVE_BUILD PARTON_BUILD CASE_ID FRESH_REPORT")
    sidis,sia,native,parton,case_id,output=args
    any(ispath,(output,output*".progress.json",joinpath(dirname(output),"context"))) && error("fresh storage case paths required")
    c=FS.bindings(sidis,sia,native,parton);s=FS.selected(c,case_id)
    report=FS.header(c,s);report["status"]="IN_PROGRESS_NOT_A_GATE"
    function progress(stage)
        report["current_stage"]=stage;writejson(output*".progress.json",report)
        println(JSON.json(Dict("case_id"=>case_id,"stage"=>stage)));flush(stdout)
    end
    try
        progress("cold fixed native/source preparation")
        cold=nothing;prepare_seconds=@elapsed cold=FS.prepare(c,s)
        cold_signature=fixed_leptonic_storage_signature(cold)
        progress("cold value/full13 gradient/directional HVP")
        original=nothing;cold_metric_seconds=@elapsed original=FS.metrics(c,s,cold)
        progress("write and verify complete fixed data")
        receipt=nothing;save_seconds=@elapsed receipt=save_fixed_leptonic_context(cold,joinpath(dirname(output),"context");bindings=s.bound,input_kind=:pinned_inputs)
        progress("reload owned context")
        loaded=nothing;load_seconds=@elapsed loaded=load_fixed_leptonic_context(joinpath(dirname(output),"context");
            expected_manifest_sha256=receipt.manifest_sha256,bindings=s.bound,input_kind=:pinned_inputs)
        loaded_signature=fixed_leptonic_storage_signature(loaded)
        cold_signature==loaded_signature || error("fixed data or row aliases changed")
        progress("loaded value/full13 gradient/directional HVP")
        reloaded=nothing;loaded_metric_seconds=@elapsed reloaded=FS.metrics(c,s,loaded)
        FS.metrics_valid(reloaded,original,s);FS.valid_header(report,c,s);FS.stable(c)
        report["manifest_sha256"]=receipt.manifest_sha256;report["payload_sha256"]=receipt.payload_sha256
        report["logical_context_sha256"]=receipt.logical_context_sha256;report["row_alias_sha256"]=receipt.row_alias_sha256
        report["payload_bytes"]=receipt.payload_bytes;report["stored_objects"]=receipt.objects
        report["cold_metrics"]=original;report["loaded_metrics"]=reloaded
        report["cold_preparation_and_storage_seconds_informational"]=Dict("cold_preparation"=>prepare_seconds,
            "cold_metrics_including_first_compilation"=>cold_metric_seconds,"save_and_verify"=>save_seconds,
            "load_and_verify"=>load_seconds,"loaded_metrics"=>loaded_metric_seconds)
        report["timing_scope"]="informational shared-workload samples; not production throughput or a full-fit estimate"
        report["status"]="COLD_RELOADED_FIXED_CONTEXT_BITS_RECORDED_NOT_PHYSICS_ADMISSION"
    catch error
        report["status"]="FAIL";report["error"]=sprint(showerror,error,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"case_id"=>case_id,"production_ready"=>false)))
    report["status"]=="COLD_RELOADED_FIXED_CONTEXT_BITS_RECORDED_NOT_PHYSICS_ADMISSION" ? 0 : 2
end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(main())
end
