using CollinsFit, JSON, Dates
include("FixedContextStorageRecords.jl")
const FS=FixedContextStorageRecords

function main(args=ARGS)
    length(args)==7 || error("usage: audit_fixed_context_storage_case.jl SIDIS_REQUEST SIA_REQUEST NATIVE_BUILD PARTON_BUILD CASE_ID REPORT FRESH_AUDIT")
    sidis,sia,native,parton,case_id,reportfile,output=args
    ispath(output) && error("fresh storage audit path required")
    c=FS.bindings(sidis,sia,native,parton);s=FS.selected(c,case_id)
    audit=FS.header(c,s);audit["status"]="FAIL"
    try
        report=JSON.parsefile(reportfile);report["status"]=="COLD_RELOADED_FIXED_CONTEXT_BITS_RECORDED_NOT_PHYSICS_ADMISSION" || error("unaccepted storage producer")
        FS.valid_header(report,c,s)
        context=load_fixed_leptonic_context(joinpath(dirname(reportfile),"context");
            expected_manifest_sha256=report["manifest_sha256"],bindings=s.bound,input_kind=:pinned_inputs)
        signature=fixed_leptonic_storage_signature(context)
        signature.logical_sha256==report["logical_context_sha256"] && signature.row_alias_sha256==report["row_alias_sha256"] || error("changed fixed context identity")
        metrics=FS.metrics(c,s,context)
        FS.metrics_valid(report["loaded_metrics"],report["cold_metrics"],s)
        FS.metrics_valid(metrics,report["cold_metrics"],s);FS.stable(c)
        audit["producer_report_sha256"]=FS.digest(reportfile)
        audit["manifest_sha256"]=report["manifest_sha256"]
        audit["loaded_metrics"]=metrics
        audit["loaded_observables_reexecuted"]=length(s.points)
        audit["full13_gradient_and_directional_HVP_reexecuted"]=true
        audit["cold_native_preparation_reexecuted"]=false
        audit["whole_command_resources_or_temporal_reuse_replayed"]=false
        audit["status"]="FRESH_LOADED_CONTEXT_BITS_AND_SOURCES_REVIEWED_NOT_PHYSICS_ADMISSION"
    catch error
        audit["error"]=sprint(showerror,error,catch_backtrace())
    finally
        audit["completed_at"]=string(now(UTC));writejson(output,audit)
    end
    println(JSON.json(Dict("status"=>audit["status"],"case_id"=>case_id,"production_ready"=>false)))
    audit["status"]=="FRESH_LOADED_CONTEXT_BITS_AND_SOURCES_REVIEWED_NOT_PHYSICS_ADMISSION" ? 0 : 2
end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(main())
end
