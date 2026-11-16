using CollinsFit, JSON, SHA, Dates

const PILOT_AUDIT_ROOT=normpath(joinpath(@__DIR__,"../.."))
pilot_audit_digest(path)=bytes2hex(sha256(read(path)))

function audit_sidis_pilot!(audit,requestfile,case_id,reportfile)
    policy=historical_sidis_pilot_policy();request=load_historical_sidis_point_request(requestfile)
    point=sidis_pilot_case(request,policy,case_id);report=JSON.parsefile(reportfile)
    report["status"]==policy["reporting"]["successful_status"] && report["scope"]==policy["scope"] &&
        report["case_id"]==case_id && report["theory_scope"]=="reconstructed_theory" &&
        report["observation_scope"]=="historical_benchmark_candidate" || error("incomplete or wrong-scope pilot report")
    all(get(report,k,nothing)===false for k in CollinsFit.HISTORICAL_SIDIS_PILOT_FALSE_FLAGS) &&
        report["full132_native_tail_reference_derivative_policy_completed"]===false || error("pilot report attempts diagnostic admission")
    isequal(report["point"],point) && report["request_points"]==132 &&
        isequal(report["parameters"],request.data["parameters"]) &&
        report["request_sha256"]==request.sha256 || error("changed source identity or parameter vector")
    report["policy_sha256"]==pilot_audit_digest(CollinsFit.HISTORICAL_SIDIS_PILOT_POLICY) &&
        report["parent_contract_sha256"]==request.data["contract_sha256"] &&
        isequal(report["numerical_setting"],policy["numerics"]) || error("changed frozen pilot settings")
    sources=sidis_pilot_current_sources(PILOT_AUDIT_ROOT)
    report["sources"]==sources && report["source_fingerprint"]==request.data["source_fingerprint"]==
        sidis_pilot_source_fingerprint(sources) || error("stale or incomplete pilot source bindings")
    native=report["native_build_receipt"];inputs=report["parton_build_receipt"]
    Set(keys(report["dependency_sha256"]))==Set((requestfile,native,inputs)) &&
        all(isfile(p) && pilot_audit_digest(p)==h for (p,h) in report["dependency_sha256"]) || error("changed pilot dependency bindings")
    CollinsFit.verified_hoppet_build(native)==report["native_build"] &&
        CollinsFit.verify_parton_build(inputs)==report["input_build"] || error("changed native or parton build receipt")
    checks=sidis_pilot_value_checks(report["values"],report["independent_reference"],point,policy)
    isequal(checks,report["checks"]) && all(c["passed"] for c in checks) || error("pilot saved scalar arithmetic failed")
    report["source_factor"]==sidis_pilot_source_factor(point) || error("changed pilot source factor")
    refscope=report["reference_scope"]
    refscope["native_rows_shared"]===true && refscope["initial_functions_shared"]===true &&
        refscope["independent_evolution"]===false && refscope["derivative_reference_executed"]===false &&
        refscope["family"]=="generic_fixed_row_scalar_assembly_and_adaptive_P_reference" &&
        CollinsFit.sidis_source_number(refscope["preparation_seconds_informational"]) &&
        refscope["preparation_seconds_informational"]>=0 || error("misstated pilot reference independence")
    samples=report["AD_samples"]
    samples["independent_derivative_reference_executed"]===samples["derivative_accuracy_accepted"]===false || error("pilot AD cannot self-certify")
    parent=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT)
    samples["direction"]==parent["derivative_policy"]["direction"] &&
        all(length(samples[k])==13 && all(CollinsFit.sidis_source_number,samples[k]) for k in ("gradient","HVP")) || error("invalid all13 AD timing samples")
    repetitions=report["evaluation_times"]["warm_repetitions"]
    length(repetitions)==policy["numerics"]["derivative_evaluations"]["warm_repetitions"] &&
        all(r["index"]==i && r["value"]==report["values"]["amplitude"] && r["gradient"]==samples["gradient"] &&
            r["HVP"]==samples["HVP"] for (i,r) in enumerate(repetitions)) || error("unrepeatable pilot samples")
    report["finite_AD_samples"]===true && report["repeatable_samples"]===true || error("missing finite/repeatable pilot samples")
    for times in vcat([report["evaluation_times"]["first_call_seconds"]],[r["seconds"] for r in repetitions])
        Set(keys(times))==Set(("value","gradient","HVP")) &&
            all(v->CollinsFit.sidis_source_number(v) && v>=0,values(times)) || error("invalid timing observation")
    end
    resource=report["resources"]
    all(resource[k] isa Integer && !(resource[k] isa Bool) && resource[k]>0 for k in ("radial_nodes","row_coefficient_payload_bytes")) &&
        CollinsFit.sidis_source_number(resource["cold_preparation_seconds_informational"]) &&
        resource["cold_preparation_seconds_informational"]>=0 || error("invalid producer resource observations")
    provenance=resource["context_provenance"];target=request.data["target_weights"][point["target"]]
    provenance["scope"]=="reconstructed_theory" && provenance["production_ready"]===false &&
        provenance["physical_extraction"]===false && provenance["process"]=="SIDIS" &&
        provenance["radial_bound_inverse_GeV"]==policy["numerics"]["radial_bound_inverse_GeV"] &&
        provenance["range_policy"]=="reference" && provenance["target"]["name"]==point["target"] &&
        provenance["target"]["unpolarized"]==target["unpolarized"] &&
        provenance["target"]["polarized"]==target["polarized"] || error("changed producer context provenance")
    historical_sidis_request_valid(request)
    audit["status"]=policy["reporting"]["audit_status"];audit["source_fingerprint"]=request.data["source_fingerprint"]
    audit["request_sha256"]=request.sha256;audit["report_sha256"]=pilot_audit_digest(reportfile)
    audit["policy_sha256"]=report["policy_sha256"];audit["case_id"]=case_id;audit["point_id"]=point["id"]
    audit["checks"]=checks;audit["saved_scalar_checks_recomputed"]=length(checks)
    audit["AD_samples_checked_for_shape_finiteness_and_repeatability_only"]=true
    audit["unreexecuted_producer_measurements"]=["native rows/evolution","independent scalar assembly",
        "AD/derivative accuracy","coefficient payload","phase timings","RSS/resource measurements"]
    audit["primary_raw_tables_reparsed"]=false
    audit
end

function main()
    length(ARGS)==4 || error("usage: audit_historical_sidis_pilot.jl REQUEST CASE REPORT FRESH_AUDIT")
    request,case_id,report,output=ARGS;ispath(output) && error("use a fresh pilot audit output")
    audit=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_PILOT_SCOPE,
        (k=>false for k in CollinsFit.HISTORICAL_SIDIS_PILOT_FALSE_FLAGS)...,
        "reader_scope"=>"current source/build bindings and saved scalar/repeatability arithmetic only")
    try
        audit_sidis_pilot!(audit,request,case_id,report)
    catch exception
        audit["status"]="FAIL";audit["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        audit["completed_at"]=string(now(UTC));writejson(output,audit)
    end
    println(JSON.json(Dict("status"=>audit["status"],"scope"=>audit["scope"],"error"=>get(audit,"error",nothing))));flush(stdout)
    audit["status"]=="PILOT_REVIEWED_NOT_DIAGNOSTIC_GATE" || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
