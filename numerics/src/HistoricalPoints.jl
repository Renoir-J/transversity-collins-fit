"Source-bound requests for fixed-parameter historical point diagnostics, without likelihood admission."
const HISTORICAL_POINT_SCOPE="historical_reconstructed_point_prediction_diagnostic"
const HISTORICAL_POINT_CONTRACT=normpath(joinpath(@__DIR__,"../../config/historical-point-contract.json"))
const HISTORICAL_INITIAL_CONTRACT=normpath(joinpath(@__DIR__,"../../config/initial-contract.json"))

struct HistoricalPointRequest
    data::Dict{String,Any}
    path::String
    sha256::String
    signature::String
end

historical_geometry_key(point)=Tuple(point["kinematics"][key] for key in ("x","y","z","Q2_GeV2","pT_GeV"))

function historical_request_valid(request::HistoricalPointRequest)
    numerical_snapshot(request.data)==request.signature || error("mutated historical point request")
    bytes2hex(sha256(read(request.path)))==request.sha256 || error("changed historical point request bytes")
    for (path,digest) in request.data["source_inputs"]
        isfile(path) && bytes2hex(sha256(read(path)))==digest || error("changed historical point source input: $path")
    end
    request
end

function load_historical_point_request(path)
    data=JSON.parsefile(path);contract=JSON.parsefile(HISTORICAL_POINT_CONTRACT)
    initial=JSON.parsefile(HISTORICAL_INITIAL_CONTRACT)
    data isa AbstractDict && get(data,"schema_version",nothing) isa Integer &&
        !(data["schema_version"] isa Bool) && data["schema_version"]==1 &&
        get(data,"scope",nothing)==HISTORICAL_POINT_SCOPE && get(data,"theory_scope",nothing)=="reconstructed_theory" &&
        get(data,"observation_scope",nothing)=="historical_benchmark_candidate" || error("wrong historical diagnostic request scope")
    for key in ("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted","fit_performed")
        get(data,key,nothing)===false && get(contract,key,nothing)===false || error("historical diagnostic cannot approve fitting or production")
    end
    data["contract_id"]==contract["id"] && data["contract_sha256"]==bytes2hex(sha256(read(HISTORICAL_POINT_CONTRACT))) &&
        data["initial_contract_sha256"]==bytes2hex(sha256(read(HISTORICAL_INITIAL_CONTRACT))) || error("stale historical point contract")
    parameters=data["parameters"]
    parameters["order"]==initial["parameter_order"] && parameters["values"]==initial["published_parameters"] &&
        parameters["values"]==collect(kpsy_published_parameters().values) &&
        all(x->x isa Real && !(x isa Bool) && isfinite(x),parameters["values"]) || error("changed fixed published parameter source")
    convention=data["model_convention"]
    all(convention[key]==contract["prediction"][key] for key in
        ("phase","normalization","model_y","source_mapping_status")) || error("changed conditional model convention")
    convention["phase"]=="sin_phi_sum" && convention["normalization"]=="fixed_y_DNN_included" || error("unreviewed historical model convention")
    inputs=data["source_inputs"]
    inputs isa AbstractDict && !isempty(inputs) || error("historical request requires input bindings")
    source=data["source_preparation"]
    source["status"]=="CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED" &&
        source["scope"]=="historical_benchmark_candidate" &&
        source["verification_kind"]=="Python_source_preparation_attestation_not_Julia_raw_table_replay" ||
        error("source preparation attestation required")
    candidate_path=joinpath(source["candidates"],"candidates.json")
    report_path=joinpath(source["candidates"],"report.json")
    get(inputs,candidate_path,nothing)==source["candidate_sha256"] &&
        get(inputs,report_path,nothing)==source["report_sha256"] || error("missing historical candidate/report bindings")
    candidate=JSON.parsefile(candidate_path);report=JSON.parsefile(report_path)
    candidate["scope"]=="historical_benchmark_candidate" && report["status"]=="CANDIDATES_VALIDATED" &&
        candidate["source_fingerprint"]==report["source_fingerprint"]==data["source_fingerprint"] || error("stale historical candidate source")
    for item in (candidate,report),key in ("production_ready","physical_extraction","historical_fit_ready")
        item[key]===false || error("historical candidate cannot self-admit")
    end
    views=Dict(view["id"]=>view for view in candidate["views"])
    length(views)==length(candidate["views"]) || error("duplicate historical candidate view")
    points=data["points"];selection=contract["selection"]
    length(points)==selection["total_points"] || error("incomplete historical point request")
    expected_ids=String[]
    for selected in selection["source_views"]
        view=get(views,selected["id"],nothing)
        view!==nothing && view["source_id"]==selected["id"] && view["process"]=="SIDIS" && view["projection"]=="x" &&
            view["publication"]==selection["publication"] && view["target"]=="neutron" &&
            view["target_material"]=="3He" && view["target_result"]==selection["target_result"] &&
            view["source_sha256"]==selected["sha256"] && view["hadron"]==selected["hadron"] &&
            view["source_observable"]=="AUTcollins" && view["output_units"]=="fraction" &&
            view["production_ready"]===false && view["historical_fit_ready"]===false && !isempty(view["blockers"]) &&
            selected["depolarization_removed"]===false && view["depolarization_removed"]===selected["depolarization_removed"] &&
            length(view["points"])==selected["points"] &&
            isempty(view["excluded"]) || error("mismatched primary HallA selection")
        for (offset,row) in enumerate(view["points"])
            row["id"]==selected["id"]*":r"*string(offset+1) && row["locator"]["excel_row"]==offset+1 ||
                error("changed source row identity or order")
            push!(expected_ids,row["id"])
            index=length(expected_ids);point=points[index]
            point["id"]==row["id"] && point["dataset_id"]==selected["id"] &&
                point["charge"] isa Integer && !(point["charge"] isa Bool) && point["charge"]==selected["charge"] &&
                point["hadron"]==selected["hadron"] && point["target"]=="neutron" &&
                point["target_material"]=="3He" && point["target_result"]==selection["target_result"] &&
                point["publication"]==view["publication"] && point["source_sha256"]==selected["sha256"] &&
                point["kinematics"]==row["averages"] && point["source_locator"]==row["locator"] &&
                point["source_remap"]==row["source_remap"] && point["source_blockers"]==view["blockers"] &&
                point["original_depolarization_removed"]===view["depolarization_removed"] || error("changed source-derived historical point")
            all(point["measured"][key]==row[key] for key in
                ("central_value","statistical_error","systematic_components")) || error("changed measured historical scalar/error")
            k=point["kinematics"]
            Set(keys(k))==Set(("x","y","z","Q2_GeV2","pT_GeV")) || error("complete historical means required")
            all(value->value isa Real && !(value isa Bool) && isfinite(value),values(k)) &&
                all(0<k[key]<1 for key in ("x","y","z")) && k["Q2_GeV2"]>0 && k["pT_GeV"]>=0 || error("invalid historical means")
            m=point["measured"]
            remap=point["source_remap"]
            remap["kind"]=="reviewed_legacy_column_remap" && remap["primary_tables"]==["table:kine","table:results"] &&
                length(remap["primary_lines"])==2 || error("primary source remap required")
            Set(keys(m["systematic_components"]))==Set(("experimental","model")) &&
                m["central_value"] isa Real && !(m["central_value"] isa Bool) && isfinite(m["central_value"]) &&
                m["statistical_error"] isa Real && !(m["statistical_error"] isa Bool) && isfinite(m["statistical_error"]) && m["statistical_error"]>0 &&
                all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>=0,values(m["systematic_components"])) || error("invalid separate source errors")
        end
    end
    length(unique(expected_ids))==length(expected_ids) || error("duplicate historical point identity")
    owned=deepcopy(Dict{String,Any}(data))
    request=HistoricalPointRequest(owned,abspath(path),bytes2hex(sha256(read(path))),numerical_snapshot(owned))
    historical_request_valid(request)
end

"Compose the existing joint prediction interface; source data never enter the parameter calculation."
function historical_point_model(request::HistoricalPointRequest,contexts)
    historical_request_valid(request)
    observations=JointObservation[]
    for point in request.data["points"]
        key=historical_geometry_key(point)
        haskey(contexts,key) || error("missing exact historical geometry context")
        context=contexts[key];k=point["kinematics"];plan=context.leptonic
        plan.process==:SIDIS && plan.Q2==k["Q2_GeV2"] && plan.P==k["pT_GeV"] && plan.z1==k["z"] &&
            collect(context.first_operator.queries)==[k["x"]] && collect(context.second_operator.queries)==[k["z"]] &&
            context.target.name=="neutron" && context.target.unpolarized==(0.,1.) &&
            context.target.polarized==(0.,1.) || error("historical context geometry/target mismatch")
        charged=with_pion_charges(context,point["charge"])
        convention=SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:fixed_y_DNN_included,y=k["y"],
            source="HallA Eq1 and KPSY conditional massless fixed-mean model; request SHA256="*request.sha256)
        observable=sidis_observable(charged,convention)
        push!(observations,JointObservation(point["id"],observable;
            source=point["publication"]*" source row; conditional diagnostic, no experimental admission"))
    end
    JointPredictionModel(observations;source="source-bound historical fixed-mean diagnostic; request SHA256="*request.sha256)
end
