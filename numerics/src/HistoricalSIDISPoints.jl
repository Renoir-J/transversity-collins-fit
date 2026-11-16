"Source-bound full132 SIDIS requests; model composition does not admit a likelihood."
const HISTORICAL_SIDIS_SCOPE="historical_reconstructed_sidis_point_prediction_diagnostic"
const HISTORICAL_SIDIS_CONTRACT=normpath(joinpath(@__DIR__,"../../config/historical-sidis-point-contract.json"))
const HISTORICAL_SIDIS_FALSE_FLAGS=("production_ready","physical_extraction","historical_fit_ready",
    "experimental_data_admitted","fit_performed","G2PERF_pass","interval_coverage_verified")

struct HistoricalSIDISPointRequest
    data::Dict{String,Any}
    path::String
    sha256::String
    signature::String
end

historical_sidis_geometry_key(point)=(point["target"],(point["kinematics"][k] for k in ("x","z","Q2_GeV2","pT_GeV"))...)

function historical_sidis_request_valid(request::HistoricalSIDISPointRequest)
    request.data["scope"]==HISTORICAL_SIDIS_SCOPE &&
        all(get(request.data,k,nothing)===false for k in HISTORICAL_SIDIS_FALSE_FLAGS) &&
        get(request.data,"accepted_fitting_gate",nothing)===false &&
        get(request.data,"likelihood_constructed",nothing)===false || error("SIDIS request cannot self-admit")
    numerical_snapshot(request.data)==request.signature || error("mutated full SIDIS request")
    bytes2hex(sha256(read(request.path)))==request.sha256 || error("changed full SIDIS request bytes")
    for (path,digest) in request.data["source_inputs"]
        isfile(path) && bytes2hex(sha256(read(path)))==digest || error("changed SIDIS source input: $path")
    end
    request
end

sidis_source_number(value)=value isa Real && !(value isa Bool) && isfinite(value)

function load_historical_sidis_point_request(path)
    data=JSON.parsefile(path);contract=JSON.parsefile(HISTORICAL_SIDIS_CONTRACT)
    root=normpath(joinpath(@__DIR__,"../.."));initialpath=joinpath(root,"config/initial-contract.json")
    initial=JSON.parsefile(initialpath);digest(p)=bytes2hex(sha256(read(p)))
    data isa AbstractDict && get(data,"schema_version",nothing) isa Integer &&
        !(data["schema_version"] isa Bool) && data["schema_version"]==1 &&
        get(data,"scope",nothing)==HISTORICAL_SIDIS_SCOPE && contract["scope"]==HISTORICAL_SIDIS_SCOPE &&
        get(data,"theory_scope",nothing)=="reconstructed_theory" &&
        get(data,"observation_scope",nothing)=="historical_benchmark_candidate" || error("wrong full SIDIS request scope")
    all(get(data,k,nothing)===false && get(contract,k,nothing)===false for k in HISTORICAL_SIDIS_FALSE_FLAGS) ||
        error("full SIDIS diagnostic cannot approve fitting, performance or production")
    data["contract_id"]==contract["id"]=="KPSY-HISTORICAL-SIDIS-POINT-001" &&
        data["contract_sha256"]==digest(HISTORICAL_SIDIS_CONTRACT) &&
        data["initial_contract_sha256"]==digest(initialpath) || error("stale full SIDIS contract")
    params=data["parameters"]
    params["order"]==initial["parameter_order"]==contract["prediction"]["parameter_order"] &&
        params["values"]==initial["published_parameters"]==contract["prediction"]["parameter_values"]==
            collect(kpsy_published_parameters().values) && all(sidis_source_number,params["values"]) &&
        params["role"]==initial["published_parameter_role"] || error("changed rounded thirteen-parameter source")
    isequal(data["model_conventions"],contract["model_conventions"]) &&
        isequal(data["target_weights"],contract["prediction"]["target_weights"]) || error("changed source convention or target weights")
    inputs=data["source_inputs"]
    inputs isa AbstractDict && !isempty(inputs) || error("full SIDIS source bindings required")
    primarypath=joinpath(root,contract["primary_manifest"]);primary=JSON.parsefile(primarypath)
    data["primary_manifest_sha256"]==contract["primary_manifest_sha256"]==digest(primarypath) || error("stale primary manifest")
    for p in (HISTORICAL_SIDIS_CONTRACT,initialpath,primarypath,joinpath(root,contract["source_manifest"]))
        get(inputs,p,nothing)==digest(p) || error("missing current SIDIS contract/manifest binding")
    end
    for (p,h) in contract["parent_sha256"]
        parent=joinpath(root,p)
        digest(parent)==h && get(inputs,parent,nothing)==h || error("changed or unbound SIDIS parent")
    end
    source=data["source_preparation"]
    source["status"]=="CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED" &&
        source["scope"]=="historical_benchmark_candidate" &&
        source["verification_kind"]=="Python_source_preparation_attestation_not_Julia_raw_table_replay" &&
        source["primary_verification"]=="pinned_convention_bytes_not_numerical_transcription_or_dataset_admission" ||
        error("explicit source preparation attestation required")
    expected_primary=Dict(joinpath(source["primary_cache"],s["filename"])=>s["sha256"] for s in primary["sources"])
    length(expected_primary)==9 && expected_primary==source["primary_inputs"] &&
        all(get(inputs,p,nothing)==h for (p,h) in expected_primary) || error("incomplete pinned primary bytes")
    candidatepath=joinpath(source["candidates"],"candidates.json");reportpath=joinpath(source["candidates"],"report.json")
    get(inputs,candidatepath,nothing)==source["candidate_sha256"]==digest(candidatepath) &&
        get(inputs,reportpath,nothing)==source["report_sha256"]==digest(reportpath) || error("changed or missing source candidate/report binding")
    candidate=JSON.parsefile(candidatepath);report=JSON.parsefile(reportpath)
    candidate["scope"]=="historical_benchmark_candidate" && report["status"]=="CANDIDATES_VALIDATED" &&
        candidate["source_fingerprint"]==report["source_fingerprint"]==data["source_fingerprint"] || error("stale historical source preparation")
    for item in (candidate,report),k in HISTORICAL_SIDIS_FALSE_FLAGS[1:3]
        get(item,k,nothing)===false || error("historical source cannot self-admit")
    end
    all(get(inputs,p,nothing)==h for (p,h) in report["source_hashes"]) || error("unbound historical raw source")
    _validate_historical_sidis_layout(data,contract,candidate)
    owned=deepcopy(Dict{String,Any}(data))
    request=HistoricalSIDISPointRequest(owned,abspath(path),digest(path),numerical_snapshot(owned))
    historical_sidis_request_valid(request)
end

"Representation checks only; raw sources and pinned primary bytes are verified by the loader."
function _validate_historical_sidis_layout(data,contract,candidate)
    views=Dict(v["id"]=>v for v in candidate["views"])
    length(views)==length(candidate["views"]) || error("duplicate candidate view")
    points=data["points"];selection=contract["selection"];expected_ids=String[];excluded=0
    length(points)==selection["points"]==132 && length(selection["source_views"])==selection["views"]==18 || error("incomplete full SIDIS selection")
    Set(keys(data["view_evidence"]))==Set(v["id"] for v in selection["source_views"]) || error("incomplete source-view evidence")
    for selected in selection["source_views"]
        view=get(views,selected["id"],nothing);convention=contract["model_conventions"][selected["publication"]]
        view!==nothing && view["source_id"]==selected["id"] && view["process"]=="SIDIS" &&
            all(isequal(view[k],selected[k]) for k in ("publication","target","hadron","projection","source_sha256")) &&
            view["depolarization_removed"]===selected["original_depolarization_removed"] &&
            convention["id"]==selected["source_convention_id"] &&
            convention["original_depolarization_removed"]===selected["original_depolarization_removed"] &&
            view["source_observable"]=="AUTcollins" && view["output_units"]=="fraction" &&
            view["production_ready"]===false && view["historical_fit_ready"]===false && !isempty(view["blockers"]) &&
            length(view["points"])==selected["points"] || error("mismatched full SIDIS source view")
        evidence=Dict(k=>v for (k,v) in view if k!="points")
        isequal(data["view_evidence"][selected["id"]],evidence) || error("changed source labels, exclusions or source-view provenance")
        excluded+=length(view["excluded"])
        for row in view["points"]
            n=row["locator"]["excel_row"]
            n isa Integer && !(n isa Bool) && n>=2 && row["id"]==selected["id"]*":r"*string(n) || error("changed SIDIS row identity")
            push!(expected_ids,row["id"]);point=points[length(expected_ids)]
            point["id"]==row["id"] && point["dataset_id"]==selected["id"] &&
                all(isequal(point[k],view[k]) for k in ("publication","target","hadron","projection","source_sha256")) &&
                point["charge"] isa Integer && !(point["charge"] isa Bool) &&
                point["charge"]==(view["hadron"]=="pi+" ? 1 : view["hadron"]=="pi-" ? -1 : 0) &&
                point["original_depolarization_removed"]===view["depolarization_removed"] &&
                isequal(point["source_row"],row) && isequal(point["kinematics"],row["averages"]) &&
                isequal(point["source_blockers"],view["blockers"]) && isequal(point["model_convention"],convention) ||
                error("changed source-derived SIDIS point or convention")
            isequal(point["measured"],Dict(k=>row[k] for k in ("central_value","statistical_error","systematic_components"))) ||
                error("changed source measurement/error composition")
            k=point["kinematics"]
            Set(keys(k))==Set(("x","y","z","Q2_GeV2","pT_GeV")) && all(sidis_source_number,values(k)) &&
                all(0<k[name]<1 for name in ("x","y","z")) && k["Q2_GeV2"]>0 && k["pT_GeV"]>=0 &&
                (view["projection"]!="pt" || k["pT_GeV"]<.8) || error("invalid fixed means or changed reported-P selection")
            sidis_source_number(row["central_value"]) || error("invalid central value")
            if convention["id"]=="COMPASS2014"
                row["statistical_error"]===nothing && isempty(row["systematic_components"]) &&
                    row["uncertainty_status"]=="reported_error_u_composition_unresolved" &&
                    sidis_source_number(row["reported_error_u"]) && row["reported_error_u"]>0 || error("proton Sigma composition must stay unresolved")
            else
                sidis_source_number(row["statistical_error"]) && row["statistical_error"]>0 &&
                    Set(keys(row["systematic_components"]))==(convention["id"]=="HERMES2010" ? Set(["experimental"]) : Set{String}()) &&
                    all(v->sidis_source_number(v) && v>=0,values(row["systematic_components"])) || error("changed separate source errors")
            end
        end
    end
    length(unique(expected_ids))==132 && excluded==selection["excluded_points"]==14 &&
        132+excluded==selection["raw_points"]==146 &&
        length(unique(historical_sidis_geometry_key.(points)))==contract["prediction"]["distinct_source_tuples"]==132 ||
        error("incomplete full SIDIS identities, exclusions or exact geometry inventory")
    nothing
end

"Declared conditional target; sum the p/n structure functions before division."
function historical_sidis_target(name::AbstractString)
    name=="proton" && return NucleonTarget(:proton)
    name=="deuteron" || error("unsupported historical SIDIS target")
    NucleonTarget("deuteron",(.5,.5),(.5,.5);
        source="conditional unbound isoscalar p+n per-nucleon sum; no D-state, binding or nuclear admission")
end

function _historical_sidis_observable(request,point,context)
    k=point["kinematics"];plan=context.leptonic;target=request.data["target_weights"][point["target"]]
    plan.process==:SIDIS && plan.Q2==k["Q2_GeV2"] && plan.P==k["pT_GeV"] && plan.z1==k["z"] &&
        collect(context.first_operator.queries)==[k["x"]] && collect(context.second_operator.queries)==[k["z"]] &&
        context.target.name==point["target"] && collect(context.target.unpolarized)==target["unpolarized"] &&
        collect(context.target.polarized)==target["polarized"] || error("SIDIS context exact geometry/target mismatch")
    c=point["model_convention"];y=c["normalization"]=="fixed_y_DNN_included" ? k["y"] : nothing
    convention=SIDISAsymmetryConvention(phase=Symbol(c["phase"]),normalization=Symbol(c["normalization"]),y=y,
        source=c["source_binding_kind"]*"; conditional source mapping; request SHA256="*request.sha256)
    sidis_observable(with_pion_charges(context,point["charge"]),convention)
end

"One identified source point for a complete Julia pilot/shard job, never full coverage."
function historical_sidis_point_observable(request::HistoricalSIDISPointRequest,point_id::AbstractString,context)
    historical_sidis_request_valid(request)
    index=findfirst(p->p["id"]==point_id,request.data["points"])
    index===nothing && error("point is outside the full132 source request")
    _historical_sidis_observable(request,request.data["points"][index],context)
end

"Compose all132 predictions using existing Julia plans, without measurement/error arithmetic."
function historical_sidis_point_model(request::HistoricalSIDISPointRequest,contexts)
    historical_sidis_request_valid(request);observations=JointObservation[]
    for point in request.data["points"]
        key=historical_sidis_geometry_key(point)
        haskey(contexts,key) || error("missing exact full SIDIS geometry")
        observable=_historical_sidis_observable(request,point,contexts[key])
        push!(observations,JointObservation(point["id"],observable;
            source=point["publication"]*" conditional fixed-mean diagnostic; no experimental admission"))
    end
    JointPredictionModel(observations;source="all132 source-bound historical SIDIS points; request SHA256="*request.sha256)
end
