"Source-bound fixed-parameter BaBar point diagnostics; no likelihood or data admission."
const HISTORICAL_SIA_POINT_SCOPE="historical_reconstructed_sia_point_prediction_diagnostic"
const HISTORICAL_SIA_POINT_CONTRACT=normpath(joinpath(@__DIR__,"../../config/historical-sia-point-contract.json"))
const HISTORICAL_SIA_MODEL_CONVENTION=Dict(
    "frame"=>"hadron2","harmonic"=>"cos_2phi0","double_ratio_approximation"=>"linearized",
    "angular_definition"=>"mean_of_angular_ratio","kinematic_prescription"=>"fixed_reported_means",
    "legacy_Q2_status"=>"source_literal_not_checked_by_selected_primary_table",
    "source_mapping_status"=>"conditional_fixed_mean_factorized_angular_model_not_acceptance_admission")
const HISTORICAL_SIA_FALSE_FLAGS=("production_ready","physical_extraction","historical_fit_ready",
    "experimental_data_admitted","fit_performed","G2PERF_pass","interval_coverage_verified")

struct HistoricalSIAPointRequest
    data::Dict{String,Any}
    path::String
    sha256::String
    signature::String
end

historical_sia_geometry_key(point)=Tuple(point["kinematics"][k] for k in ("z1","z2","Q2_GeV2","pT_GeV"))

function historical_sia_request_valid(request::HistoricalSIAPointRequest)
    numerical_snapshot(request.data)==request.signature || error("mutated historical SIA request")
    bytes2hex(sha256(read(request.path)))==request.sha256 || error("changed historical SIA request bytes")
    for (path,digest) in request.data["source_inputs"]
        isfile(path) && bytes2hex(sha256(read(path)))==digest || error("changed historical SIA source input: $path")
    end
    request
end

function historical_sia_source_bindings()
    root=normpath(joinpath(@__DIR__,"../.."))
    paths=[joinpath(root,p) for p in ("config/historical-sia-point-contract.json","config/initial-contract.json",
        "config/historical-sources.json","fortran/hoppet_bridge.f90","bridges/lhapdf.cpp",
        "numerics/Project.toml","numerics/Manifest.toml","requirements.lock")]
    for (directory,suffix) in (("numerics/src",".jl"),("numerics/bin",".jl"),("src/collinsflow",".py"),("config","-contract.json"))
        append!(paths,filter(p->endswith(p,suffix),readdir(joinpath(root,directory);join=true)))
    end
    Dict(replace(relpath(p,root),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in unique(paths))
end

function load_historical_sia_point_request(path)
    data=JSON.parsefile(path);contract=JSON.parsefile(HISTORICAL_SIA_POINT_CONTRACT)
    initial=JSON.parsefile(HISTORICAL_INITIAL_CONTRACT)
    data isa AbstractDict && get(data,"schema_version",nothing) isa Integer && !(data["schema_version"] isa Bool) &&
        data["schema_version"]==1 && get(data,"scope",nothing)==HISTORICAL_SIA_POINT_SCOPE &&
        get(data,"theory_scope",nothing)=="reconstructed_theory" &&
        get(data,"observation_scope",nothing)=="historical_benchmark_candidate" || error("wrong historical SIA request scope")
    for key in HISTORICAL_SIA_FALSE_FLAGS
        get(data,key,nothing)===false && get(contract,key,nothing)===false || error("SIA diagnostic cannot approve fitting or production")
    end
    digest(p)=bytes2hex(sha256(read(p)))
    data["contract_id"]==contract["id"]=="KPSY-HISTORICAL-SIA-POINT-001" &&
        data["contract_sha256"]==digest(HISTORICAL_SIA_POINT_CONTRACT) &&
        data["initial_contract_sha256"]==digest(HISTORICAL_INITIAL_CONTRACT) || error("stale historical SIA contract")
    parameters=data["parameters"]
    parameters["order"]==initial["parameter_order"] && parameters["values"]==initial["published_parameters"] &&
        parameters["values"]==collect(kpsy_published_parameters().values) &&
        parameters["role"]==initial["published_parameter_role"] &&
        all(v->v isa Real && !(v isa Bool) && isfinite(v),parameters["values"]) || error("changed fixed SIA parameter source")
    data["model_convention"]==HISTORICAL_SIA_MODEL_CONVENTION &&
        data["limitations"]==contract["limitations"] || error("changed conditional SIA convention or limitations")
    source=data["source_preparation"];inputs=data["source_inputs"]
    source["status"]=="CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED" && source["scope"]=="historical_benchmark_candidate" &&
        source["verification_kind"]=="Python_source_preparation_attestation_not_Julia_raw_table_replay" ||
        error("source preparation attestation required")
    candidate_path=joinpath(source["candidates"],"candidates.json");report_path=joinpath(source["candidates"],"report.json")
    candidate=JSON.parsefile(candidate_path);report=JSON.parsefile(report_path)
    candidate["scope"]=="historical_benchmark_candidate" && report["status"]=="CANDIDATES_VALIDATED" &&
        candidate["source_fingerprint"]==report["source_fingerprint"]==data["source_fingerprint"] || error("stale SIA candidate source")
    for obj in (candidate,report),key in ("production_ready","physical_extraction","historical_fit_ready")
        get(obj,key,nothing)===false || error("historical source cannot self-admit")
    end
    source["historical_counts"]==report["counts"]==Dict("views"=>26,"points"=>262,"SIDIS_views"=>20,
        "SIDIS_points"=>140,"SIA_views"=>6,"SIA_points"=>122) &&
        report["primary_BaBar_raw_points_crosschecked"]==90 && report["primary_BaBar_numeric_cells_crosschecked"]==954 ||
        error("incomplete historical/BaBar source replay")
    for (process,nviews,npoints) in (("SIDIS",20,140),("SIA",6,122))
        selected_views=filter(v->v["process"]==process,candidate["views"])
        length(selected_views)==nviews && sum(length(v["points"]) for v in selected_views)==npoints || error("candidate counts contradict source attestation")
    end
    manifest=normpath(joinpath(@__DIR__,"../../config/historical-sources.json"))
    source["manifest"]==manifest && get(report["source_hashes"],manifest,nothing)==digest(manifest) || error("missing current historical manifest binding")
    expected_inputs=merge(Dict{String,Any}(report["source_hashes"]),Dict(candidate_path=>digest(candidate_path),
        report_path=>digest(report_path),HISTORICAL_SIA_POINT_CONTRACT=>digest(HISTORICAL_SIA_POINT_CONTRACT),
        HISTORICAL_INITIAL_CONTRACT=>digest(HISTORICAL_INITIAL_CONTRACT)))
    inputs==expected_inputs && inputs[candidate_path]==source["candidate_sha256"] &&
        inputs[report_path]==source["report_sha256"] || error("incomplete or changed source-input bindings")
    views=Dict(v["id"]=>v for v in candidate["views"])
    length(views)==length(candidate["views"]) || error("duplicate SIA candidate view")
    points=data["points"];selection=contract["selection"];primary=contract["primary_definition"]
    length(points)==selection["total_points"]==18 || error("incomplete BaBar point request")
    ids=String[]
    finite(v)=v isa Real && !(v isa Bool) && isfinite(v)
    for selected in selection["source_views"]
        view=get(views,selected["id"],nothing)
        view!==nothing && view["source_id"]==selected["id"] && view["source_sha256"]==selected["sha256"] &&
            view["process"]=="SIA" && view["projection"]=="pT" && view["ratio"]==selected["ratio"] &&
            selected["ratio"] in ("UL","UC") && view["frame"]=="hadron2" && view["publication"]=="1309.5278" &&
            view["source_observable"]=="A"*selected["ratio"]*"-0-PT" && view["input_units"]=="percent" &&
            view["output_units"]=="fraction" && view["overlap_group"]=="babar_1309.5278" &&
            view["production_ready"]===false && view["historical_fit_ready"]===false && !isempty(view["blockers"]) &&
            length(view["points"])==selected["points"] && isempty(view["excluded"]) || error("mismatched primary RF0 selection")
        cross=view["primary_numeric_crosscheck"]
        cross["publication"]==selection["publication"] && cross["frame"]=="hadron2" &&
            cross["sha256"]==primary["archive_sha256"] && cross["member"]==primary["member"] &&
            cross["member_sha256"]==primary["member_sha256"] && cross["table"]=="tab:totPt" &&
            cross["raw_points"]==selected["points"] && cross["numeric_cells"]==9*selected["points"] &&
            Set(cross["compared_columns"])==Set(primary["checked_columns"]) || error("incomplete primary P-table identity")
        for (offset,row) in enumerate(view["points"])
            row["id"]==selected["id"]*":r"*string(offset+1) && row["locator"]["excel_row"]==offset+1 || error("changed source point order")
            push!(ids,row["id"])
            expected=Dict{String,Any}("id"=>row["id"],"dataset_id"=>view["id"],"publication"=>view["publication"],
                "primary_publication"=>cross["publication"],"source_sha256"=>view["source_sha256"],
                "ratio"=>view["ratio"],"frame"=>view["frame"],"harmonic"=>"cos_2phi0","kinematics"=>row["averages"],
                "measured"=>Dict(k=>row[k] for k in ("central_value","statistical_error","systematic_components")),
                "source_locator"=>row["locator"],"primary_locator"=>row["primary_locator"],
                "source_bins"=>row["bins"],"primary_bins"=>row["primary_bins"],"source_angular"=>row["angular_source"],
                "source_primary_crosscheck"=>cross,"error_source_columns"=>row["uncertainty_source_columns"],
                "overlap_group"=>view["overlap_group"],"source_blockers"=>view["blockers"],
                "legacy_Q2_status"=>HISTORICAL_SIA_MODEL_CONVENTION["legacy_Q2_status"])
            point=points[length(ids)];point==expected || error("changed source-derived SIA point")
            k=point["kinematics"];m=point["measured"];a=point["source_angular"]
            Set(keys(k))==Set(("z1","z2","Q2_GeV2","pT_GeV")) && all(finite,values(k)) &&
                all(0<k[key]<1 for key in ("z1","z2")) && k["Q2_GeV2"]==112.36 && k["pT_GeV"]>=0 || error("invalid SIA reported means")
            finite(m["central_value"]) && finite(m["statistical_error"]) && m["statistical_error"]>0 &&
                Set(keys(m["systematic_components"]))==Set(("reported_systematic",)) &&
                all(v->finite(v) && v>=0,values(m["systematic_components"])) &&
                point["error_source_columns"]==Dict("statistical"=>"stat_u","reported_systematic"=>"sys_u") || error("invalid source error components")
            a["column"]=="S2/1+C2" && a["average_type"]=="mean_of_angular_ratio" && a["frame"]=="hadron2" &&
                a["measure_status"]=="joint_acceptance_unresolved" && !isempty(a["averaging_measure"]) &&
                finite(a["value"]) && 0<=a["value"]<=1 || error("unreviewed angular source")
            loc=point["primary_locator"];bins=point["primary_bins"]
            loc["publication"]==selection["publication"] && loc["member"]==primary["member"] && loc["table"]=="tab:totPt" &&
                loc["table_row"]==offset && loc["source_line"] isa Integer && !(loc["source_line"] isa Bool) && loc["source_line"]>0 &&
                Set(keys(bins))==Set(("pT",)) && bins["pT"]==point["source_bins"]["pT_GeV"] && length(bins["pT"])==2 &&
                all(v->finite(v) && v>=0,bins["pT"]) && bins["pT"][1]<bins["pT"][2] &&
                bins["pT"][1]<=k["pT_GeV"]<=bins["pT"][2] || error("invalid primary source point/bin locator")
        end
    end
    length(unique(ids))==length(ids) && length(unique(historical_sia_geometry_key.(points)))==selection["exact_shared_geometry_count"]==9 ||
        error("duplicate identities or mismatched exact SIA geometries")
    owned=deepcopy(Dict{String,Any}(data))
    historical_sia_request_valid(HistoricalSIAPointRequest(owned,abspath(path),digest(path),numerical_snapshot(owned)))
end

"Build shared thirteen-coordinate prediction plans; measured data never enter prediction."
function historical_sia_point_model(request::HistoricalSIAPointRequest,contexts)
    historical_sia_request_valid(request)
    observations=JointObservation[]
    for point in request.data["points"]
        key=historical_sia_geometry_key(point);haskey(contexts,key) || error("missing exact SIA geometry context")
        context=contexts[key];plan=context.leptonic;k=point["kinematics"]
        plan.process==:SIA && plan.Q2==k["Q2_GeV2"] && plan.P==k["pT_GeV"] &&
            plan.z1==k["z1"] && plan.z2==k["z2"] && collect(context.first_operator.queries)==[k["z1"],k["z2"]] &&
            context.second_operator===nothing || error("SIA context geometry or operator mismatch")
        a=point["source_angular"]
        angular=SIAMeanOfRatio(a["value"];source="BaBar RF0 primary table; request SHA256="*request.sha256,
            averaging_measure=a["averaging_measure"])
        observable=sia_observable(context;ratio=Symbol(point["ratio"]),angular=angular,
            source="conditional fixed-mean BaBar RF0 diagnostic, no experimental admission")
        push!(observations,JointObservation(point["id"],observable;source=point["publication"]*" RF0 source row"))
    end
    JointPredictionModel(observations;source="source-bound historical SIA point diagnostic; request SHA256="*request.sha256)
end
