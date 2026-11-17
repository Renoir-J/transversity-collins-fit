"Native value computations and independent leaf arithmetic; no physical dispatcher."
module SIAGeometryNativeJobs
using CollinsFit, JSON, SHA
include("SIAGeometryModel.jl")
const M=SIAGeometryModel
const S=M.S
const FLAGS=(S.FLAGS...,"accepted_fitting_gate","accepted_prediction_diagnostic","numerical_dispatch_ready")
const COMPUTED="GEOMETRY_NATIVE_VALUES_COMPUTED_NOT_DIAGNOSTIC_GATE"
const REVIEWED="GEOMETRY_NATIVE_SAVED_ARITHMETIC_REVIEWED_NOT_COMPLETE_JOB_OR_GATE"
flags()=Dict{String,Any}(key=>false for key in FLAGS)
require(value,message)=value || error(message)
digest(path)=bytes2hex(sha256(read(path)))

function definition(data,parent,phase,index)
    S.validate_parent(parent);groups,keys=S.geometry_groups(data)
    require(phase in ("native123","native4") && index isa Int && 1<=index<=36,"undeclared native phase or geometry")
    Dict{String,Any}("scope"=>"mathematical_job_shape_not_dispatch_authorization","phase"=>phase,
        "geometry_index"=>index,"geometry"=>collect(keys[index]),"source_point_indices"=>groups[index],
        "point_ids"=>[data["points"][i]["id"] for i in groups[index]],
        "native_levels"=>phase=="native123" ? [1,2,3] : [4],"bound_inverse_GeV"=>24,
        "momentum_orders"=>copy(parent["numerical_policy"]["momentum_orders"]))
end

function literal_theta(data)
    theta=data["parameters"]["values"]
    require(length(theta)==13 && all(S.finite_number,theta) &&
        Float64.(theta)==collect(kpsy_published_parameters().values),"literal rounded13 parameters required")
    Float64.(theta)
end

"Verify retained builds once per complete job; no new response cache is added."
function native_inputs(nativefile,inputfile)
    root=dirname(dirname(dirname(pathof(CollinsFit))))
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    require(any(basename(p)=="hoppet_bridge.f90" && h==digest(joinpath(root,"fortran/hoppet_bridge.f90"))
        for (p,h) in native["source_hashes"]),"stale native bridge")
    require(any(basename(p)=="lhapdf.cpp" && h==digest(joinpath(root,"bridges/lhapdf.cpp"))
        for (p,h) in inputs["source_hashes"]),"stale parton bridge")
    (nativefile=abspath(nativefile),inputfile=abspath(inputfile),native=native,inputs=inputs,
        native_sha256=digest(nativefile),input_sha256=digest(inputfile),ff=load_dss14(inputs["input_root"]),css=KPSYCSS())
end

function inputs_stable(inputs)
    require(digest(inputs.nativefile)==inputs.native_sha256 && digest(inputs.inputfile)==inputs.input_sha256 &&
        CollinsFit.verified_hoppet_build(inputs.nativefile)==inputs.native &&
        CollinsFit.verify_parton_build(inputs.inputfile)==inputs.inputs,"native or parton build changed during job")
end

"The parent's exact preparation for one geometry/native setting at B24."
function prepare_native(view,parent,inputs,level)
    M.valid(view);S.validate_parent(parent)
    require(level isa Int && level in 1:4,"undeclared native level")
    setting=parent["numerical_policy"]["refinements"][level];k=first(view.data["points"])["kinematics"]
    radial=RadialRule(parton_radial_edges(inputs.css,inputs.ff,inputs.inputs["input_root"],24);order=setting["radial_order"])
    context=nothing
    elapsed=@elapsed context=prepare_polarized_sia(inputs.ff,inputs.nativefile,inputs.css,radial,k["Q2_GeV2"],0.,k["z1"],k["z2"],1,1;
        matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference)
    require(context.first_samples.input_kind==:pinned_inputs && context.first_samples.kind==:KPSY_Hhat3_homogeneous &&
        context.second_samples===nothing,"pinned fragmentation-only initial samples required")
    resource=Dict{String,Any}("geometry"=>collect(historical_sia_integral_geometry_key(first(view.data["points"]))),
        "radial_nodes"=>length(radial.b),"row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes,
        "cold_seconds_informational"=>elapsed,"context_provenance"=>deepcopy(context.provenance))
    context,resource
end

"Julia callbacks select fixed test or physical preparation outside prediction loops."
function compute(view,parent,spec,binding,legacy,prepare;progress=stage->nothing,checkpoint=report->nothing)
    M.valid(view)
    expected=definition(view.full.data,parent,spec["phase"],view.index)
    require(S.same_json(spec,expected),"changed native mathematical job definition")
    theta=literal_theta(view.full.data);theory=S.declaration(view.full.data)
    report=merge(flags(),Dict{String,Any}("schema_version"=>1,"status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>S.SCOPE,
        "theory_scope"=>theory,"source_context"=>deepcopy(binding),"definition"=>deepcopy(spec),
        "request_points"=>deepcopy(view.data["points"]),"parameters"=>theta,"native_values"=>Any[],
        "native_evolution_reexecuted_by_independent_reader"=>false,
        "producer_attestation_scope"=>"Temporal preparation/object identity attestations are not independently repeated native evolution."))
    checkpoint(report)
    for level in spec["native_levels"]
        progress(Dict("native_level"=>level,"operation"=>"prepare one fixed geometry"))
        context,resource=prepare(view,parent,level)
        require(context.first_samples.input_kind==(theory=="synthetic_closure" ? :analytic : :pinned_inputs),
            "synthetic rows and pinned physical inputs must remain distinct")
        key=historical_sia_integral_geometry_key(first(view.data["points"]))
        value=merge(flags(),Dict{String,Any}("status"=>S.VALUES_STATUS,"scope"=>S.SCOPE,"theory_scope"=>theory,
            "geometry_index"=>view.index,"geometry"=>collect(key),"native_level"=>level,"bound_inverse_GeV"=>24,
            "setting"=>deepcopy(parent["numerical_policy"]["refinements"][level]),"momentum_levels"=>Any[],
            "resources"=>[deepcopy(resource)],"initial_source_kind"=>string(context.first_samples.input_kind),
            "same_native_rows_across_ratios_producer_attestation"=>true))
        push!(report["native_values"],value)
        checkpoint(report)
        for order in spec["momentum_orders"]
            progress(Dict("native_level"=>level,"momentum_order"=>order,"operation"=>"evaluate both source ratios"))
            result=legacy.predictions(view,Dict(key=>context),theta,order)
            require(all(observation.plan.context.first_operator===context.first_operator && observation.plan.context.second_operator===nothing
                for observation in result.model.observations),"expected unchanged native row objects")
            for (row,point) in zip(result.points,view.data["points"])
                S.validate_point(row,point,order,parent["numerical_policy"])
            end
            require(S.same_json(result.points[1]["pairs"],result.points[2]["pairs"]),"both source ratios must retain the same absolute pair integrals")
            push!(value["momentum_levels"],Dict("order"=>order,"points"=>result.points))
            checkpoint(report)
        end
    end
    M.valid(view);report["status"]=COMPUTED;checkpoint(report)
    report
end

function audit(data,parent,spec,binding,report)
    expected=definition(data,parent,spec["phase"],spec["geometry_index"])
    require(S.same_json(spec,expected),"changed native mathematical job definition")
    theory=S.declaration(data);points=data["points"][spec["source_point_indices"]]
    require(report["schema_version"]===1 && report["status"]==COMPUTED && report["scope"]==S.SCOPE &&
        report["theory_scope"]==theory && all(get(report,key,nothing)===false for key in FLAGS) &&
        report["native_evolution_reexecuted_by_independent_reader"]===false,"unaccepted or self-admitting native values")
    require(S.same_json(report["source_context"],binding) && S.same_json(report["definition"],spec) &&
        S.same_json(report["request_points"],points) && S.same_json(report["parameters"],literal_theta(data)),"changed source/native job identity")
    values=report["native_values"]
    require(length(values)==length(spec["native_levels"]),"all declared native levels required")
    count=0
    for (record,level) in zip(values,spec["native_levels"])
        require(record["status"]==S.VALUES_STATUS && record["scope"]==S.SCOPE && record["theory_scope"]==theory &&
            all(get(record,key,nothing)===false for key in FLAGS),"wrong native value scope or admission")
        require(record["geometry_index"]===spec["geometry_index"] && S.same_json(record["geometry"],spec["geometry"]) &&
            record["native_level"]===level && record["bound_inverse_GeV"]===24 &&
            S.same_json(record["setting"],parent["numerical_policy"]["refinements"][level]),"changed geometry/native setting")
        require(record["initial_source_kind"]==(theory=="synthetic_closure" ? "analytic" : "pinned_inputs") &&
            record["same_native_rows_across_ratios_producer_attestation"]===true,"changed source kind or preparation attestation")
        resources=record["resources"]
        require(length(resources)==1 && S.same_json(resources[1]["geometry"],spec["geometry"]) &&
            all(resources[1][key] isa Int && resources[1][key]>0 for key in ("radial_nodes","row_coefficient_payload_bytes")) &&
            S.finite_number(resources[1]["cold_seconds_informational"]) && resources[1]["cold_seconds_informational"]>=0 &&
            resources[1]["context_provenance"] isa AbstractDict,"incomplete preparation observations")
        orders=record["momentum_levels"]
        require(length(orders)==5 && S.same_json([row["order"] for row in orders],spec["momentum_orders"]),"all five momentum orders required")
        for (row,order) in zip(orders,spec["momentum_orders"])
            require(length(row["points"])==2,"both ordered source ratios required")
            for (point,source) in zip(row["points"],points)
                S.validate_point(point,source,order,parent["numerical_policy"]);count+=1
            end
            require(S.same_json(row["points"][1]["pairs"],row["points"][2]["pairs"]),"absolute pair yields changed between source ratios")
        end
    end
    merge(flags(),Dict{String,Any}("schema_version"=>1,"status"=>REVIEWED,"scope"=>S.SCOPE,"theory_scope"=>theory,
        "source_context"=>deepcopy(binding),"definition"=>deepcopy(spec),"saved_value_points_recomputed"=>count,
        "native_evolution_reexecuted_by_independent_reader"=>false,"global_common_order_or_native_selection_performed"=>false,
        "job_resource_and_complete_provenance_ledger_verified"=>false))
end
end
