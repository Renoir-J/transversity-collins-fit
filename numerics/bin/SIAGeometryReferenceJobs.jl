"Finite-tail and representative derivative operands; no complete-job dispatcher."
module SIAGeometryReferenceJobs
using CollinsFit, JSON
include("SIAGeometryNativeJobs.jl")
const N=SIAGeometryNativeJobs
const M=N.M
const S=N.S
const TAIL_COMPUTED="GEOMETRY_TAILREF_OPERANDS_COMPUTED_NOT_DIAGNOSTIC_GATE"
const TAIL_REVIEWED="GEOMETRY_TAILREF_SAVED_ARITHMETIC_REVIEWED_NOT_COMPLETE_JOB_OR_GATE"
const DERIVATIVE_COMPUTED="GEOMETRY_DERIVATIVE_OPERANDS_COMPUTED_NOT_DIAGNOSTIC_GATE"
const DERIVATIVE_REVIEWED="GEOMETRY_DERIVATIVE_SAVED_ARITHMETIC_REVIEWED_NOT_COMPLETE_JOB_OR_GATE"
const REFERENCE_SCOPE="generic_fixed_row_scalar_assembly_and_adaptive_qT_reference"
require(value,message)=value || error(message)

function tail_definition(data,parent,index,level)
    spec=N.definition(data,parent,"native123",index)
    require(level isa Int && level in (3,4),"one globally selected native3 or native4 is required")
    Dict{String,Any}("scope"=>"mathematical_job_shape_not_dispatch_authorization","phase"=>"tailref",
        "geometry_index"=>index,"geometry"=>spec["geometry"],"point_ids"=>spec["point_ids"],
        "source_point_indices"=>spec["source_point_indices"],"selected_native_level"=>level,
        "bounds_inverse_GeV"=>[12,32],"momentum_orders"=>spec["momentum_orders"],
        "B24_policy"=>"Exact selected native outputs are bound by the complete global ledger, never recomputed or relabeled here.")
end

function derivative_definition(data,parent,index,level,order)
    spec=tail_definition(data,parent,index,level)
    source_id=first(spec["point_ids"])
    require(source_id in parent["representative_derivative_policy"]["source_geometry_ids"] &&
        order isa Int && order in parent["numerical_policy"]["momentum_orders"] &&
        order>=parent["numerical_policy"]["minimum_momentum_order"],"declared representative geometry and common momentum order required")
    Dict{String,Any}("scope"=>"mathematical_job_shape_not_dispatch_authorization","phase"=>"derivative",
        "geometry_index"=>index,"geometry"=>spec["geometry"],"source_point_indices"=>spec["source_point_indices"],
        "point_ids"=>spec["point_ids"],"geometry_source_id"=>source_id,"selected_native_level"=>level,
        "bound_inverse_GeV"=>32,"selected_momentum_order"=>order)
end

function report_base(view,spec,binding)
    merge(N.flags(),Dict{String,Any}("schema_version"=>1,"status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>S.SCOPE,
        "theory_scope"=>S.declaration(view.full.data),"source_context"=>deepcopy(binding),"definition"=>deepcopy(spec),
        "request_points"=>deepcopy(view.data["points"]),"parameters"=>N.literal_theta(view.full.data),
        "global_common_order_or_native_selection_performed"=>false,"native_evolution_reexecuted_by_independent_reader"=>false))
end

function shared_reference(legacy,context,angular)
    legacy.HSP.ObservableReference.prepare_generic_reference(context;momentum_interval=(0.,3.5),momentum_coordinate=:qT,angular=angular)
end

function compute_tail(view,parent,spec,binding,legacy,prepare;progress=stage->nothing,checkpoint=report->nothing)
    M.valid(view)
    require(S.same_json(spec,tail_definition(view.full.data,parent,view.index,spec["selected_native_level"])),"changed tail mathematical job definition")
    report=report_base(view,spec,binding);theta=report["parameters"];floor=parent["numerical_policy"]["relative_floor"]
    report["tail_values"]=Any[];report["independent_references"]=Any[];report["fixed_value_controls"]=Any[];checkpoint(report)
    for B in spec["bounds_inverse_GeV"]
        progress(Dict("bound_inverse_GeV"=>B,"operation"=>"prepare one fixed geometry"))
        context,resource=prepare(view,parent,spec["selected_native_level"],B)
        value=N.value_record(view,parent,spec["selected_native_level"],B,context,resource)
        push!(report["tail_values"],value);checkpoint(report)
        models=N.fill_values!(value,view,parent,legacy,context,theta;progress=progress,checkpoint=()->checkpoint(report))
        if B==32
            for (i,point) in enumerate(view.data["points"])
                progress(Dict("point_id"=>point["id"],"operation"=>"independent adaptive qT reference"))
                angular=point["source_angular"]["value"];reference=shared_reference(legacy,context,angular)
                pol=reference.polarized_integrals(theta);den=reference.denominators
                pairs=Dict(name=>Dict("UU"=>den[j],"Collins"=>pol[j]) for (j,name) in enumerate(S.PAIRS))
                ref=Dict("pairs"=>pairs,"categories"=>legacy.HSP.category_values(pairs),
                    "amplitude"=>reference.predict(theta)[point["ratio"]=="UL" ? 1 : 2])
                comparisons=[Dict("order"=>row["order"],"max_scaled"=>S.scaled(S.packed_point(row["points"][i]),S.packed_point(ref),floor))
                    for row in value["momentum_levels"]]
                push!(report["independent_references"],Dict("id"=>point["id"],"reference"=>ref,"angular_value"=>angular,
                    "scope"=>reference.scope,"qT_bounds_GeV"=>[0.,3.5],"native_rows_shared"=>true,
                    "initial_functions_shared"=>true,"independent_evolution"=>false,"order_comparisons"=>comparisons))
                checkpoint(report)
            end
            flipped=copy(theta);flipped[[7,8]].*=-1
            changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
            for row in value["momentum_levels"]
                order=row["order"];original=[point["amplitude"] for point in row["points"]]
                flip_values=predict_joint(models[order],flipped);changed_values=predict_joint(models[order],changed)
                push!(report["fixed_value_controls"],Dict("order"=>order,"theta"=>theta,"flipped_theta"=>flipped,
                    "changed_transversity_theta"=>changed,"original"=>original,"flipped"=>flip_values,
                    "changed_transversity"=>changed_values,"flipped_max_scaled"=>S.scaled(flip_values,original,floor),
                    "transversity_max_scaled"=>S.scaled(changed_values,original,floor)))
                checkpoint(report)
            end
        end
    end
    M.valid(view);report["status"]=TAIL_COMPUTED;checkpoint(report);report
end

function compute_derivative(view,parent,spec,binding,legacy,prepare;progress=stage->nothing,checkpoint=report->nothing)
    M.valid(view)
    require(S.same_json(spec,derivative_definition(view.full.data,parent,view.index,spec["selected_native_level"],spec["selected_momentum_order"])),
        "changed derivative mathematical job definition")
    report=report_base(view,spec,binding);checkpoint(report)
    context,resource=prepare(view,parent,spec["selected_native_level"],32)
    require(context.first_samples.input_kind==(report["theory_scope"]=="synthetic_closure" ? :analytic : :pinned_inputs),"changed derivative source kind")
    report["resources"]=[resource];report["initial_source_kind"]=string(context.first_samples.input_kind);checkpoint(report)
    theta=report["parameters"];key=historical_sia_integral_geometry_key(first(view.data["points"]))
    result=legacy.predictions(view,Dict(key=>context),theta,spec["selected_momentum_order"])
    report["value_points"]=result.points;checkpoint(report)
    indices=[only(findall(point->point["ratio"]==ratio,view.data["points"])) for ratio in ("UL","UC")]
    plans=[result.model.observations[i].plan for i in indices]
    angulars=[view.data["points"][i]["source_angular"]["value"] for i in indices]
    reference=shared_reference(legacy,context,1.)
    candidate=t->[predict_observable(plan,t) for plan in plans];independent=t->angulars.*reference.predict(t)
    record=legacy.HIR.derivative_reference(candidate,independent,theta,parent["representative_derivative_policy"];
        progress=()->progress(Dict("geometry_source_id"=>spec["geometry_source_id"],"operation"=>"inward all13 reference")))
    record["geometry_source_id"]=spec["geometry_source_id"];record["point_ids"]=[view.data["points"][i]["id"] for i in indices]
    record["angular_values"]=angulars;record["native_rows_shared"]=true;record["independent_evolution"]=false
    report["derivative_record"]=record;M.valid(view);report["status"]=DERIVATIVE_COMPUTED;checkpoint(report);report
end

function audit_identity(data,parent,spec,binding,report,status)
    S.validate_parent(parent);groups,_=S.geometry_groups(data)
    require(report["schema_version"]===1 && report["status"]==status && report["scope"]==S.SCOPE &&
        report["theory_scope"]==S.declaration(data) && all(get(report,key,nothing)===false for key in N.FLAGS) &&
        report["global_common_order_or_native_selection_performed"]===false && report["native_evolution_reexecuted_by_independent_reader"]===false,
        "unaccepted or self-admitting reference job operands")
    points=data["points"][groups[spec["geometry_index"]]]
    require(S.same_json(report["source_context"],binding) && S.same_json(report["definition"],spec) &&
        S.same_json(report["request_points"],points) && S.same_json(report["parameters"],N.literal_theta(data)),"changed reference job source identity")
    points
end

"A scalar adaptive reference has no Gauss momentum-order label."
function validate_reference(ref,point,policy)
    require(Set(keys(ref["pairs"]))==Set(S.PAIRS) && Set(keys(ref["categories"]))==Set(("U","L","C")),"incomplete scalar pair/category operands")
    for value in vcat(collect(values(ref["pairs"])),collect(values(ref["categories"])))
        require(Set(keys(value))==Set(("UU","Collins")) && all(S.finite_number,values(value)) && value["UU"]>0,"invalid scalar signed yields or denominator")
    end
    for key in ("UU","Collins")
        u=ref["pairs"]["pm"][key]+ref["pairs"]["mp"][key]
        l=ref["pairs"]["pp"][key]+ref["pairs"]["mm"][key]
        require(ref["categories"]["U"][key]==u && ref["categories"]["L"][key]==l && ref["categories"]["C"][key]==u+l,
            "changed scalar category sum")
    end
    u=ref["categories"]["U"];right=ref["categories"][point["ratio"]=="UL" ? "L" : "C"]
    algebra=point["source_angular"]["value"]*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"])
    require(S.finite_number(ref["amplitude"]) && S.scaled([ref["amplitude"]],[algebra],policy["relative_floor"])<=policy["observable_algebra_max_scaled"],
        "invalid scalar reference ratio algebra")
    ref
end

function audit_tail(data,parent,spec,binding,report)
    require(S.same_json(spec,tail_definition(data,parent,spec["geometry_index"],spec["selected_native_level"])),"changed tail mathematical job definition")
    points=audit_identity(data,parent,spec,binding,report,TAIL_COMPUTED);policy=parent["numerical_policy"];floor=policy["relative_floor"]
    values=report["tail_values"];require(length(values)==2,"both finite bounds required; B24 belongs to the exact native ledger")
    count=sum(N.audit_value(data,parent,value,spec["geometry_index"],spec["selected_native_level"],B) for (value,B) in zip(values,(12,32)))
    references=report["independent_references"];require(length(references)==2,"both source scalar references required")
    for (i,(row,point)) in enumerate(zip(references,points))
        require(row["id"]==point["id"] && row["angular_value"]===point["source_angular"]["value"] && row["scope"]==REFERENCE_SCOPE &&
            S.same_json(row["qT_bounds_GeV"],[0.,3.5]) && row["native_rows_shared"]===true && row["initial_functions_shared"]===true &&
            row["independent_evolution"]===false,"changed scalar reference definition")
        ref=validate_reference(row["reference"],point,policy)
        expected=[Dict("order"=>entry["order"],"max_scaled"=>S.scaled(S.packed_point(entry["points"][i]),S.packed_point(ref),floor))
            for entry in values[2]["momentum_levels"]]
        require(S.same_json(row["order_comparisons"],expected),"changed saved adaptive reference comparisons")
    end
    require(S.same_json(references[1]["reference"]["pairs"],references[2]["reference"]["pairs"]),"absolute scalar pair yields changed between ratios")
    controls=report["fixed_value_controls"];require(length(controls)==5,"all fixed-value control orders required")
    theta=N.literal_theta(data);flipped=copy(theta);flipped[[7,8]].*=-1;changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
    for (control,value) in zip(controls,values[2]["momentum_levels"])
        original=[point["amplitude"] for point in value["points"]]
        require(control["order"]===value["order"] && S.same_json(control["theta"],theta) && S.same_json(control["flipped_theta"],flipped) &&
            S.same_json(control["changed_transversity_theta"],changed) && S.same_json(control["original"],original),"changed fixed-value operands")
        require(all(length(control[key])==2 && all(S.finite_number,control[key]) for key in ("flipped","changed_transversity")),"invalid fixed-value outputs")
        fe=S.scaled(control["flipped"],original,floor);te=S.scaled(control["changed_transversity"],original,floor)
        require(control["flipped_max_scaled"]==fe && control["transversity_max_scaled"]==te &&
            fe<=policy["observable_algebra_max_scaled"] && te<=policy["observable_algebra_max_scaled"],"failed or changed fixed-value comparison")
    end
    merge(N.flags(),Dict{String,Any}("schema_version"=>1,"status"=>TAIL_REVIEWED,"scope"=>S.SCOPE,"theory_scope"=>S.declaration(data),
        "source_context"=>deepcopy(binding),"definition"=>deepcopy(spec),"saved_value_points_recomputed"=>count,
        "saved_scalar_reference_comparisons_recomputed"=>10,"fixed_value_control_comparisons_recomputed"=>10,
        "global_common_order_or_native_selection_performed"=>false,"native_evolution_reexecuted_by_independent_reader"=>false,
        "job_resource_and_complete_provenance_ledger_verified"=>false,"infinite_tail_proved"=>false))
end

function audit_derivative(data,parent,spec,binding,report,inward_reader)
    require(S.same_json(spec,derivative_definition(data,parent,spec["geometry_index"],spec["selected_native_level"],spec["selected_momentum_order"])),
        "changed derivative mathematical job definition")
    points=audit_identity(data,parent,spec,binding,report,DERIVATIVE_COMPUTED)
    require(report["initial_source_kind"]==(S.declaration(data)=="synthetic_closure" ? "analytic" : "pinned_inputs"),"changed derivative input kind")
    resources=report["resources"]
    require(length(resources)==1 && S.same_json(resources[1]["geometry"],spec["geometry"]) && resources[1]["bound_inverse_GeV"]===32 &&
        all(resources[1][key] isa Int && resources[1][key]>0 for key in ("radial_nodes","row_coefficient_payload_bytes")) &&
        S.finite_number(resources[1]["cold_seconds_informational"]) && resources[1]["cold_seconds_informational"]>=0 &&
        resources[1]["context_provenance"] isa AbstractDict,"incomplete derivative preparation observations")
    require(length(report["value_points"])==2,"both derivative value points required")
    for (row,point) in zip(report["value_points"],points)
        S.validate_point(row,point,spec["selected_momentum_order"],parent["numerical_policy"])
    end
    require(S.same_json(report["value_points"][1]["pairs"],report["value_points"][2]["pairs"]),"derivative value pair yields differ between ratios")
    record=report["derivative_record"];ordered=[only(filter(p->p["ratio"]==ratio,points)) for ratio in ("UL","UC")]
    require(record["scope"]=="three_source_selected_geometries_inward_all13_reference_not_production_certificate" &&
        S.same_json(record["theta"],N.literal_theta(data)) && S.same_json(record["direction"],parent["representative_derivative_policy"]["direction"]) &&
        record["geometry_source_id"]==spec["geometry_source_id"] && record["point_ids"]==[p["id"] for p in ordered] &&
        S.same_json(record["angular_values"],[p["source_angular"]["value"] for p in ordered]) && record["native_rows_shared"]===true &&
        record["independent_evolution"]===false,"changed representative derivative source definition")
    checks=inward_reader(record,parent["representative_derivative_policy"])
    require(length(checks)==9 && all(isfinite(error) && error<=tolerance for (_,error,tolerance) in checks),"representative all13 derivative policy failed")
    merge(N.flags(),Dict{String,Any}("schema_version"=>1,"status"=>DERIVATIVE_REVIEWED,"scope"=>S.SCOPE,"theory_scope"=>S.declaration(data),
        "source_context"=>deepcopy(binding),"definition"=>deepcopy(spec),"stored_derivative_comparisons_recomputed"=>9,
        "independent_cold_AD_reexecuted_by_reader"=>false,"native_evolution_reexecuted_by_independent_reader"=>false,
        "global_common_order_or_native_selection_performed"=>false,"job_resource_and_complete_provenance_ledger_verified"=>false,
        "all72_production_derivative_certificate"=>false))
end
end
