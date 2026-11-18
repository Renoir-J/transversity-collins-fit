"Complete72 selection arithmetic, with a separate saved-operand reducer."
module SIAGeometryGlobalAudit
using CollinsFit, JSON
include("SIAGeometryReferenceJobs.jl")
const R=SIAGeometryReferenceJobs
const N=R.N
const S=R.S
require(value,message)=value || error(message)

# This reducer deliberately does not call the producer's packed/scaled or
# select_momentum/select_native/select_tails functions. Leaf validation is shared.
function operands(point)
    result=Float64[]
    for pair in ("pp","pm","mp","mm"), field in ("UU","Collins")
        push!(result,point["pairs"][pair][field])
    end
    for category in ("U","L","C"), field in ("UU","Collins")
        push!(result,point["categories"][category][field])
    end
    push!(result,point["amplitude"])
    require(length(result)==15 && all(isfinite,result),"all15 finite absolute/ratio operands required")
    result
end

function difference(current,previous,floor)
    require(length(current)==length(previous)>0,"complete matching point coverage required")
    maximum(maximum(abs.(operands(a).-operands(b))./max.(abs.(operands(b)),floor))
        for (a,b) in zip(current,previous))
end

function momentum(data,records,level,B,parent)
    groups,_=S.geometry_groups(data);policy=parent["numerical_policy"]
    require(length(records)==36,"independent selection requires all36 geometries")
    orders=policy["momentum_orders"];rows=[Vector{Any}(undef,72) for _ in orders]
    for (index,record) in enumerate(records)
        N.audit_value(data,parent,record,index,level,B)
        for (oi,order) in enumerate(orders)
            entry=record["momentum_levels"][oi]
            require(entry["order"]===order,"source order changed")
            for (point,source_index) in zip(entry["points"],groups[index])
                rows[oi][source_index]=point
            end
        end
    end
    selected=nothing;errors=Any[]
    for oi in eachindex(orders)
        residual=oi==1 ? nothing : difference(rows[oi],rows[oi-1],policy["relative_floor"])
        push!(errors,residual)
        if selected===nothing && oi>1 && orders[oi]>=policy["minimum_momentum_order"] && residual<=policy["momentum_max_scaled"]
            selected=oi
        end
    end
    (theory=S.declaration(data),level=level,bound=B,order=selected===nothing ? nothing : orders[selected],errors=errors,allpoints=rows,
        points=selected===nothing ? Any[] : rows[selected])
end

function native(data,records,parent)
    require(length(records) in (3,4),"all native1/2/3 records, optionally all native4, required")
    policy=parent["numerical_policy"];levels=Any[];errors=Any[];selected=nothing
    for (level,values) in enumerate(records)
        row=momentum(data,values,level,24,parent)
        require(row.order!==nothing,"global all72 momentum orders exhausted")
        residual=level==1 ? nothing : difference(row.points,last(levels).points,policy["relative_floor"])
        push!(levels,row);push!(errors,residual)
        if level>=policy["minimum_refinement_level"] && residual<=policy["physical_value_max_scaled"]
            selected=level
            require(level==length(records),"native4 cannot be run after native3 already passes")
        end
    end
    (theory=S.declaration(data),selected=selected,levels=levels,errors=errors)
end

function compare_momentum(actual,expected)
    require(actual["status"]=="MOMENTUM_SELECTION_OBSERVATION_NOT_GATE" && actual["scope"]==S.SCOPE &&
        actual["theory_scope"]==expected.theory &&
        actual["native_level"]===expected.level && actual["bound_inverse_GeV"]===expected.bound &&
        all(get(actual,key,nothing)===false for key in (S.FLAGS...,"accepted_fitting_gate")),"changed momentum scope or setting")
    require(actual["selected_order"]===expected.order && actual["selection_complete"]===(expected.order!==nothing),
        "changed first common all72 momentum decision")
    require(length(actual["momentum_levels"])==5,"all5 momentum observations required")
    for (index,row) in enumerate(actual["momentum_levels"])
        require(row["order"]===expected.allpoints[index][1]["momentum_order"],"changed declared momentum observation order")
        require(S.same_json(row["points"],expected.allpoints[index]) &&
            row["previous_order_max_scaled"]===expected.errors[index],"changed complete momentum operands or residual")
        met=expected.errors[index]!==nothing && expected.errors[index]<=1e-7
        require(row["target_met"]===met && row["accepted"]===(row["order"]>=64 && met),"changed momentum acceptance arithmetic")
    end
    require(S.same_json(actual["points"],expected.points),"selected momentum operands changed")
end

function compare_native(actual,expected)
    require(actual["scope"]==S.SCOPE && actual["theory_scope"]==expected.theory &&
        all(get(actual,key,nothing)===false for key in (S.FLAGS...,"accepted_fitting_gate")),
        "changed native scope or premature admission")
    require(length(actual["levels"])==length(expected.levels) && actual["selected_level"]===expected.selected &&
        actual["selection_complete"]===(expected.selected!==nothing),"changed common native decision")
    for (index,(row,level)) in enumerate(zip(actual["levels"],expected.levels))
        compare_momentum(row,level)
        residual=expected.errors[index]
        require(row["previous_native_level_max_scaled"]===residual &&
            row["native_target_met"]===(residual!==nothing && residual<=1e-4),"changed native refinement residual")
    end
    expected_status=expected.selected!==nothing ? "COMMON_NATIVE_SELECTION_OBSERVATION_NOT_GATE" :
        length(expected.levels)==3 ? "NATIVE4_REQUIRED_NOT_DIAGNOSTIC_GATE" : "NATIVE_REFINEMENT_EXHAUSTED_NOT_GATE"
    require(actual["status"]==expected_status,"changed native exhaustion/continuation status")
end

function native_records(reports)
    require(length(reports) in (1,2) && length(reports[1])==36 &&
        (length(reports)==1 || length(reports[2])==36),"complete native producer reports required")
    values=[[report["native_values"][level] for report in reports[1]] for level in 1:3]
    length(reports)==2 && push!(values,[only(report["native_values"]) for report in reports[2]])
    values
end

function produce_native(data,parent,reports)
    S.select_native(data,native_records(reports),parent)
end

function audit_native(data,parent,reports,result)
    expected=native(data,native_records(reports),parent)
    compare_native(result,expected)
    expected
end

function produce_tails(data,parent,reports,tails)
    S.select_tails(data,native_records(reports),[r["tail_values"][1] for r in tails],
        [r["tail_values"][2] for r in tails],parent)
end

function audit_tails(data,parent,reports,tails,result)
    chosen=native(data,native_records(reports),parent)
    require(chosen.selected!==nothing && length(tails)==36,"selected common native level and all36 tail jobs required")
    lo=momentum(data,[r["tail_values"][1] for r in tails],chosen.selected,12,parent)
    hi=momentum(data,[r["tail_values"][2] for r in tails],chosen.selected,32,parent)
    require(lo.order!==nothing && hi.order!==nothing,"finite-bound common momentum exhausted")
    require(result["selected_level"]===chosen.selected && length(result["bounds"])==3 &&
        result["scope"]==S.SCOPE && result["theory_scope"]==S.declaration(data) &&
        all(get(result,key,nothing)===false for key in (S.FLAGS...,"accepted_fitting_gate")) &&
        result["status"]=="FINITE_TAIL_SELECTION_OBSERVATION_NOT_GATE" && result["infinite_tail_proved"]===false,
        "changed finite-tail scope")
    compare_momentum(result["bounds"][1],lo)
    compare_momentum(result["bounds"][2],last(chosen.levels))
    require(result["bounds"][2]["previous_native_level_max_scaled"]===last(chosen.errors) &&
        result["bounds"][2]["native_target_met"]===(last(chosen.errors)<=1e-4),"B24 native observation differs from selected source operands")
    compare_momentum(result["bounds"][3],hi)
    residual=difference(hi.points,last(chosen.levels).points,parent["numerical_policy"]["relative_floor"])
    require(result["finite24to32_max_scaled"]===residual && result["finite_tail_target_met"]===(residual<=1e-4),
        "changed all72 finite-tail comparison")
    (native=chosen,lo=lo,hi=hi,error=residual)
end

function reference_errors(data,parent,tails,selected)
    groups,_=S.geometry_groups(data);errors=Vector{Float64}(undef,72)
    require(length(tails)==36 && length(selected.points)==72,"all72 scalar reference comparisons required")
    for (report,indices) in zip(tails,groups)
        require(length(report["independent_references"])==2,"both source scalar references required")
        for (record,index) in zip(report["independent_references"],indices)
            reference=R.validate_reference(record["reference"],data["points"][index],parent["numerical_policy"])
            errors[index]=difference([selected.points[index]],[reference],parent["numerical_policy"]["relative_floor"])
            comparison=only(filter(x->x["order"]==selected.order,record["order_comparisons"]))
            require(comparison["max_scaled"]===errors[index],"changed selected-order scalar reference residual")
        end
    end
    errors
end
end
