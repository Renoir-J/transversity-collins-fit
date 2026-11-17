"A typed geometry view keeps the complete source request intact."
module SIAGeometryModel
using CollinsFit, JSON
include("SIAGeometrySelection.jl")
const S=SIAGeometrySelection
const VIEW_SCOPE="source_bound_two_ratio_geometry_view_not_full72_request"
require(value,message)=value || error(message)

function selected_data(full,index,indices)
    Dict{String,Any}("scope"=>VIEW_SCOPE,"theory_scope"=>S.declaration(full.data),
        "full_request_sha256"=>full.sha256,"geometry_index"=>index,"source_point_indices"=>collect(indices),
        "points"=>deepcopy(full.data["points"][collect(indices)]),
        "momentum_window"=>deepcopy(full.data["momentum_window"]))
end

"An owned two-row view, explicitly distinct from HistoricalSIAIntegralRequest."
struct SelectedSIAGeometry
    full::HistoricalSIAIntegralRequest
    index::Int
    indices::NTuple{2,Int}
    data::Dict{String,Any}
    sha256::String
    signature::String
    function SelectedSIAGeometry(full::HistoricalSIAIntegralRequest,index)
        historical_sia_integral_request_valid(full)
        groups,_=S.geometry_groups(full.data)
        require(index isa Integer && !(index isa Bool) && 1<=index<=36,"select one declared geometry index")
        i=Int(index);indices=Tuple(groups[i]);data=selected_data(full,i,indices)
        signature=CollinsFit.numerical_snapshot((full.signature,i,indices,data,full.sha256))
        new(full,i,indices,data,full.sha256,signature)
    end
end

function valid(view::SelectedSIAGeometry)
    historical_sia_integral_request_valid(view.full)
    groups,_=S.geometry_groups(view.full.data)
    require(view.indices==Tuple(groups[view.index]) && view.sha256==view.full.sha256 &&
        S.same_json(view.data,selected_data(view.full,view.index,view.indices)) &&
        CollinsFit.numerical_snapshot((view.full.signature,view.index,view.indices,view.data,view.sha256))==view.signature,
        "mutated or substituted selected geometry view")
    view
end

"Use the inherited observable construction for this owned view, without changing its full-request method."
function CollinsFit.historical_sia_integral_model(view::SelectedSIAGeometry,contexts;order=64)
    valid(view)
    contract=JSON.parsefile(CollinsFit.HISTORICAL_SIA_INTEGRAL_CONTRACT)
    require(order isa Integer && !(order isa Bool) && order in contract["numerical_policy"]["momentum_orders"],
        "undeclared momentum order")
    key=historical_sia_integral_geometry_key(first(view.data["points"]))
    require(Set(keys(contexts))==Set((key,)),"one exact selected geometry context required")
    w=view.data["momentum_window"];lo,hi=w["qT_bounds_GeV"]
    window=SIAQTWindow(lo,hi;source=w["source"],lower_limit_source=w["lower_limit_source"])
    observations=JointObservation[]
    for point in view.data["points"]
        context=contexts[key];plan=context.leptonic;k=point["kinematics"]
        require(historical_sia_integral_geometry_key(point)==key && plan.process==:SIA &&
            plan.Q2==k["Q2_GeV2"] && plan.P==contract["prediction"]["preparation_reference_P_GeV"] &&
            plan.z1==k["z1"] && plan.z2==k["z2"] && collect(context.first_operator.queries)==[k["z1"],k["z2"]] &&
            context.second_operator===nothing,"SIA context geometry or operator mismatch")
        a=point["source_angular"]
        angular=SIAMeanOfRatio(a["value"];source="BaBar RF0 primary table; request SHA256="*view.sha256,
            averaging_measure=a["averaging_measure"])
        observable=sia_qt_observable(context,window;ratio=Symbol(point["ratio"]),angular=angular,order=order)
        push!(observations,JointObservation(point["id"],observable;source=point["publication"]*" RF0 source row"))
    end
    JointPredictionModel(observations;source="source-bound conditional historical SIA integral diagnostic; request SHA256="*view.sha256)
end
end
