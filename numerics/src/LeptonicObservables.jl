"Fixed finite momentum quadrature including the physical P*dP measure."
struct MomentumRule
    nodes::Vector{Float64}
    weights::Vector{Float64}
    lower::Float64
    upper::Float64
    order::Int
    signature::UInt
end

function MomentumRule(lower,upper;order=32)
    low,high=fixed_scalar(lower,"Pmin"),fixed_scalar(upper,"Pmax")
    0<=low<high || error("finite momentum interval 0<=Pmin<Pmax required")
    order isa Integer && !(order isa Bool) || error("fixed integer momentum order required")
    nodes,weights=gauss_rule(order)
    half=(high-low)/2
    P=(low+high)/2 .+ half.*nodes
    w=half.*weights.*P
    MomentumRule(P,w,low,high,Int(order),hash((P,w,low,high,Int(order))))
end

momentum_signature(rule::MomentumRule)=hash((rule.nodes,rule.weights,rule.lower,rule.upper,rule.order))
momentum_valid(rule::MomentumRule)=momentum_signature(rule)==rule.signature || error("mutated momentum quadrature")

"Source-labeled ratio of angular averages for hadron-2-frame SIA."
struct SIAAngularFactor
    value::Float64
    source::String
    function SIAAngularFactor(value::Float64,source::String)
        isfinite(value) && 0<=value<=1 && !isempty(source) || error("finite ratio of angular means in [0,1] and source required")
        new(value,source)
    end
end
SIAAngularFactor(value;source)=SIAAngularFactor(fixed_scalar(value,"ratio of angular means"),String(source))

"Source-labeled mean of pointwise angular ratios, with its declared averaging measure."
struct SIAMeanOfRatio
    value::Float64
    source::String
    averaging_measure::String
    function SIAMeanOfRatio(value::Float64,source::String,averaging_measure::String)
        isfinite(value) && 0<=value<=1 || error("finite mean of angular ratios in [0,1] required")
        !isempty(strip(source)) && !isempty(strip(averaging_measure)) || error("explicit angular source and averaging measure required")
        new(value,source,averaging_measure)
    end
end
SIAMeanOfRatio(value;source,averaging_measure)=SIAMeanOfRatio(
    fixed_scalar(value,"mean of angular ratios"),String(source),String(averaging_measure))

sia_angular_metadata(a::SIAAngularFactor)=Dict{String,Any}(
    "angular_definition"=>"ratio_of_angular_means","angular_source"=>a.source,"angular_factor"=>a.value)
sia_angular_metadata(a::SIAMeanOfRatio)=Dict{String,Any}(
    "angular_definition"=>"mean_of_angular_ratio","angular_source"=>a.source,"angular_factor"=>a.value,
    "angular_averaging_measure"=>a.averaging_measure)

function depolarization_factor(y)
    value=fixed_scalar(y,"fixed SIDIS y")
    0<=value<=1 || error("SIDIS y must lie in [0,1]")
    2(1-value)/(1+(1-value)^2)
end

function physical_denominator(value)
    value isa Real && !(value isa Bool) && isfinite(value) && value>0 || error("finite positive final unpolarized denominator required")
    value
end

"Sum structure functions before forming any unlike/like/charged ratio."
function sia_charge_sums(pp,pm,mp,mm)
    all(r->all(v->v isa Real && !(v isa Bool) && isfinite(v),(r.UU,r.Collins)),(pp,pm,mp,mm)) || error("finite real charge-pair structure functions required")
    U=(UU=pm.UU+mp.UU,Collins=pm.Collins+mp.Collins)
    L=(UU=pp.UU+mm.UU,Collins=pp.Collins+mm.Collins)
    (U=U,L=L,C=(UU=U.UU+L.UU,Collins=U.Collins+L.Collins))
end

function sia_double_ratio(sums,ratio,angular::Union{SIAAngularFactor,SIAMeanOfRatio})
    ratio in (:UL,:UC) || error("SIA double ratio must be UL or UC")
    sums.C.UU==sums.U.UU+sums.L.UU && sums.C.Collins==sums.U.Collins+sums.L.Collins || error("charged category must be the sum of unlike and like categories")
    denominator=ratio==:UL ? sums.L : sums.C
    angular.value*(sums.U.Collins/physical_denominator(sums.U.UU)-denominator.Collins/physical_denominator(denominator.UU))
end

struct LeptonicObservablePlan
    kind::Symbol
    context::PolarizedLeptonicPlan
    momentum::Union{Nothing,MomentumRule}
    uu_weights::Vector{Float64}
    collins_weights::Vector{Float64}
    denominators::Vector{Float64}
    angular::Float64
    provenance::Dict{String,Any}
    signature::UInt
end

function observable_signature(kind,context,momentum,uu,collins,denominators,angular,metadata)
    hash((kind,context.signature,momentum===nothing ? nothing : momentum_signature(momentum),uu,collins,denominators,angular,metadata))
end

function observable_weights(context,momentum)
    plan=context.leptonic
    momentum===nothing && return copy(plan.uu_weights),copy(plan.collins_weights)
    momentum isa MomentumRule || error("explicit momentum quadrature required")
    momentum_valid(momentum)
    uu=zeros(length(plan.radial.b));collins=copy(uu)
    for (P,w) in zip(momentum.nodes,momentum.weights)
        k=P/plan.z1
        uu .+= w.*hankel_weights(plan.radial,0,k)./plan.z1^2
        collins .+= plan.process==:SIDIS ? -w.*hankel_weights(plan.radial,1,k)./(2plan.z1^3) :
            w.*hankel_weights(plan.radial,2,k)./(4plan.z1^3*plan.z2)
    end
    uu,collins
end

function observable_unpolarized(context,weights,c1,c2)
    plan=context.leptonic
    first=context.unpolarized_first;second=context.unpolarized_second
    if plan.process==:SIDIS
        t=context.target
        A=t.unpolarized[1].*first+t.unpolarized[2].*first[isospin_indices(plan.flavors),:]
        D=pion_charge_matrix(second,c1,plan.antiparticle)
    else
        A=pion_charge_matrix(first,c1,plan.antiparticle)
        D=pion_charge_matrix(second,c2,plan.antiparticle)
    end
    finite_transform(weights,flavor_product(plan,A,D).*plan.uu_damping)
end

function assemble_observable(kind,context,momentum,angular,source;extra=Dict{String,Any}())
    polarized_context_valid(context)
    allowed=context.leptonic.process==:SIDIS ? (:structure_function_ratio,:sin_phi_sum_amplitude) : (:UL,:UC)
    kind in allowed || error("observable definition/process mismatch")
    source isa AbstractString && !isempty(source) || error("explicit observable source required")
    uu,collins=observable_weights(context,momentum)
    denominators=if context.leptonic.process==:SIDIS
        [observable_unpolarized(context,uu,context.charges...)]
    else
        U=observable_unpolarized(context,uu,1,-1)+observable_unpolarized(context,uu,-1,1)
        L=observable_unpolarized(context,uu,1,1)+observable_unpolarized(context,uu,-1,-1)
        [U,L]
    end
    foreach(physical_denominator,denominators)
    metadata=Dict{String,Any}("contract_id"=>"KPSY-LEPTONIC-OBS-002","scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"process"=>string(context.leptonic.process),
        "definition"=>string(kind),"source"=>String(source),"context_provenance"=>deepcopy(context.provenance),
        "kinematics"=>(momentum===nothing ? "fixed_kinematics" : "fixed_other_kinematics_finite_P_integral"),
        "acceptance"=>"uniform; explicit reconstruction assumption, not an experimental correction")
    if momentum!==nothing
        metadata["momentum_measure"]="P*dP"
        metadata["momentum_bounds_GeV"]=[momentum.lower,momentum.upper]
        metadata["momentum_order"]=momentum.order
    else
        metadata["P_GeV"]=context.leptonic.P
    end
    merge!(metadata,deepcopy(extra))
    metadata["scope"]=="reconstructed_theory" && metadata["production_ready"]===false && metadata["physical_extraction"]===false || error("observable cannot claim production admission")
    signature=observable_signature(kind,context,momentum,uu,collins,denominators,angular,metadata)
    LeptonicObservablePlan(kind,context,momentum,uu,collins,denominators,angular,metadata,signature)
end

function sidis_observable(context::PolarizedLeptonicPlan;definition,y=nothing,momentum=nothing,source)
    context.leptonic.process==:SIDIS || error("SIDIS context required")
    definition in (:structure_function_ratio,:sin_phi_sum_amplitude) || error("explicit SIDIS depolarization convention required")
    value=y===nothing ? nothing : fixed_scalar(y,"fixed SIDIS y")
    value!==nothing && !(0<value<=1) && error("physical SIDIS observable requires 0<y<=1")
    definition==:sin_phi_sum_amplitude && value===nothing && error("raw SIDIS amplitude requires y")
    angular=definition==:sin_phi_sum_amplitude ? depolarization_factor(value) : 1.
    assemble_observable(definition,context,momentum,angular,source;
        extra=Dict{String,Any}("y"=>value,"angle_convention"=>"sin(phi_h+phi_s)",
            "depolarization_removed"=>definition==:structure_function_ratio))
end

function sia_observable(context::PolarizedLeptonicPlan;ratio,angular::Union{SIAAngularFactor,SIAMeanOfRatio},momentum=nothing,source)
    context.leptonic.process==:SIA || error("SIA hadron2 context required")
    ratio in (:UL,:UC) || error("SIA double ratio must be UL or UC")
    extra=Dict{String,Any}("frame"=>"hadron2","double_ratio_approximation"=>"linearized")
    merge!(extra,sia_angular_metadata(angular))
    assemble_observable(ratio,context,momentum,angular.value,source;extra=extra)
end

function observable_valid(observable)
    observable.momentum!==nothing && momentum_valid(observable.momentum)
    polarized_context_valid(observable.context)
    observable_signature(observable.kind,observable.context,observable.momentum,observable.uu_weights,
        observable.collins_weights,observable.denominators,observable.angular,observable.provenance)==observable.signature || error("mutated observable settings, weights or metadata")
end

function predict_observable(observable::LeptonicObservablePlan,p::KPSYParameters)
    observable_valid(observable)
    context=observable.context;plan=context.leptonic
    first,second=evolved_polarized_columns(context,p)
    damping=exp.(-plan.fixed_exponent.-collins_width(plan.model,p.values[13]).*plan.fragmentation_shape)
    function polarized(c1,c2)
        A,D=polarized_flavor_legs(context,first,second;charges=(c1,c2))
        finite_transform(observable.collins_weights,flavor_product(plan,A,D).*damping)
    end
    if plan.process==:SIDIS
        observable.angular*polarized(context.charges...)/observable.denominators[1]
    else
        U=polarized(1,-1)+polarized(-1,1)
        L=polarized(1,1)+polarized(-1,-1)
        right=observable.kind==:UL ? L/observable.denominators[2] : (U+L)/sum(observable.denominators)
        observable.angular*(U/observable.denominators[1]-right)
    end
end
predict_observable(observable::LeptonicObservablePlan,theta::AbstractVector)=predict_observable(observable,KPSYParameters(theta;gh=observable.context.leptonic.model.gh))
