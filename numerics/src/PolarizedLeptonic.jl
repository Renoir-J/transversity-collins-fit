"Explicit impulse-model weights; nuclear prescriptions require a supplied source."
struct NucleonTarget
    name::String
    unpolarized::NTuple{2,Float64}
    polarized::NTuple{2,Float64}
    source::String
    function NucleonTarget(name::String,u::NTuple{2,Float64},p::NTuple{2,Float64},source::String)
        !isempty(name) && !isempty(source) && all(isfinite,u) && all(isfinite,p) &&
            all(>=(0),u) && sum(u)>0 && all(abs(p[i])<=u[i] for i in 1:2) || error("invalid explicit nucleon target")
        new(name,u,p,source)
    end
end

function NucleonTarget(name::AbstractString,unpolarized,polarized;source)
    !isempty(name) && !isempty(source) && length(unpolarized)==length(polarized)==2 ||
        error("named target, source and two explicit nucleon weights required")
    u=Tuple(fixed_scalar(v,"unpolarized target weight") for v in unpolarized)
    p=Tuple(fixed_scalar(v,"polarized target weight") for v in polarized)
    all(>=(0),u) && sum(u)>0 && all(abs(p[i])<=u[i] for i in 1:2) || error("invalid nucleon target weights")
    NucleonTarget(String(name),u,p,String(source))
end
function NucleonTarget(name::Symbol)
    name in (:proton,:neutron) || error("nuclear targets require explicit weights and a source")
    weights=name==:proton ? (1.,0.) : (0.,1.)
    NucleonTarget(string(name),weights,weights;source="isolated nucleon with exact u/d isospin interchange")
end

"Fixed source samples and logarithms; parameter-dependent factors remain generic."
struct PolarizedInitialSamples
    kind::Symbol
    grid::Vector{Float64}
    logarithm::Vector{Float64}
    complement_logarithm::Vector{Float64}
    bases::Matrix{Float64}
    input_kind::Symbol
    signature::UInt
end

function PolarizedInitialSamples(kind,grid,bases;input_kind)
    kind in (:transversity,:KPSY_Hhat3_homogeneous) || error("unsupported initial sample kind")
    input_kind in (:pinned_inputs,:analytic) || error("explicit initial source kind required")
    x=collinear_grid(grid);b=Matrix{Float64}(bases)
    expected_rows=kind==:transversity ? 2 : 4
    size(b)==(expected_rows,length(x)) && all(isfinite,b) && all(iszero,b[:,end]) ||
        error("invalid polarized initial bases or endpoint")
    kind==:transversity && !all(>=(0),b) && error("negative Soffer source; no clipping is allowed")
    lx=log.(x);lc=[xx==1 ? 0. : log1p(-xx) for xx in x]
    PolarizedInitialSamples(kind,x,lx,lc,b,input_kind,hash((kind,x,lx,lc,b,input_kind)))
end

function polarized_initial_columns(samples::PolarizedInitialSamples,p::KPSYParameters)
    hash((samples.kind,samples.grid,samples.logarithm,samples.complement_logarithm,samples.bases,samples.input_kind))==samples.signature ||
        error("mutated polarized initial samples")
    result=zeros(typeof(p.values[1]),length(samples.grid),size(samples.bases,1))
    if samples.kind==:transversity
        for (column,(ni,ai,bi)) in enumerate(((1,3,5),(2,4,6)))
            n,a,b=p.values[ni],p.values[ai],p.values[bi]
            # Reuse parameter-only entropy terms, preserving the original order.
            ab,aa,bb=entropy_term(a+b),entropy_term(a),entropy_term(b)
            for j in 1:length(samples.grid)-1
                result[j,column]=n*exp(a*samples.logarithm[j]+b*samples.complement_logarithm[j]+ab-aa-bb)*samples.bases[column,j]
            end
        end
    else
        for j in 1:length(samples.grid)-1
            favored=p.values[7]*exp(p.values[9]*samples.logarithm[j]+p.values[11]*samples.complement_logarithm[j])
            unfavored=p.values[8]*exp(p.values[10]*samples.logarithm[j]+p.values[12]*samples.complement_logarithm[j])
            result[j,1]=favored*samples.bases[1,j]
            for column in 2:4
                result[j,column]=unfavored*samples.bases[column,j]
            end
        end
    end
    result
end

struct PolarizedLeptonicPlan
    leptonic::LeptonicPlan
    unpolarized_first::Matrix{Float64}
    unpolarized_second::Matrix{Float64}
    unpolarized::Float64
    first_operator::HoppetMatchingPlan
    second_operator::Union{Nothing,HoppetMatchingPlan}
    first_samples::PolarizedInitialSamples
    second_samples::Union{Nothing,PolarizedInitialSamples}
    target::Union{Nothing,NucleonTarget}
    charges::NTuple{2,Int}
    provenance::Dict{String,Any}
    signature::UInt
end

pion_charge(charge)=charge isa Integer && !(charge isa Bool) && charge in (-1,0,1) ? Int(charge) : error("fixed pion charge required")
isospin_flavor(pid)=abs(pid)==1 ? 2sign(pid) : abs(pid)==2 ? sign(pid) : pid
isospin_indices(flavors)=[something(findfirst(==(isospin_flavor(pid)),flavors)) for pid in flavors]

function initial_model_settings(model)
    model.Q02==2.4 && model.coupling.nf==5 && model.coupling.lambda_GeV==.225 ||
        error("initial family requires Q02=2.4 and the fixed nf5/Lambda=.225 coupling")
end

function polarized_operator_valid(op,samples,model,scales,process,quantity,queries)
    op.process==process && op.quantity==quantity && op.nf==5 && op.Q02==model.Q02 &&
        op.Q2==scales && op.queries==queries && op.grid==samples.grid && samples.kind==quantity &&
        op.alpha==alpha_s.(Ref(model.coupling),scales) && op.tau==evolution_time.(Ref(model.coupling),scales,model.Q02) ||
        error("polarized operator/source/process/scale mismatch")
end

function polarized_context_signature(plan,first,second,target,charges,provenance)
    model=plan.model
    hash((plan.process,plan.Q2,plan.P,plan.z1,plan.z2,plan.flavors,plan.antiparticle,plan.charge2,
        plan.radial.b,plan.radial.db,plan.radial.edges,plan.radial.order,
        plan.uu_weights,plan.collins_weights,plan.uu_damping,plan.fixed_exponent,plan.fragmentation_shape,
        model.bmax,model.Q02,model.g2,model.gq,model.gh,model.coupling.lambda_GeV,model.coupling.nf,model.nodes,model.weights,
        first,second,target,charges,provenance))
end

function pion_charge_matrix(basis,charge,antiparticle)
    charge==1 ? basis : charge==-1 ? basis[antiparticle,:] : (basis+basis[antiparticle,:])/2
end

"Assemble pi+ unpolarized bases and prepared polarized operators."
function polarized_leptonic_plan(plan::LeptonicPlan,first,second,first_operator,first_samples;
                                 second_operator=nothing,second_samples=nothing,target=nothing,
                                 charges=(1,0),provenance=Dict{String,Any}())
    initial_model_settings(plan.model)
    plan.flavors==INITIAL_FLAVORS || error("complete ordered historical flavor basis required")
    u1,u2=Matrix{Float64}(first),Matrix{Float64}(second)
    matched_matrix(plan,u1);matched_matrix(plan,u2)
    length(charges)==2 || error("two charge slots required")
    cs=(pion_charge(charges[1]),pion_charge(charges[2]))
    scales=css_mu_b2.(Ref(plan.model),plan.radial.b)
    if plan.process==:SIDIS
        cs[2]==0 || error("SIDIS has one pion charge and a zero second charge slot")
        target isa NucleonTarget && second_operator isa HoppetMatchingPlan && second_samples isa PolarizedInitialSamples ||
            error("SIDIS requires explicit target and two polarized legs")
        length(first_operator.queries)==1 || error("SIDIS has a single x query")
        polarized_operator_valid(first_operator,first_samples,plan.model,scales,:SIDIS,:transversity,first_operator.queries)
        polarized_operator_valid(second_operator,second_samples,plan.model,scales,:SIDIS,:KPSY_Hhat3_homogeneous,[plan.z1])
        tfirst=target.unpolarized[1].*u1+target.unpolarized[2].*u1[isospin_indices(plan.flavors),:]
        tsecond=pion_charge_matrix(u2,cs[1],plan.antiparticle)
    else
        isnothing(target) && isnothing(second_operator) && isnothing(second_samples) || error("SIA uses two queries in one fragmentation operator")
        polarized_operator_valid(first_operator,first_samples,plan.model,scales,:SIA,:KPSY_Hhat3_homogeneous,[plan.z1,plan.z2])
        tfirst=pion_charge_matrix(u1,cs[1],plan.antiparticle)
        tsecond=pion_charge_matrix(u2,cs[2],plan.antiparticle)
    end
    get(provenance,"scope","reconstructed_theory")=="reconstructed_theory" && get(provenance,"production_ready",false)===false &&
        get(provenance,"physical_extraction",false)===false ||
        error("polarized component cannot consume a production admission claim")
    meta=deepcopy(provenance)
    meta["scope"]="reconstructed_theory";meta["production_ready"]=false;meta["physical_extraction"]=false
    meta["contract_id"]="KPSY-POLARIZED-PROCESS-001"
    meta["process"]=string(plan.process);meta["charges"]=collect(cs)
    meta["unpolarized_fragmentation_basis"]="pi+; charge conjugation or neutral average applied in assembly"
    meta["initial_source_kind"]=string(first_samples.input_kind)
    !isnothing(second_samples) && (meta["second_initial_source_kind"]=string(second_samples.input_kind))
    !isnothing(target) && (meta["target"]=Dict("name"=>target.name,"unpolarized"=>collect(target.unpolarized),"polarized"=>collect(target.polarized),"source"=>target.source))
    denominator=leptonic_unpolarized(plan,tfirst,tsecond)
    signature=polarized_context_signature(plan,u1,u2,target,cs,meta)
    PolarizedLeptonicPlan(plan,u1,u2,denominator,first_operator,second_operator,first_samples,second_samples,target,cs,meta,signature)
end

function with_pion_charges(context::PolarizedLeptonicPlan,charge1,charge2=context.charges[2])
    polarized_context_valid(context)
    context.leptonic.process==:SIDIS && charge2!=0 && error("SIDIS has one pion charge")
    polarized_leptonic_plan(context.leptonic,context.unpolarized_first,context.unpolarized_second,context.first_operator,context.first_samples;
        second_operator=context.second_operator,second_samples=context.second_samples,target=context.target,
        charges=(charge1,charge2),provenance=context.provenance)
end

function polarized_context_valid(context)
    polarized_context_signature(context.leptonic,context.unpolarized_first,context.unpolarized_second,
        context.target,context.charges,context.provenance)==context.signature || error("mutated polarized process context")
end

function collins_flavor_columns(bases,charge,flavors)
    result=zeros(eltype(bases),length(flavors),size(bases,2))
    for (f,pid) in enumerate(flavors)
        if abs(pid)==3
            result[f,:]=bases[pid==3 ? 3 : 4,:]
        elseif abs(pid)<3
            result[f,:]=charge==0 ? (bases[1,:]+bases[2,:])/2 :
                bases[(charge==1 ? pid in (2,-1) : pid in (-2,1)) ? 1 : 2,:]
        end
    end
    result
end

"Transient responses for one parameter value; never cached across parameter calls."
function evolved_polarized_columns(context::PolarizedLeptonicPlan,p::KPSYParameters)
    polarized_context_valid(context)
    p.gh==context.leptonic.model.gh || error("initial/CSS width mismatch")
    first=apply_operator(context.first_operator,polarized_initial_columns(context.first_samples,p))
    second=context.leptonic.process==:SIDIS ?
        apply_operator(context.second_operator,polarized_initial_columns(context.second_samples,p)) : nothing
    first,second
end

function polarized_flavor_legs(context,first,second;charges=context.charges)
    length(charges)==2 || error("two charge slots required")
    cs=(pion_charge(charges[1]),pion_charge(charges[2]))
    flavors=context.leptonic.flavors
    if context.leptonic.process==:SIDIS
        cs[2]==0 || error("SIDIS has one pion charge")
        hp=zeros(eltype(first),length(flavors),size(first,3))
        for (f,pid) in enumerate(flavors)
            pid in (2,1) && (hp[f,:]=first[1,pid==2 ? 1 : 2,:])
        end
        target=context.target
        pol1=target.polarized[1].*hp+target.polarized[2].*hp[isospin_indices(flavors),:]
        pol2=collins_flavor_columns(second[1,:,:],cs[1],flavors)
    else
        pol1=collins_flavor_columns(first[1,:,:],cs[1],flavors)
        pol2=collins_flavor_columns(first[2,:,:],cs[2],flavors)
    end
    pol1,pol2
end

"Complete reconstructed structure functions; angular factors and bin adapters are separate."
function polarized_structure(context::PolarizedLeptonicPlan,p::KPSYParameters)
    first,second=evolved_polarized_columns(context,p)
    pol1,pol2=polarized_flavor_legs(context,first,second)
    (UU=context.unpolarized,Collins=leptonic_collins(context.leptonic,pol1,pol2,p.values[13]))
end
polarized_structure(context::PolarizedLeptonicPlan,theta::AbstractVector)=polarized_structure(context,KPSYParameters(theta;gh=context.leptonic.model.gh))

"Reuse prepared evolution/input/CSS arrays while changing only fixed Fourier momentum."
function with_momentum(context::PolarizedLeptonicPlan,P)
    polarized_context_valid(context)
    old=context.leptonic;p=fixed_scalar(P,"P");p>=0 || error("nonnegative momentum required")
    k=p/old.z1
    uu=hankel_weights(old.radial,0,k)./old.z1^2
    pol=old.process==:SIDIS ? -hankel_weights(old.radial,1,k)./(2old.z1^3) :
        hankel_weights(old.radial,2,k)./(4old.z1^3*old.z2)
    plan=LeptonicPlan(old.process,old.model,old.radial,old.Q2,p,old.z1,old.z2,copy(old.flavors),copy(old.antiparticle),copy(old.charge2),
        uu,pol,copy(old.uu_damping),copy(old.fixed_exponent),copy(old.fragmentation_shape))
    polarized_leptonic_plan(plan,context.unpolarized_first,context.unpolarized_second,context.first_operator,context.first_samples;
        second_operator=context.second_operator,second_samples=context.second_samples,target=context.target,charges=context.charges,provenance=context.provenance)
end

function prepare_polarized_sidis(pdf::CT10Input,hel::LegacyPartonInput,ff::LegacyPartonInput,buildfile,
                                  model,radial,Q2,P,x,z,charge;matching_order=8,intervals=1024,
                                  range_policy=:strict,target=NucleonTarget(:proton),construction=:batched_probes)
    initial_model_settings(model);pion_charge(charge)
    target isa NucleonTarget || error("explicit nucleon target prescription required")
    unpol=prepare_unpolarized_sidis(pdf,ff,model,radial,Q2,P,x,z,1;matching_order=matching_order,range_policy=range_policy)
    scales=css_mu_b2.(Ref(model),radial.b)
    h=HoppetMatchingPlan(buildfile,[x],scales,model.Q02,model.coupling;process=:SIDIS,quantity=:transversity,intervals=intervals,construction=construction)
    H=HoppetMatchingPlan(buildfile,[z],scales,model.Q02,model.coupling;process=:SIDIS,intervals=intervals,construction=construction)
    hx,hb=prepare_initial_soffer(pdf,hel,h.grid);Hz,Hb=prepare_initial_ff(ff,H.grid)
    hs=PolarizedInitialSamples(:transversity,hx,hb;input_kind=:pinned_inputs)
    Hs=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,Hz,Hb;input_kind=:pinned_inputs)
    polarized_leptonic_plan(unpol.plan,unpol.first,unpol.second,h,hs;
        second_operator=H,second_samples=Hs,target=target,charges=(charge,0),provenance=unpol.provenance)
end

function prepare_polarized_sia(ff::LegacyPartonInput,buildfile,model,radial,Q2,P,z1,z2,charge1,charge2;
                              matching_order=8,intervals=1024,range_policy=:strict,frame=:hadron2,
                              construction=:batched_probes)
    initial_model_settings(model);pion_charge(charge1);pion_charge(charge2)
    frame==:hadron2 || error("SIA structure function requires the hadron-2 frame")
    unpol=prepare_unpolarized_sia(ff,model,radial,Q2,P,z1,z2,1,1;matching_order=matching_order,range_policy=range_policy)
    scales=css_mu_b2.(Ref(model),radial.b)
    H=HoppetMatchingPlan(buildfile,[z1,z2],scales,model.Q02,model.coupling;process=:SIA,intervals=intervals,construction=construction)
    z,bases=prepare_initial_ff(ff,H.grid)
    samples=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,z,bases;input_kind=:pinned_inputs)
    polarized_leptonic_plan(unpol.plan,unpol.first,unpol.second,H,samples;charges=(charge1,charge2),provenance=unpol.provenance)
end
