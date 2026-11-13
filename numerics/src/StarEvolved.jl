"STAR C11 gluon coefficients: log(v) removes the splitting-function plus singularity."
function star_gluon_regular(v)
    0<v<=1 || error("matching fraction must be in (0,1]")
    L=v==1 ? -one(v) : (v>.5 ? log1p(v-1) : log(v))/(1-v)
    gg=6*(v*L+(1-v)*log(v)/v+v*(1-v)*log(v))
    qg=(v*(1-v)+(v^2+(1-v)^2)*log(v))/2
    gg,qg
end

function star_gluon_matching_rule(z,alpha;order=8,knots=Float64[])
    zz=css_z(z);a=fixed_scalar(alpha,"matching alpha_s")
    a>=0 || error("matching alpha_s must be nonnegative")
    grid=[fixed_scalar(x,"matching knot") for x in knots]
    all(x->0<x<=1,grid) && all(diff(grid).>0) || error("matching knots must increase in (0,1]")
    nodes,weights=gauss_rule(order)
    coordinates,wg,wq=[zz],[1.],[0.]
    if a>0 && zz<1
        edges=vcat(zz,filter(x->zz<x<1,grid),1.)
        for i in 1:length(edges)-1
            lo,hi=log(edges[i]),log(edges[i+1]);half=(hi-lo)/2
            for j in eachindex(nodes)
                y=exp((lo+hi)/2+half*nodes[j]);gg,qg=star_gluon_regular(zz/y)
                push!(coordinates,y);push!(wg,a/pi*half*weights[j]*gg);push!(wq,a/pi*half*weights[j]*qg)
            end
        end
    end
    (coordinates=coordinates,gluon_weights=wg,singlet_weights=wq)
end

"Single gluon Sudakov, with the paper's signed integration limits and nf5 coefficients."
function star_gluon_sudakov(model::KPSYCSS,b,Q2)
    coupling_log(model.coupling,Q2)
    function coefficients(q2)
        a=alpha_s(model.coupling,q2);A,_=css_coefficients(a,model.coupling.nf)
        ((9/4)*A,-(11-2model.coupling.nf/3)*a/(2pi))
    end
    css_integral(css_mu_b2(model,b),Q2,coefficients,model.nodes,model.weights)/2
end

function star_species(pids,fragmentation,gluon_model)
    all(p->p isa Integer && !(p isa Bool),pids) || error("explicit integer PDG basis required")
    ids=Tuple(Int(p) for p in pids)
    length(ids)==11 && Set(ids)==Set((INITIAL_FLAVORS...,21)) || error("complete massless nf5 PDG basis required")
    fragmentation in (:quark_only,:quark_and_gluon) || error("explicit outgoing fragmentation choice required")
    expected=fragmentation==:quark_only ? :none : :casimir_common_same_intrinsic_1705
    gluon_model==expected || error("explicit gluon prescription does not match fragmentation choice")
    ids
end

star_sample_payload(s)=(s.kind,s.grid,s.logarithm,s.complement_logarithm,s.bases,s.input_kind)
star_model_payload(m)=(m.bmax,m.Q02,m.g2,m.gq,m.gh,m.coupling.lambda_GeV,m.coupling.nf,m.nodes,m.weights)
star_radial_payload(r)=(r.b,r.db,r.edges,r.order)

"Reject native precision errors even when a retained ABI reports status zero."
function star_native_checked(f)
    mktemp() do _,stream
        result=redirect_stderr(stream) do
            f()
        end
        flush(stream);seekstart(stream);diagnostic=read(stream,String)
        occursin("TOO HIGH ACCURACY REQUIRED",diagnostic) &&
            error("unqualified HOPPET preparation: native integrator returned zero on a precision error; captured stderr: $diagnostic")
        result
    end
end

"One stable common native beam grid; detach the full evolution/interpolation map."
function star_beam_operator(buildfile,queries,Q2,model::KPSYCSS;intervals=1024,construction=:batched_probes)
    initial_model_settings(model)
    intervals isa Integer && !(intervals isa Bool) && 11<=intervals<=20000 || error("fixed beam intervals11..20000 required")
    qs=[fixed_scalar(x,"STAR beam query") for x in queries]
    !isempty(qs) && all(x->0<x<=1,qs) && minimum(qs)<1 || error("STAR beam queries must include an interior point")
    q2=positive_fixed(Q2,"STAR beam Q2")
    verification=HoppetBuildVerification(verified_hoppet_build(buildfile))
    contract=JSON.parsefile(normpath(joinpath(@__DIR__,"../../config/hoppet-point-contract.json")))
    star_native_checked() do
        native=HoppetPlan(verification,[q2],model.Q02,model.coupling;xmin=minimum(qs),
            dy=-log(minimum(qs))/intervals,layout=:single,construction=construction)
        try
            n=length(native.grid);basis=Matrix{Float64}(I,n,n);basis[end,:].=0.
            response=hoppet_evaluate(native,apply_operator(native,basis),qs)
            blocks=[x==1 ? HoppetPointBlock(Int[],zeros(1,0)) :
                HoppetPointBlock(collect(1:n),permutedims(response[i,:,:])) for (i,x) in enumerate(qs)]
            probe=hcat(native.grid.*(1 .-native.grid).^3,native.grid.^2 .* (1 .-native.grid).^2)
            reference=hoppet_evaluate(native,apply_operator(native,probe),qs)
            actual=cat([x==1 ? zeros(1,2) : b.rows*probe for (x,b) in zip(qs,blocks)]...;dims=1)
            expected=reference[:,:,1]
            row_error=maximum(abs.(actual-expected)./max.(abs.(expected),.01maximum(abs,expected),1e-15))
            isfinite(row_error) && row_error<=contract["validation_policy"]["row_native_max_scaled"] || error("STAR common-row/native mismatch")
            tag=construction==:batched_probes ? :STAR_common_grid_batched_probes : :STAR_common_grid_independent
            data=(copy(native.grid),qs,[q2],model.Q02,copy(native.tau),:transversity,5,Int(intervals),tag,
                blocks,sum(sizeof(b.rows) for b in blocks),row_error)
            verified_hoppet_build(buildfile)==verification.receipt || error("native build changed during STAR beam preparation")
            HoppetPointPlan(data...,hoppet_point_snapshot(data...))
        finally
            close(native)
        end
    end
end

function star_provenance(provenance)
    provenance isa AbstractDict || error("STAR provenance must be a dictionary")
    get(provenance,"scope","reconstructed_theory")=="reconstructed_theory" &&
        get(provenance,"production_ready",false)===false && get(provenance,"physical_extraction",false)===false ||
        error("reconstructed STAR component cannot consume production or extraction claims")
    meta=Dict{String,Any}(deepcopy(provenance))
    merge!(meta,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "contract_id"=>"STAR-EVOLVED-C11-001","historical_author_code_recovered"=>false,
        "finite_R_matching_implemented"=>false,"full_twist3_evolution"=>false))
    meta
end

"One fixed fragmentation leg, parameter-independent FFs and detached Hhat3 rows."
struct StarFragmentationPlan
    model::KPSYCSS
    radial::RadialRule
    Q2::Float64
    z::Float64
    jT::Float64
    pids::Tuple{Vararg{Int}}
    fragmentation::Symbol
    gluon_model::Symbol
    charge::Int
    matched::Matrix{Float64}
    operator::HoppetMatchingPlan
    samples::PolarizedInitialSamples
    uu_weights::Vector{Float64}
    collins_weights::Vector{Float64}
    quark_exponent::Vector{Float64}
    gluon_exponent::Vector{Float64}
    shape::Vector{Float64}
    unpolarized::Vector{Float64}
    provenance::Dict{String,Any}
    signature::String
end

function star_fragmentation_snapshot(model,radial,Q2,z,jT,pids,fragmentation,gluon_model,charge,matched,
        operator,samples,uu,pol,qe,ge,shape,unpolarized,provenance)
    numerical_snapshot((star_model_payload(model),star_radial_payload(radial),Q2,z,jT,pids,fragmentation,
        gluon_model,charge,matched,hoppet_point_signature(operator),star_sample_payload(samples),
        uu,pol,qe,ge,shape,unpolarized,provenance))
end
star_fragmentation_signature(p::StarFragmentationPlan)=star_fragmentation_snapshot(
    p.model,p.radial,p.Q2,p.z,p.jT,p.pids,p.fragmentation,p.gluon_model,p.charge,p.matched,p.operator,
    p.samples,p.uu_weights,p.collins_weights,p.quark_exponent,p.gluon_exponent,p.shape,p.unpolarized,p.provenance)
function star_fragmentation_valid(p::StarFragmentationPlan)
    hoppet_point_signature(p.operator)==p.operator.signature || error("mutated STAR Hhat3 operator")
    numerical_snapshot(star_sample_payload(p.samples))==p.samples.signature || error("mutated STAR Hhat3 samples")
    star_fragmentation_signature(p)==p.signature || error("mutated STAR fragmentation context")
    nothing
end

function star_fragmentation_plan(model::KPSYCSS,radial::RadialRule,Q2,z,jT,pids,matched,
        operator::HoppetMatchingPlan,samples::PolarizedInitialSamples;
        fragmentation,gluon_model,charge,provenance=Dict{String,Any}())
    initial_model_settings(model)
    ids=star_species(pids,fragmentation,gluon_model);q2=positive_fixed(Q2,"STAR Q2")
    zz=css_z(z);jt=fixed_scalar(jT,"STAR jT");jt>=0 || error("nonnegative STAR jT required")
    c=pion_charge(charge);m=deepcopy(model);r=deepcopy(radial);D=Matrix{Float64}(matched)
    size(D)==(length(ids),length(r.b)) && all(isfinite,D) || error("STAR matched FF must be flavor by radial node")
    fragmentation==:quark_only && !all(iszero,D[findfirst(==(21),ids),:]) && error("quark-only outgoing gluon leg must be zero")
    scales=css_mu_b2.(Ref(m),r.b)
    polarized_operator_valid(operator,samples,m,scales,:STAR_C11,:KPSY_Hhat3_homogeneous,[zz])
    hoppet_point_signature(operator)==operator.signature && numerical_snapshot(star_sample_payload(samples))==samples.signature ||
        error("invalid retained STAR fragmentation inputs")
    uu=hankel_weights(r,0,jt/zz)./zz^2;pol=-hankel_weights(r,1,jt/zz)./(2zz^3)
    qe=[(sudakov_perturbative(m,b,q2)+css_np_common(m,b,q2))/2 for b in r.b]
    ge=fragmentation==:quark_and_gluon ?
        [star_gluon_sudakov(m,b,q2)+(9/8)*css_np_common(m,b,q2) for b in r.b] : zeros(length(r.b))
    shape=r.b.^2 ./ zz^2
    U=[finite_transform(uu,view(D,f,:).*exp.(-(pid==21 ? ge : qe).-m.gh.*shape)) for (f,pid) in enumerate(ids)]
    all(isfinite,U) || error("nonfinite STAR unpolarized transform")
    meta=star_provenance(provenance)
    meta["initial_source_kind"]=string(samples.input_kind)
    signature=star_fragmentation_snapshot(m,r,q2,zz,jt,ids,fragmentation,gluon_model,c,D,operator,samples,uu,pol,qe,ge,shape,U,meta)
    StarFragmentationPlan(m,r,q2,zz,jt,ids,fragmentation,gluon_model,c,D,operator,samples,uu,pol,qe,ge,shape,U,meta,signature)
end

function star_collins(p::StarFragmentationPlan,theta::KPSYParameters)
    star_fragmentation_valid(p)
    theta.gh==p.model.gh || error("STAR initial/CSS width mismatch")
    evolved=apply_operator(p.operator,polarized_initial_columns(p.samples,theta))
    bases=collins_flavor_columns(evolved[1,:,:],p.charge,p.pids)
    damping=exp.(-p.quark_exponent.-collins_width(p.model,theta.values[13]).*p.shape)
    [finite_transform(p.collins_weights,view(bases,f,:).*damping) for f in eachindex(p.pids)]
end

"Fixed tree/PDF contractions for a complete differential STAR point."
struct StarPointPlan
    born::Tuple{Vararg{JetBornPlan}}
    weights::Vector{Float64}
    couplings::Vector{Float64}
    fA::Matrix{Float64}
    fB::Matrix{Float64}
    operator::HoppetPointPlan
    samples::PolarizedInitialSamples
    fragmentation::StarFragmentationPlan
    unpolarized_coefficients::Vector{Float64}
    transfer_coefficients::Array{Float64,3}
    jet_radius::Float64
    axis::Symbol
    scale_prescription::Symbol
    measure_source::String
    provenance::Dict{String,Any}
    signature::String
end

star_born_payload(b)=(b.pids,b.fragmentation,
    Tuple(getfield(b.kinematics,n) for n in fieldnames(JetBornKinematics)),
    Tuple((t.a,t.b,t.c,t.unpolarized_kind,t.transfer_kind,t.unpolarized,t.transfer) for t in b.terms))
function star_point_snapshot(born,weights,couplings,fA,fB,op,samples,ff,U,T,R,axis,scale,source,meta)
    numerical_snapshot((Tuple(star_born_payload(b) for b in born),weights,couplings,fA,fB,
        hoppet_point_signature(op),star_sample_payload(samples),star_fragmentation_signature(ff),U,T,R,axis,scale,source,meta))
end
star_point_signature(p::StarPointPlan)=star_point_snapshot(p.born,p.weights,p.couplings,p.fA,p.fB,p.operator,
    p.samples,p.fragmentation,p.unpolarized_coefficients,p.transfer_coefficients,p.jet_radius,p.axis,
    p.scale_prescription,p.measure_source,p.provenance)
function star_point_valid(p::StarPointPlan)
    star_fragmentation_valid(p.fragmentation)
    hoppet_point_signature(p.operator)==p.operator.signature || error("mutated STAR beam operator")
    numerical_snapshot(star_sample_payload(p.samples))==p.samples.signature || error("mutated STAR beam samples")
    star_point_signature(p)==p.signature || error("mutated STAR point context")
    nothing
end

function star_point_plan(born,weights,couplings,fA,fB,operator::HoppetPointPlan,samples::PolarizedInitialSamples,
        ff::StarFragmentationPlan;jet_radius,axis,scale_prescription,measure_source,provenance=Dict{String,Any}())
    star_fragmentation_valid(ff);n=length(born);nf=length(ff.pids)
    n>0 && all(b->b isa JetBornPlan,born) || error("nonempty Born-node collection required")
    k=first(born).kinematics
    all(b->b.pids==ff.pids && b.fragmentation==ff.fragmentation &&
        (b.kinematics.sqrt_s,b.kinematics.pT,b.kinematics.eta)==(k.sqrt_s,k.pT,k.eta),born) ||
        error("STAR point must share fixed jet kinematics, flavors and fragmentation choice")
    ff.Q2==k.pT^2 || error("historical STAR point requires Q=pT")
    R=positive_fixed(jet_radius,"jet radius");axis==:standard || error("only standard jet axis implemented")
    scale_prescription==:historical_pT_C11 || error("explicit historical_pT_C11 scale prescription required")
    measure_source isa AbstractString && !isempty(strip(measure_source)) || error("explicit STAR xa measure source required")
    w=[fixed_scalar(x,"STAR xa weight") for x in weights];a=[positive_fixed(x,"STAR hard alpha_s") for x in couplings]
    length(w)==length(a)==n && all(>=(0),w) && sum(w)>0 || error("invalid STAR xa measure or hard couplings")
    all(==(alpha_s(ff.model.coupling,k.pT^2)),a) || error("STAR hard coupling must match the fixed physical model at pT")
    A,B=Matrix{Float64}(fA),Matrix{Float64}(fB)
    size(A)==size(B)==(nf,n) && all(isfinite,A) && all(isfinite,B) || error("finite flavor-by-xa PDFs required")
    operator.quantity==:transversity && operator.nf==5 && operator.Q02==ff.model.Q02 &&
        operator.Q2==[k.pT^2] && operator.queries==[b.kinematics.x_a for b in born] &&
        operator.tau==evolution_time.(Ref(ff.model.coupling),[k.pT^2],ff.model.Q02) &&
        operator.grid==samples.grid && samples.kind==:transversity || error("STAR collinear beam operator/sample mismatch")
    operator.construction in (:STAR_common_grid_batched_probes,:STAR_common_grid_independent,:synthetic_fixture) ||
        error("STAR beam requires the separately qualified common native grid")
    hoppet_point_signature(operator)==operator.signature && numerical_snapshot(star_sample_payload(samples))==samples.signature ||
        error("invalid retained STAR beam inputs")
    U=zeros(nf);T=zeros(nf,nf,n)
    for (i,b) in enumerate(born)
        factor=w[i]*jet_born_prefactor(b.kinematics,a[i])
        for t in b.terms
            U[t.c]+=factor*t.unpolarized*A[t.a,i]*B[t.b,i]
            t.transfer_kind==:zero && continue
            T[t.a,t.c,i]+=factor*t.transfer*B[t.b,i]
        end
    end
    samples.input_kind==ff.samples.input_kind || error("STAR initial source kinds must agree")
    meta=star_provenance(provenance);nodes=Tuple(born);source=String(measure_source)
    meta["initial_source_kind"]=string(samples.input_kind)
    meta["beam_row_construction"]=string(operator.construction)
    meta["uniform_near_endpoint_h1_accuracy_certified"]=false
    signature=star_point_snapshot(nodes,w,a,A,B,operator,samples,ff,U,T,R,axis,scale_prescription,source,meta)
    StarPointPlan(nodes,w,a,A,B,operator,samples,ff,U,T,R,axis,scale_prescription,source,meta,signature)
end

"Only generic Julia arithmetic and detached row products occur during evaluation."
function star_structure(p::StarPointPlan,theta::KPSYParameters)
    star_point_valid(p);ff=p.fragmentation
    theta.gh==ff.model.gh || error("STAR initial/CSS width mismatch")
    h=apply_operator(p.operator,polarized_initial_columns(p.samples,theta))[:,:,1]
    C=star_collins(ff,theta);T=zero(theta.values[1]);ids=ff.pids
    for (i,a) in enumerate(ids)
        a in (2,1) || continue
        column=a==2 ? 1 : 2
        for c in eachindex(ids),j in axes(h,1)
            T+=p.transfer_coefficients[i,c,j]*h[j,column]*C[c]
        end
    end
    U=dot(p.unpolarized_coefficients,ff.unpolarized)
    isfinite(U) && isfinite(T) || error("nonfinite evolved STAR structure")
    (UU=U,UT=T)
end
star_structure(p::StarPointPlan,theta::AbstractVector)=star_structure(p,KPSYParameters(theta;gh=p.fragmentation.model.gh))
predict_star(p::StarPointPlan,theta;phase)=jet_born_asymmetry(star_structure(p,theta);phase=phase)

function with_momentum(p::StarPointPlan,jT)
    star_point_valid(p);f=p.fragmentation
    ff=star_fragmentation_plan(f.model,f.radial,f.Q2,f.z,jT,f.pids,f.matched,f.operator,f.samples;
        fragmentation=f.fragmentation,gluon_model=f.gluon_model,charge=f.charge,provenance=f.provenance)
    star_point_plan(p.born,p.weights,p.couplings,p.fA,p.fB,p.operator,p.samples,ff;
        jet_radius=p.jet_radius,axis=p.axis,scale_prescription=p.scale_prescription,measure_source=p.measure_source,provenance=p.provenance)
end

struct StarBinPlan
    points::Tuple{Vararg{StarPointPlan}}
    weights::Tuple{Vararg{Float64}}
    phase::Symbol
    measure_source::String
    signature::String
end
star_bin_signature(p::StarBinPlan)=numerical_snapshot((Tuple(star_point_signature(x) for x in p.points),p.weights,p.phase,p.measure_source))
function star_bin(points,weights;phase,measure_source)
    !isempty(points) && all(x->x isa StarPointPlan,points) || error("nonempty STAR point collection required")
    foreach(star_point_valid,points)
    phase in (:sin_phiS_minus_phiH,:sin_phiH_minus_phiS) || error("explicit supported STAR phase required")
    measure_source isa AbstractString && !isempty(strip(measure_source)) || error("explicit STAR bin measure source required")
    w=Tuple(fixed_scalar(x,"STAR bin weight") for x in weights)
    length(w)==length(points) && all(>=(0),w) && sum(w)>0 || error("nonnegative STAR bin weights with positive sum required")
    firstp=first(points);f=firstp.fragmentation
    all(p->p.fragmentation.pids==f.pids && p.fragmentation.fragmentation==f.fragmentation &&
        p.fragmentation.gluon_model==f.gluon_model && p.fragmentation.charge==f.charge &&
        star_model_payload(p.fragmentation.model)==star_model_payload(f.model) &&
        (p.axis,p.scale_prescription,p.jet_radius,first(p.born).kinematics.sqrt_s)==
        (firstp.axis,firstp.scale_prescription,firstp.jet_radius,first(firstp.born).kinematics.sqrt_s),points) ||
        error("STAR bin mixes incompatible physical model choices")
    ps=Tuple(points);source=String(measure_source)
    signature=numerical_snapshot((Tuple(star_point_signature(x) for x in ps),w,phase,source))
    StarBinPlan(ps,w,phase,source,signature)
end
function star_structure(p::StarBinPlan,theta)
    foreach(star_point_valid,p.points)
    star_bin_signature(p)==p.signature || error("mutated STAR bin context")
    values=map(x->star_structure(x,theta),p.points)
    (UU=sum(p.weights[i]*values[i].UU for i in eachindex(values)),UT=sum(p.weights[i]*values[i].UT for i in eachindex(values)))
end
predict_star(p::StarBinPlan,theta)=jet_born_asymmetry(star_structure(p,theta);phase=p.phase)
function star_momentum_bin(p::StarPointPlan,rule::MomentumRule;phase,measure,measure_source)
    momentum_valid(rule)
    measure_source isa AbstractString && !isempty(strip(measure_source)) || error("explicit STAR momentum measure source required")
    measure in (:jT_djT,:d2jT) || error("explicit jT_djT or d2jT measure required")
    factor=measure==:d2jT ? 2pi : 1.
    star_bin([with_momentum(p,j) for j in rule.nodes],factor.*rule.weights;phase=phase,
        measure_source="$(measure_source); $(measure) on [$(rule.lower),$(rule.upper)]")
end

function matched_star_ff_inputs(ff::LegacyPartonInput,model,radial,z,charge,pids,order,policy;fragmentation,gluon_model)
    ids=star_species(pids,fragmentation,gluon_model);quarks=filter(!=(21),collect(ids))
    output=zeros(length(ids),length(radial.b));qi=[findfirst(==(p),ids) for p in quarks]
    output[qi,:]=matched_ff_inputs(ff,:STAR_C11,model,radial,z,charge,quarks,order,policy)
    if fragmentation==:quark_and_gluon
        g=findfirst(==(21),ids)
        for (i,b) in enumerate(radial.b)
            q2=css_mu_b2(model,b);rule=star_gluon_matching_rule(z,alpha_s(model.coupling,q2);order=order,knots=ff.x)
            values=hcat([pion_ff(ff,collect(ids),charge,y,q2;policy=policy) for y in rule.coordinates]...)
            output[g,i]=dot(rule.gluon_weights,view(values,g,:))+dot(rule.singlet_weights,vec(sum(values[qi,:];dims=1)))
        end
    end
    output
end

"All provider/native preparation is performed once, outside parameter evaluation."
function prepare_star_point(pdf::CT10Input,hel::LegacyPartonInput,ff::LegacyPartonInput,buildfile,
        model::KPSYCSS,radial::RadialRule,sqrt_s,pT,eta,z,jT,charge;
        fragmentation,gluon_model,jet_radius,axis,scale_prescription,outer_order=64,matching_order=8,
        intervals=1024,range_policy=:strict,construction=:batched_probes)
    initial_model_settings(model);ids=(INITIAL_FLAVORS...,21);star_species(ids,fragmentation,gluon_model)
    pdf.policy==range_policy || error("CT10 and requested range policies must agree explicitly")
    pt=positive_fixed(pT,"jet pT");xmin=jet_born_xmin(sqrt_s,pt,eta)
    nodes,weights=gauss_rule(outer_order);xa=(1+xmin)/2 .+(1-xmin)/2 .* nodes;w=(1-xmin)/2 .* weights
    born=[JetBornPlan(ids,JetBornKinematics(sqrt_s,pt,eta,x);fragmentation=fragmentation) for x in xa]
    xb=[b.kinematics.x_b for b in born];q2=fill(pt^2,length(xa))
    A=ct10_batch(pdf,collect(ids),xa,q2)./permutedims(xa)
    B=ct10_batch(pdf,collect(ids),xb,q2)./permutedims(xb)
    hop=star_beam_operator(buildfile,xa,pt^2,model;intervals=intervals,construction=construction)
    x,bases=prepare_initial_soffer(pdf,hel,hop.grid)
    hs=PolarizedInitialSamples(:transversity,x,bases;input_kind=:pinned_inputs)
    scales=css_mu_b2.(Ref(model),radial.b)
    H=star_native_checked(()->HoppetMatchingPlan(buildfile,[z],scales,model.Q02,model.coupling;
        process=:STAR_C11,intervals=intervals,construction=construction))
    zs,ffs=prepare_initial_ff(ff,H.grid)
    Hs=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,zs,ffs;input_kind=:pinned_inputs)
    D=matched_star_ff_inputs(ff,model,radial,z,charge,ids,matching_order,range_policy;fragmentation=fragmentation,gluon_model=gluon_model)
    audit=Dict("CT10_beam_nodes_outside_strict"=>count(t->!parton_domain(pdf,t,pt^2)["within_grid"],vcat(xa,xb)),
        "DSS14_radial_nodes_outside_strict"=>count(q->!parton_domain(ff,z,q;policy=range_policy)["within_grid"],scales),
        "radial_nodes"=>length(scales),"beam_nodes"=>2length(xa))
    meta=prepared_input_metadata(:STAR_C11,model,radial,matching_order,range_policy,audit)
    meta["outer_order"]=outer_order;meta["input_kind"]="pinned_CT10_DSSV08_DSS14"
    f=star_fragmentation_plan(model,radial,pt^2,z,jT,ids,D,H,Hs;
        fragmentation=fragmentation,gluon_model=gluon_model,charge=charge,provenance=meta)
    star_point_plan(born,w,fill(alpha_s(model.coupling,pt^2),length(xa)),A,B,hop,hs,f;
        jet_radius=jet_radius,axis=axis,scale_prescription=scale_prescription,
        measure_source="Gauss-Legendre dxa on the exact open physical interval; no jet or jT bin measure",provenance=meta)
end
