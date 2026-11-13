# Explicit analytic provider functions and native retained operators. These are
# reconstruction fixtures, not fitted distributions or experiment admission.
const STAR_IDS=(CollinsFit.INITIAL_FLAVORS...,21)
star_edges(B)=vcat([0.,.001,.005,.02,.1,.4,1.,2.,4.,8.,12.,16.],B==24 ? [24.] : [24.,B])
star_gluon_model(fragmentation)=fragmentation==:quark_only ? :none : :casimir_common_same_intrinsic_1705
function analytic_star_fragmentation(buildfile,model,radial,pt,z,jT,charge,fragmentation;intervals=1024)
    scales=css_mu_b2.(Ref(model),radial.b)
    H=CollinsFit.star_native_checked(()->HoppetMatchingPlan(buildfile,[z],scales,model.Q02,model.coupling;process=:STAR_C11,intervals=intervals))
    Hs=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,H.grid,
        [StarEvolvedReference.analytic_ff_bases[i]*(1-x)^3 for i in 1:4,x in H.grid];input_kind=:analytic)
    D=zeros(11,length(radial.b));qi=findall(!=(21),STAR_IDS);g=findfirst(==(21),STAR_IDS)
    for (j,q2) in enumerate(scales)
        rule=matching_rule(:STAR_C11,:ff,z,alpha_s(model.coupling,q2);order=32)
        for i in qi
            D[i,j]=match_leg(rule,StarEvolvedReference.analytic_ff.(STAR_IDS[i],rule.coordinates);
                gluon=StarEvolvedReference.analytic_ff.(21,rule.coordinates))
        end
        if fragmentation==:quark_and_gluon
            gr=CollinsFit.star_gluon_matching_rule(z,alpha_s(model.coupling,q2);order=32)
            D[g,j]=dot(gr.gluon_weights,StarEvolvedReference.analytic_ff.(21,gr.coordinates))+
                dot(gr.singlet_weights,[sum(StarEvolvedReference.analytic_ff(pid,x) for pid in STAR_IDS if pid!=21) for x in gr.coordinates])
        end
    end
    star_fragmentation_plan(model,radial,pt^2,z,jT,STAR_IDS,D,H,Hs;fragmentation=fragmentation,
        gluon_model=star_gluon_model(fragmentation),charge=charge,provenance=Dict("input_kind"=>"analytic_polynomial"))
end
function star_on_outer_grid(buildfile,f,sqrt_s,pt,eta;outer_order=64,intervals=1024,pdf=nothing,hel=nothing)
    model=f.model;xmin=jet_born_xmin(sqrt_s,pt,eta);nodes,weights=CollinsFit.gauss_rule(outer_order)
    xa=(1+xmin)/2 .+(1-xmin)/2 .* nodes;w=(1-xmin)/2 .* weights
    born=[JetBornPlan(f.pids,JetBornKinematics(sqrt_s,pt,eta,x);fragmentation=f.fragmentation) for x in xa]
    xb=[b.kinematics.x_b for b in born]
    h=CollinsFit.star_beam_operator(buildfile,xa,pt^2,model;intervals=intervals)
    if pdf===nothing
        bases=[(i==1 ? 1.3 : .8)*(1-x) for i in 1:2,x in h.grid]
        kind=:analytic
        A=[StarEvolvedReference.analytic_pdf(pid,x) for pid in f.pids,x in xa]
        B=[StarEvolvedReference.analytic_pdf(pid,x) for pid in f.pids,x in xb]
    else
        _,bases=CollinsFit.prepare_initial_soffer(pdf,hel,h.grid);kind=:pinned_inputs
        A=ct10_batch(pdf,collect(f.pids),xa,fill(pt^2,length(xa)))./permutedims(xa)
        B=ct10_batch(pdf,collect(f.pids),xb,fill(pt^2,length(xb)))./permutedims(xb)
    end
    hs=PolarizedInitialSamples(:transversity,h.grid,bases;input_kind=kind)
    meta=deepcopy(f.provenance);meta["outer_order"]=outer_order
    if pdf!==nothing
        audit=deepcopy(get(meta,"source_range_audits",Dict{String,Any}()))
        audit["beam_nodes"]=2length(xa)
        audit["CT10_beam_nodes_outside_strict"]=count(x->!parton_domain(pdf,x,pt^2)["within_grid"],vcat(xa,xb))
        meta["source_range_audits"]=audit
    end
    star_point_plan(born,w,fill(alpha_s(model.coupling,pt^2),length(xa)),A,B,h,hs,f;
        jet_radius=.6,axis=:standard,scale_prescription=:historical_pT_C11,
        measure_source="Gauss-Legendre dxa; exact massless physical interval; no measured acceptance",provenance=meta)
end
function analytic_star_context(buildfile,model,radial;pt=12.9,sqrt_s=200.,eta=.3,z=.3,jT=.3,charge=1,
        fragmentation=:quark_only,outer_order=64,intervals=1024)
    f=analytic_star_fragmentation(buildfile,model,radial,pt,z,jT,charge,fragmentation;intervals=intervals)
    star_on_outer_grid(buildfile,f,sqrt_s,pt,eta;outer_order=outer_order,intervals=intervals)
end
