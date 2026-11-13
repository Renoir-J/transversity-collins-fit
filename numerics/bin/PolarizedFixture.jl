# Analytic polynomial initial inputs for reconstructed-theory checks only.
# Include PolarizedReference.jl first; this helper is not a physical provider.
function analytic_context(buildfile,process,model,radial;intervals=1024,target=NucleonTarget(:proton))
    q2=css_mu_b2.(Ref(model),radial.b)
    Q2=process==:SIDIS ? 4. : 110.
    plan=leptonic_plan(process,model,radial,Q2,.3,.3;z2=process==:SIA ? .5 : nothing)
    unpol1=[1+.1pid for pid in plan.flavors,b in radial.b]
    unpol2=[.8-.05pid for pid in plan.flavors,b in radial.b]
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    op=HoppetMatchingPlan(buildfile,process==:SIDIS ? [.1] : [.3,.5],q2,2.4,model.coupling;
        process=process,quantity=kind,intervals=intervals)
    bases=process==:SIDIS ? [(i==1 ? 1.3 : .8)*(1-x) for i in 1:2,x in op.grid] :
        [PolarizedReference.FF_BASES[i]*(1-x)^3 for i in 1:4,x in op.grid]
    samples=PolarizedInitialSamples(kind,op.grid,bases;input_kind=:analytic)
    if process==:SIDIS
        op2=HoppetMatchingPlan(buildfile,[.3],q2,2.4,model.coupling;process=process,intervals=intervals)
        samples2=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,
            [PolarizedReference.FF_BASES[i]*(1-x)^3 for i in 1:4,x in op2.grid];input_kind=:analytic)
        polarized_leptonic_plan(plan,unpol1,unpol2,op,samples;second_operator=op2,second_samples=samples2,target=target,charges=(1,0))
    else
        polarized_leptonic_plan(plan,unpol1,unpol2,op,samples;charges=(1,-1))
    end
end
