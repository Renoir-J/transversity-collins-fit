# These are deliberately synthetic row fixtures. They check assembly and generic
# differentiation, not HOPPET accuracy or physical process admission.
function synthetic_polarized_operator(process,kind,queries,radial,model)
    grid=sort!(unique(vcat(.05,.1,queries,.75,1.)))
    q2=css_mu_b2.(Ref(model),radial.b)
    tau=evolution_time.(Ref(model.coupling),q2,model.Q02)
    alpha=alpha_s.(Ref(model.coupling),q2)
    blocks=CollinsFit.HoppetPointBlock[]
    for query in queries
        rows=zeros(length(q2),length(grid))
        rows[:,findfirst(==(query),grid)].=1 .+ 0.01 .* radial.b
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),rows))
    end
    HoppetMatchingPlan(grid,collect(queries),q2,model.Q02,tau,kind,5,32,:synthetic_fixture,
        blocks,sum(sizeof(b.rows) for b in blocks),0.,process,alpha,
        hash((grid,collect(queries),q2,tau,process,alpha)))
end

@testset "cached initial columns preserve the physical parameter family" begin
    p=kpsy_published_parameters();theta=collect(p.values)
    grid=[.05,.2,.4,.9,1.]
    soffer=[(1+.2i)*x*(1-x)^2 for i in 1:2,x in grid]
    ff=[(1+.3i)*x*(1-x)^3 for i in 1:4,x in grid]
    hs=PolarizedInitialSamples(:transversity,grid,soffer;input_kind=:analytic)
    Hs=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,grid,ff;input_kind=:analytic)
    reference_h(t)=permutedims([initial_transversity(KPSYParameters(t),pid,x,soffer[i,j]) for (i,pid) in enumerate((2,1)),(j,x) in enumerate(grid)])
    reference_H(t)=permutedims([CollinsFit.initial_shape(x,t[i==1 ? 7 : 8],t[i==1 ? 9 : 10],t[i==1 ? 11 : 12],ff[i,j]) for i in 1:4,(j,x) in enumerate(grid)])
    for (samples,reference) in ((hs,reference_h),(Hs,reference_H))
        actual(t)=vec(polarized_initial_columns(samples,KPSYParameters(t)))
        @test actual(theta) ≈ vec(reference(theta)) rtol=1e-14
        @test ForwardDiff.jacobian(actual,theta) ≈ ForwardDiff.jacobian(t->vec(reference(t)),theta) rtol=1e-13 atol=1e-15
        @test ForwardDiff.hessian(t->sum(abs2,actual(t)),theta) ≈ ForwardDiff.hessian(t->sum(abs2,reference(t)),theta) rtol=1e-12 atol=1e-14
        @test all(iszero,polarized_initial_columns(samples,p)[end,:])
    end
    bases=permutedims(polarized_initial_columns(Hs,p))
    for charge in (-1,0,1)
        actual=CollinsFit.collins_flavor_columns(bases,charge,CollinsFit.INITIAL_FLAVORS)
        expected=[initial_collins(p,pid,charge,x,ff[:,j]) for pid in CollinsFit.INITIAL_FLAVORS,(j,x) in enumerate(grid)]
        @test actual ≈ expected rtol=1e-14
    end
    singular=copy(theta);singular[3]=0.
    @test polarized_initial_columns(hs,KPSYParameters(singular)) ≈ reference_h(singular)
    @test_throws ErrorException ForwardDiff.jacobian(t->vec(polarized_initial_columns(hs,KPSYParameters(t))),singular)
    hs.bases[1,1]+=1
    @test_throws ErrorException polarized_initial_columns(hs,p)
    @test_throws ErrorException PolarizedInitialSamples(:transversity,grid,-soffer;input_kind=:analytic)
    @test_throws ErrorException PolarizedInitialSamples(:transversity,grid,soffer;input_kind=:unknown)
end

@testset "synthetic process assembly, targets, charges and nested derivatives" begin
    model=KPSYCSS();radial=RadialRule([0.,.25,.75,2.];order=4)
    theta=collect(kpsy_published_parameters().values)
    first=[(1+.07i)*(1+.03b) for i in 1:10,b in radial.b]
    second=[(1+.11i)*(1+.05b) for i in 1:10,b in radial.b]
    anti=[findfirst(==(-pid),CollinsFit.INITIAL_FLAVORS) for pid in CollinsFit.INITIAL_FLAVORS]
    swap=[findfirst(==(abs(pid)==1 ? 2sign(pid) : abs(pid)==2 ? sign(pid) : pid),CollinsFit.INITIAL_FLAVORS) for pid in CollinsFit.INITIAL_FLAVORS]
    contexts=PolarizedLeptonicPlan[]
    for process in (:SIDIS,:SIA)
        plan=leptonic_plan(process,model,radial,process==:SIDIS ? 4. : 110.,.2,.3;z2=process==:SIA ? .5 : nothing)
        kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
        op=synthetic_polarized_operator(process,kind,process==:SIDIS ? [.15] : [.3,.5],radial,model)
        n=kind==:transversity ? 2 : 4
        samples=PolarizedInitialSamples(kind,op.grid,[(1+.3i)*x*(1-x)^3 for i in 1:n,x in op.grid];input_kind=:analytic)
        op2=process==:SIDIS ? synthetic_polarized_operator(process,:KPSY_Hhat3_homogeneous,[.3],radial,model) : nothing
        samples2=op2===nothing ? nothing : PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,[(1+.3i)*x*(1-x)^3 for i in 1:4,x in op2.grid];input_kind=:analytic)
        target=process==:SIDIS ? NucleonTarget("synthetic mixture",(.6,.4),(.5,-.2);source="algebra test fixture") : nothing
        context=polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=samples2,target=target,charges=process==:SIDIS ? (-1,0) : (-1,1))
        push!(contexts,context)
        @test !context.provenance["production_ready"] && !context.provenance["physical_extraction"]
        @test context.provenance["contract_id"]=="KPSY-POLARIZED-PROCESS-001"
        @test iszero(polarized_structure(with_momentum(context,0.),theta).Collins)
        @test with_momentum(context,.4).first_operator === op
        @test with_pion_charges(context,0).first_samples === samples
        @test with_momentum(context,.2).unpolarized ≈ context.unpolarized
        charge_matrix(m,c)=c==1 ? m : c==-1 ? m[anti,:] : (m+m[anti,:])/2
        f=process==:SIDIS ? .6first+.4first[swap,:] : charge_matrix(first,-1)
        g=charge_matrix(second,process==:SIDIS ? -1 : 1)
        @test context.unpolarized ≈ leptonic_unpolarized(plan,f,g)
        for charge in (-1,0,1)
            changed=with_pion_charges(context,charge)
            f=process==:SIDIS ? .6first+.4first[swap,:] : charge_matrix(first,charge)
            g=charge_matrix(second,process==:SIDIS ? charge : 1)
            @test changed.unpolarized ≈ leptonic_unpolarized(plan,f,g)
        end
        prediction(t)=polarized_structure(context,t).Collins
        if process==:SIDIS
            nucleons=[polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=samples2,
                target=NucleonTarget(nucleon),charges=(-1,0)) for nucleon in (:proton,:neutron)]
            @test prediction(theta) ≈ .5polarized_structure(nucleons[1],theta).Collins-.2polarized_structure(nucleons[2],theta).Collins
            @test context.unpolarized ≈ .6nucleons[1].unpolarized+.4nucleons[2].unpolarized
            p=KPSYParameters(theta)
            ix=findfirst(==(.15),samples.grid);iz=findfirst(==(.3),samples2.grid)
            hu=initial_transversity(p,2,.15,samples.bases[1,ix])
            hd=initial_transversity(p,1,.15,samples.bases[2,ix])
            Hu=initial_collins(p,2,-1,.3,samples2.bases[:,iz])
            Hd=initial_collins(p,1,-1,.3,samples2.bases[:,iz])
            # Scalar u/d isospin table, without the implementation's flavor map.
            scalar_transform=sum(-radial.db[j]*b^2*CollinsFit.besselj(1,.2b/.3)*
                (1+.01b)^2*exp(-sudakov_perturbative(model,b,4.)-sidis_collins_snp(model,b,4.,.3,p.values[13]))/
                (4pi*.3^3) for (j,b) in enumerate(radial.b))
            @test prediction(theta) ≈ scalar_transform*((4/9)*(.5hu-.2hd)*Hu+(1/9)*(.5hd-.2hu)*Hd)
            @test polarized_structure(nucleons[2],theta).Collins ≈ scalar_transform*((4/9)*hd*Hu+(1/9)*hu*Hd)
        end
        # BigFloat-backed Duals take the independent generic row path, while
        # Float64 primal/tangent columns use batched BLAS. Compare all 13 columns.
        big=BigFloat.(theta)
        @test ForwardDiff.gradient(prediction,theta) ≈ Float64.(ForwardDiff.gradient(prediction,big)) rtol=1e-12 atol=1e-14
        @test ForwardDiff.hessian(prediction,theta) ≈ Float64.(ForwardDiff.hessian(prediction,big)) rtol=1e-11 atol=1e-13
        @test_throws ErrorException polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=samples2,target=target,provenance=Dict("physical_extraction"=>true))
        saved=context.unpolarized_first[1,1];context.unpolarized_first[1,1]+=1
        @test_throws ErrorException polarized_structure(context,theta)
        context.unpolarized_first[1,1]=saved
        # The subsequent contexts must not depend on native global state.
    end
    @test NucleonTarget(:neutron).unpolarized==(0.,1.)
    @test_throws ErrorException NucleonTarget(:deuteron)
    @test_throws ErrorException NucleonTarget("bad",(1.,0.),(1.1,0.);source="test")
    @test_throws ErrorException NucleonTarget("bad",(1.,0.),(NaN,0.),"test")
    @test_throws ErrorException with_pion_charges(contexts[1],1,1)
end
