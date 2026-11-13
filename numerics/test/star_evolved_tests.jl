using QuadGK
include("../bin/StarEvolvedReference.jl")

# Synthetic identity rows with explicit arbitrary b dependence. These tests
# exercise assembly/integrity/AD, not native evolution or physical data.
function synthetic_star_fragmentation_operator(radial,model)
    grid=[.05,.1,.3,.75,1.];q2=css_mu_b2.(Ref(model),radial.b)
    rows=zeros(length(q2),length(grid));rows[:,3].=1 .+ .01 .* radial.b
    blocks=[CollinsFit.HoppetPointBlock(collect(eachindex(grid)),rows)]
    data=(grid,[.3],q2,model.Q02,evolution_time.(Ref(model.coupling),q2,model.Q02),
        :KPSY_Hhat3_homogeneous,5,32,:synthetic_fixture,blocks,sizeof(first(blocks).rows),0.)
    alpha=alpha_s.(Ref(model.coupling),q2)
    HoppetMatchingPlan(data...,:STAR_C11,alpha,CollinsFit.hoppet_point_snapshot(data...;process=:STAR_C11,alpha=alpha))
end
function synthetic_star_point(;pt=20.,charge=1,fragmentation=:quark_only,jT=.25)
    model=KPSYCSS();radial=RadialRule([0.,.1,.5,1.,3.,8.];order=4)
    ids=(CollinsFit.INITIAL_FLAVORS...,21);xmin=jet_born_xmin(200.,pt,.3)
    xa=xmin .+(1-xmin).*[.15,.4,.75]
    born=[JetBornPlan(ids,JetBornKinematics(200.,pt,.3,x);fragmentation=fragmentation) for x in xa]
    H=synthetic_star_fragmentation_operator(radial,model)
    Hs=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,H.grid,
        [(1+.2i)*(1-x)^3 for i in 1:4,x in H.grid];input_kind=:analytic)
    D=[(pid==21 && fragmentation==:quark_only) ? 0. : (1+.02i)*(1+.01b) for (i,pid) in enumerate(ids),b in radial.b]
    gluon_model=fragmentation==:quark_only ? :none : :casimir_common_same_intrinsic_1705
    ff=star_fragmentation_plan(model,radial,pt^2,.3,jT,ids,D,H,Hs;fragmentation=fragmentation,gluon_model=gluon_model,charge=charge)
    grid=sort!(unique(vcat(.05,xa,1.)));blocks=CollinsFit.HoppetPointBlock[]
    for x in xa
        row=zeros(1,length(grid));row[1,findfirst(==(x),grid)]=1+.02x
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),row))
    end
    data=(grid,xa,[pt^2],model.Q02,evolution_time.(Ref(model.coupling),[pt^2],model.Q02),
        :transversity,5,32,:synthetic_fixture,blocks,sum(sizeof(b.rows) for b in blocks),0.)
    h=HoppetPointPlan(data...,CollinsFit.hoppet_point_snapshot(data...))
    hs=PolarizedInitialSamples(:transversity,grid,[(1+.3i)*(1-x)^2 for i in 1:2,x in grid];input_kind=:analytic)
    A=[(1+.1i)*(1-x)^3 for i in 1:11,x in xa];B=[(1+.07i)*(1-b.kinematics.x_b)^3 for i in 1:11,b in born]
    star_point_plan(born,[.2,.3,.4],fill(alpha_s(model.coupling,pt^2),3),A,B,h,hs,ff;
        jet_radius=.6,axis=:standard,scale_prescription=:historical_pT_C11,measure_source="synthetic three-node dxa fixture")
end

@testset "C11 single fragmentation matching, gluon singlet and exponents" begin
    model=KPSYCSS()
    @test CollinsFit.star_gluon_regular(1.) == (-6.,0.)
    @test CollinsFit.matching_delta(:STAR_C11,:ff)==0.
    @test CollinsFit.matching_delta(:STAR_C11,:collins)==0.
    @test CollinsFit.matching_delta(:SIDIS,:ff)==-8/3
    @test CollinsFit.matching_delta(:SIA,:collins)==(pi^2-8)/3
    for channel in (:pdf,:transversity)
        @test_throws ErrorException matching_rule(:STAR_C11,channel,.3,.2)
    end
    for z in (.05,.3,.8),alpha in (0.,.2)
        q(x)=1.2*(1-x)^3;g(x)=.7*(1-x)^2;singlet(x)=4.3*(1-x)^3
        rule=matching_rule(:STAR_C11,:ff,z,alpha;order=32)
        @test match_leg(rule,q.(rule.coordinates);gluon=g.(rule.coordinates)) ≈
            StarEvolvedReference.reference_matching(z,alpha,:quark,q,g) rtol=1e-11
        gr=CollinsFit.star_gluon_matching_rule(z,alpha;order=32)
        actual=dot(gr.gluon_weights,g.(gr.coordinates))+dot(gr.singlet_weights,singlet.(gr.coordinates))
        @test actual ≈ StarEvolvedReference.reference_matching(z,alpha,:gluon,q,g;singlet=singlet) rtol=1e-11
    end
    for b in (.03,.3,3.),Q2 in (4.,400.,961.)
        @test sudakov_perturbative(model,b,Q2)/2 ≈ StarEvolvedReference.reference_single_sudakov(model,b,Q2) rtol=1e-11 atol=1e-13
        @test CollinsFit.star_gluon_sudakov(model,b,Q2) ≈ StarEvolvedReference.reference_single_sudakov(model,b,Q2;gluon=true) rtol=1e-11 atol=1e-13
    end
    @test_throws ErrorException CollinsFit.star_species((CollinsFit.INITIAL_FLAVORS...,21),:quark_and_gluon,:none)
    @test_throws ErrorException CollinsFit.star_species((2,-2,21),:quark_only,:none)
end

@testset "synthetic complete STAR scalar assembly, bin measures and full13 AD" begin
    t=collect(kpsy_published_parameters().values);t[2]=-.91;t[11]=.4;t[12]=.7
    for fragmentation in (:quark_only,:quark_and_gluon),charge in (-1,0,1)
        p=synthetic_star_point(;fragmentation=fragmentation,charge=charge)
        reference=StarEvolvedReference.scalar_reference(p)
        actual=star_structure(p,t);expected=reference(t)
        @test collect(actual) ≈ collect(expected) rtol=1e-11
        @test predict_star(p,t;phase=:sin_phiH_minus_phiS)==-predict_star(p,t;phase=:sin_phiS_minus_phiH)
        @test iszero(star_structure(with_momentum(p,0.),t).UT)
        @test with_momentum(p,.4).operator===p.operator
        @test with_momentum(p,.4).fragmentation.operator===p.fragmentation.operator
        @test all(iszero,star_collins(p.fragmentation,KPSYParameters(t))[findall(x->abs(x)>=4,p.fragmentation.pids)])
        f(x)=predict_star(p,x;phase=:sin_phiS_minus_phiH)
        r(x)=reference(x).UT/reference(x).UU
        @test ForwardDiff.gradient(f,t) ≈ ForwardDiff.gradient(r,t) rtol=1e-10 atol=1e-13
        @test ForwardDiff.hessian(f,t) ≈ ForwardDiff.hessian(r,t) rtol=1e-9 atol=1e-12
        @test all(abs.(ForwardDiff.gradient(f,t)).>1e-12)
        @test !p.provenance["production_ready"] && p.provenance["initial_source_kind"]=="analytic"
        high=synthetic_star_point(;fragmentation=fragmentation,charge=charge,pt=30.,jT=.5)
        bin=star_bin([p,high],[.3,1.7];phase=:sin_phiS_minus_phiH,measure_source="synthetic supplied jet/Jacobian/acceptance weights")
        hi=star_structure(high,t)
        target=(.3actual.UT+1.7hi.UT)/(.3actual.UU+1.7hi.UU)
        @test predict_star(bin,t) ≈ target rtol=1e-13
        @test abs(predict_star(bin,t)-(.3f(t)+1.7predict_star(high,t;phase=:sin_phiS_minus_phiH))/2)>1e-8
        momentum=MomentumRule(.05,.6;order=4)
        mb=star_momentum_bin(p,momentum;phase=:sin_phiS_minus_phiH,measure=:jT_djT,measure_source="synthetic finite jT interval")
        ab=star_momentum_bin(p,momentum;phase=:sin_phiS_minus_phiH,measure=:d2jT,measure_source="synthetic finite d2jT interval")
        @test collect(star_structure(ab,t)) ≈ 2pi.*collect(star_structure(mb,t)) rtol=1e-13
        @test abs(star_structure(mb,t).UT)>0
        @test predict_star(mb,t) ≈ predict_star(ab,t) rtol=1e-13
    end
end

@testset "STAR fixed cache mutation and model boundary rejection" begin
    t=collect(kpsy_published_parameters().values)
    mutations=[p->(p.fA[2,2]+=1),p->(p.fB[3,2]+=1),p->(p.weights[2]*=2),
        p->(p.couplings[2]*=2),p->(p.unpolarized_coefficients[2]+=1),p->(p.transfer_coefficients[2,2,2]+=1),
        p->(p.operator.grid[2]*=.99),p->(p.operator.blocks[2].indices[2]=1),
        p->(p.samples.bases[1,2]+=1),p->(p.samples.logarithm[2]+=1),
        p->(p.fragmentation.matched[2,2]+=1),p->(p.fragmentation.uu_weights[2]*=2),
        p->(p.fragmentation.collins_weights[2]*=2),p->(p.fragmentation.quark_exponent[2]+=1),
        p->(p.fragmentation.gluon_exponent[2]+=1),p->(p.fragmentation.shape[2]*=2),
        p->(p.fragmentation.unpolarized[2]+=1),p->(p.fragmentation.model.nodes[2]+=.01),
        p->(p.fragmentation.radial.b[2]+=.01),p->(p.fragmentation.operator.alpha[2]*=2),
        p->(p.fragmentation.operator.blocks[1].indices[2]=1),p->(p.fragmentation.samples.bases[1,2]+=1),
        p->(p.provenance["scope"]="production"),p->(p.fragmentation.provenance["production_ready"]=true)]
    for mutate in mutations
        p=synthetic_star_point();mutate(p)
        @test_throws ErrorException star_structure(p,t)
    end
    p=synthetic_star_point();f=p.fragmentation
    @test_throws ErrorException star_bin([p],[-1.];phase=:sin_phiS_minus_phiH,measure_source="bad weight")
    @test_throws ErrorException star_bin([p],[1.];phase=:unknown,measure_source="bad phase")
    @test_throws ErrorException star_bin([p],[1.];phase=:sin_phiS_minus_phiH,measure_source=" ")
    @test_throws ErrorException star_bin([p,synthetic_star_point(;fragmentation=:quark_and_gluon)],[1.,1.];phase=:sin_phiS_minus_phiH,measure_source="mixed models")
    for (axis,scale,R) in ((:winner_take_all,:historical_pT_C11,.6),(:standard,:pTR,.6),(:standard,:historical_pT_C11,0.))
        @test_throws ErrorException star_point_plan(p.born,p.weights,p.couplings,p.fA,p.fB,p.operator,p.samples,f;
            jet_radius=R,axis=axis,scale_prescription=scale,measure_source="unsupported model")
    end
    @test_throws ErrorException star_point_plan(p.born,p.weights,2p.couplings,p.fA,p.fB,p.operator,p.samples,f;
        jet_radius=.6,axis=:standard,scale_prescription=:historical_pT_C11,measure_source="wrong coupling")
    @test_throws ErrorException star_fragmentation_plan(f.model,f.radial,f.Q2,f.z,f.jT,f.pids,f.matched,f.operator,f.samples;
        fragmentation=:quark_only,gluon_model=:none,charge=1,provenance=Dict("production_ready"=>true))
    owned=copy(p.fA);owned[1,1]+=1
    @test owned!=p.fA
    @test_throws ErrorException star_momentum_bin(p,MomentumRule(0.,.5;order=4);phase=:sin_phiS_minus_phiH,measure=:unknown,measure_source="bad measure")
    @test CollinsFit.star_native_checked(()->42)==42
    @test_throws ErrorException CollinsFit.star_native_checked(()->begin
        println(stderr,"cgauss64D113.1TOO HIGH ACCURACY REQUIRED");42
    end)
end
