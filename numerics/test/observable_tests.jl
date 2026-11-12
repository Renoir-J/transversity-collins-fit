using QuadGK

function observable_fixture(process)
    model=KPSYCSS();radial=RadialRule([0.,.25,1.,4.,8.];order=8)
    plan=leptonic_plan(process,model,radial,process==:SIDIS ? 4. : 110.,.2,.3;z2=process==:SIA ? .5 : nothing)
    first=[(1+.07i)*(1+.03b) for i in 1:10,b in radial.b]
    second=[(1+.11i)*(1+.05b) for i in 1:10,b in radial.b]
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    op=synthetic_polarized_operator(process,kind,process==:SIDIS ? [.15] : [.3,.5],radial,model)
    n=process==:SIDIS ? 2 : 4
    samples=PolarizedInitialSamples(kind,op.grid,[(1+.3i)*x*(1-x)^3 for i in 1:n,x in op.grid];input_kind=:analytic)
    op2=process==:SIDIS ? synthetic_polarized_operator(process,:KPSY_Hhat3_homogeneous,[.3],radial,model) : nothing
    samples2=op2===nothing ? nothing : PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,
        [(1+.3i)*x*(1-x)^3 for i in 1:4,x in op2.grid];input_kind=:analytic)
    polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=samples2,
        target=process==:SIDIS ? NucleonTarget(:proton) : nothing,charges=process==:SIDIS ? (1,0) : (1,-1))
end

@testset "explicit SIA angular source types and averaging measures" begin
    # Explicitly synthetic angles and weights, not experimental averages.
    cos2=[0.,1.];sin2=1 .-cos2;uniform=[.5,.5]
    mean_of_ratio=sum(uniform.*sin2./(1 .+cos2))
    ratio_of_means=sum(uniform.*sin2)/(1+sum(uniform.*cos2))
    event_weights=uniform.*(1 .+cos2);event_weights./=sum(event_weights)
    @test mean_of_ratio==.5
    @test ratio_of_means ≈ 1/3 rtol=1e-14
    @test mean_of_ratio != ratio_of_means
    @test sum(event_weights.*sin2./(1 .+cos2)) ≈ ratio_of_means rtol=1e-14
    angular=SIAMeanOfRatio(.6;source="synthetic event average",averaging_measure="fixed synthetic event weights")
    @test angular.value==.6 && angular.averaging_measure=="fixed synthetic event weights"
    for value in (-.01,1.01,NaN,Inf,true)
        @test_throws ErrorException SIAMeanOfRatio(value;source="synthetic",averaging_measure="synthetic weights")
    end
    @test_throws ErrorException SIAMeanOfRatio(.6;source=" ",averaging_measure="synthetic weights")
    @test_throws ErrorException SIAMeanOfRatio(.6;source="synthetic",averaging_measure=" ")
    @test_throws UndefKeywordError SIAMeanOfRatio(.6;source="missing measure")
    @test_throws ErrorException ForwardDiff.derivative(t->SIAMeanOfRatio(.6+t;source="synthetic",averaging_measure="fixed").value,0.)
    sums=sia_charge_sums((UU=100.,Collins=10.),(UU=10.,Collins=2.),(UU=20.,Collins=6.),(UU=1.,Collins=-.1))
    ratio_input=SIAAngularFactor(.6;source="synthetic acceptance average")
    for ratio in (:UL,:UC)
        @test sia_double_ratio(sums,ratio,angular)==sia_double_ratio(sums,ratio,ratio_input)
    end
    context=observable_fixture(:SIA);theta=collect(kpsy_published_parameters().values)
    for interval in (nothing,(.05,.35)),ratio in (:UL,:UC)
        momentum=interval===nothing ? nothing : MomentumRule(interval...;order=8)
        point=sia_observable(context;ratio=ratio,angular=angular,momentum=momentum,source="synthetic source-kind fixture")
        previous=sia_observable(context;ratio=ratio,angular=ratio_input,momentum=momentum,source="synthetic source-kind fixture")
        @test predict_observable(point,theta)==predict_observable(previous,theta)
        @test point.provenance["angular_definition"]=="mean_of_angular_ratio"
        @test point.provenance["angular_averaging_measure"]==angular.averaging_measure
        @test previous.provenance["angular_definition"]=="ratio_of_angular_means"
        @test point.provenance["production_ready"]===false && point.provenance["physical_extraction"]===false
        @test ForwardDiff.gradient(t->predict_observable(point,t),theta)==ForwardDiff.gradient(t->predict_observable(previous,t),theta)
        @test ForwardDiff.hessian(t->predict_observable(point,t),theta)==ForwardDiff.hessian(t->predict_observable(previous,t),theta)
        point.provenance["angular_averaging_measure"]="changed after preparation"
        @test_throws ErrorException predict_observable(point,theta)
    end
end

@testset "momentum measure, depolarization and charge-category algebra" begin
    for order in (8,16,32)
        rule=MomentumRule(.1,.9;order=order)
        for k in 0:8
            @test sum(rule.weights.*rule.nodes.^k) ≈ (.9^(k+2)-.1^(k+2))/(k+2) rtol=1e-13 atol=1e-15
        end
    end
    rule=MomentumRule(0.,1.)
    correct=sum(rule.weights.*rule.nodes.^2)/sum(rule.weights.*(1 .+rule.nodes.^2))
    wrong=sum(rule.weights.*rule.nodes.^2 ./ (1 .+rule.nodes.^2))/sum(rule.weights)
    @test correct ≈ 1/3 rtol=1e-13
    @test wrong ≈ 1-log(2) rtol=1e-13
    @test abs(correct-wrong)>.02
    @test depolarization_factor(0.)==1.
    @test depolarization_factor(1.)==0.
    @test depolarization_factor(.5)==.8
    @test_throws ErrorException depolarization_factor(-.1)
    @test_throws ErrorException ForwardDiff.derivative(depolarization_factor,.5)
    @test_throws ErrorException MomentumRule(.2,.2)
    @test_throws ErrorException MomentumRule(-.1,1.)
    @test_throws ErrorException MomentumRule(0.,Inf)
    @test_throws ErrorException MomentumRule(0.,1.;order=true)
    @test_throws ErrorException ForwardDiff.derivative(t->sum(MomentumRule(0.,1.0+t).weights),0.)
    angular=SIAAngularFactor(.5;source="synthetic angular ratio")
    sums=sia_charge_sums((UU=100.,Collins=10.),(UU=10.,Collins=2.),(UU=20.,Collins=6.),(UU=1.,Collins=-.1))
    @test sums.U==(UU=30.,Collins=8.)
    @test sums.L==(UU=101.,Collins=9.9)
    @test sums.C==(UU=131.,Collins=17.9)
    @test sia_double_ratio(sums,:UL,angular) ≈ .5*(8/30-9.9/101) rtol=1e-13
    @test sia_double_ratio(sums,:UC,angular) ≈ .5*(8/30-17.9/131) rtol=1e-13
    @test sia_double_ratio(sums,:UC,angular) ≈ (101/131)*sia_double_ratio(sums,:UL,angular) rtol=1e-13
    @test_throws ErrorException sia_double_ratio(sums,:bad,angular)
    @test_throws ErrorException sia_double_ratio((U=sums.U,L=sums.L,C=(UU=1.,Collins=1.)),:UC,angular)
    @test_throws ErrorException sia_double_ratio((U=(UU=0.,Collins=1.),L=sums.L,C=(UU=101.,Collins=10.9)),:UL,angular)
    @test_throws ErrorException SIAAngularFactor(1.01;source="test")
    @test_throws ErrorException SIAAngularFactor(.5;source="")
    @test_throws ErrorException SIAAngularFactor(true;source="test")
    @test_throws ErrorException ForwardDiff.derivative(t->SIAAngularFactor(.5+t;source="test").value,0.)
end

@testset "complete observable point values and independently integrated momentum bins" begin
    theta=collect(kpsy_published_parameters().values)
    angular=SIAAngularFactor(.6;source="explicit ratio of angular means")
    for process in (:SIDIS,:SIA)
        context=observable_fixture(process)
        point=process==:SIDIS ? sidis_observable(context;definition=:structure_function_ratio,source="synthetic point") :
            sia_observable(context;ratio=:UL,angular=angular,source="synthetic point")
        if process==:SIDIS
            structure=polarized_structure(context,theta)
            @test predict_observable(point,theta) ≈ structure.Collins/structure.UU rtol=1e-12
            raw=sidis_observable(context;definition=:sin_phi_sum_amplitude,y=.4,source="synthetic point")
            @test predict_observable(raw,theta) ≈ depolarization_factor(.4)*predict_observable(point,theta) rtol=1e-12
            @test_throws ErrorException sidis_observable(context;definition=:sin_phi_sum_amplitude,source="missing y")
            @test_throws ErrorException sidis_observable(context;definition=:sin_phi_sum_amplitude,y=0.,source="invalid y")
        else
            pairs=[polarized_structure(with_pion_charges(context,c1,c2),theta) for (c1,c2) in ((1,1),(1,-1),(-1,1),(-1,-1))]
            sums=sia_charge_sums(pairs...)
            @test predict_observable(point,theta) ≈ sia_double_ratio(sums,:UL,angular) rtol=1e-12
            uc=sia_observable(context;ratio=:UC,angular=angular,source="synthetic point")
            @test predict_observable(uc,theta) ≈ sia_double_ratio(sums,:UC,angular) rtol=1e-12
        end
        low,high=.05,.35
        function adaptive_pair(c1,c2)
            charged=process==:SIDIS ? context : with_pion_charges(context,c1,c2)
            quadgk(low,high;rtol=1e-12,atol=1e-14) do P
                r=polarized_structure(with_momentum(charged,P),theta)
                P.*[r.UU,r.Collins]
            end[1]
        end
        reference=if process==:SIDIS
            v=adaptive_pair(context.charges...);v[2]/v[1]
        else
            pairs=[adaptive_pair(c1,c2) for (c1,c2) in ((1,1),(1,-1),(-1,1),(-1,-1))]
            sums=sia_charge_sums([(UU=p[1],Collins=p[2]) for p in pairs]...)
            sia_double_ratio(sums,:UL,angular)
        end
        for order in (8,16,32)
            rule=MomentumRule(low,high;order=order)
            observable=process==:SIDIS ? sidis_observable(context;definition=:structure_function_ratio,momentum=rule,source="synthetic P bin") :
                sia_observable(context;ratio=:UL,angular=angular,momentum=rule,source="synthetic P bin")
            @test predict_observable(observable,theta) ≈ reference rtol=1e-10 atol=1e-14
            @test observable.provenance["production_ready"]===false && observable.provenance["physical_extraction"]===false
            @test observable.provenance["momentum_measure"]=="P*dP"
            if order==8
                saved=observable.collins_weights[1];observable.collins_weights[1]=nextfloat(saved)
                @test_throws ErrorException predict_observable(observable,theta)
                observable.collins_weights[1]=saved
                rule.nodes[1]=nextfloat(rule.nodes[1])
                @test_throws ErrorException predict_observable(observable,theta)
            end
        end
        @test_throws ErrorException process==:SIDIS ? sia_observable(context;ratio=:UL,angular=angular,source="bad process") :
            sidis_observable(context;definition=:structure_function_ratio,source="bad process")
    end
end

@testset "observable 13-parameter derivatives through projected momentum weights" begin
    theta=collect(kpsy_published_parameters().values)
    angular=SIAAngularFactor(.6;source="synthetic angular ratio")
    rule=MomentumRule(.05,.35;order=8)
    for process in (:SIDIS,:SIA)
        context=observable_fixture(process)
        observable=process==:SIDIS ? sidis_observable(context;definition=:sin_phi_sum_amplitude,y=.4,momentum=rule,source="synthetic derivative bin") :
            sia_observable(context;ratio=:UL,angular=angular,momentum=rule,source="synthetic derivative bin")
        function manual(t)
            function integrate(c1,c2)
                charged=process==:SIDIS ? context : with_pion_charges(context,c1,c2)
                values=[polarized_structure(with_momentum(charged,P),t) for P in rule.nodes]
                (UU=sum(w*r.UU for (w,r) in zip(rule.weights,values)),Collins=sum(w*r.Collins for (w,r) in zip(rule.weights,values)))
            end
            if process==:SIDIS
                value=integrate(context.charges...)
                depolarization_factor(.4)*value.Collins/value.UU
            else
                pairs=[integrate(c1,c2) for (c1,c2) in ((1,1),(1,-1),(-1,1),(-1,-1))]
                sia_double_ratio(sia_charge_sums(pairs...),:UL,angular)
            end
        end
        prediction(t)=predict_observable(observable,t)
        @test prediction(theta) ≈ manual(theta) rtol=1e-12
        @test ForwardDiff.gradient(prediction,theta) ≈ ForwardDiff.gradient(manual,theta) rtol=1e-11 atol=1e-14
        @test ForwardDiff.hessian(prediction,theta) ≈ ForwardDiff.hessian(manual,theta) rtol=1e-10 atol=1e-13
        direction=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,.001]
        J=ForwardDiff.gradient(prediction,theta)
        hv=ForwardDiff.jacobian(s->ForwardDiff.gradient(prediction,theta+s[1]*direction),[0.])[:,1]
        for step in (1e-4,1e-5)
            fd,fdh=setprecision(256) do
                t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
                (Float64((-3manual(t)+4manual(t+h*d)-manual(t+2h*d))/(2h)),
                 Float64.((-3ForwardDiff.gradient(manual,t)+4ForwardDiff.gradient(manual,t+h*d)-ForwardDiff.gradient(manual,t+2h*d))/(2h)))
            end
            @test abs(fd-dot(J,direction))/max(abs(dot(J,direction)),1e-8) <= 1e-7
            @test maximum(abs.(fdh-hv)./max.(abs.(hv),1e-8)) <= 1e-7
        end
    end
end
