include("../bin/SIAQTChecks.jl")

@testset "explicit SIA qT windows preserve finite P integration and source definitions" begin
    settings=JSON.parsefile(QT_CONTRACT)["synthetic_fixture"]
    context=sia_qt_fixture(.25,settings)
    theta=Float64.(settings["parameters"])
    angular=SIAMeanOfRatio(.6;source="synthetic event sample",averaging_measure="fixed synthetic test measure")
    window=SIAQTWindow(0.,3.5;source="explicit synthetic upper cut",lower_limit_source="explicit synthetic radial origin")
    @test !ismutable(window)
    @test window.lower_GeV==0. && window.upper_GeV==3.5
    for ratio in (:UL,:UC)
        observable=sia_qt_observable(context,window;ratio=ratio,angular=angular,order=32)
        previous=sia_observable(context;ratio=ratio,angular=angular,momentum=MomentumRule(0.,.25*3.5;order=32),source="same synthetic P bounds")
        @test predict_observable(observable,theta)==predict_observable(previous,theta)
        @test ForwardDiff.gradient(t->predict_observable(observable,t),theta)==ForwardDiff.gradient(t->predict_observable(previous,t),theta)
        @test observable.provenance["input_momentum_bounds_GeV"]==[0.,3.5]
        @test observable.provenance["momentum_bounds_GeV"]==[0.,.875]
        @test observable.provenance["input_momentum_measure"]=="z1^2*qT*dqT = P_hperp*dP_hperp"
        @test observable.provenance["fixed_theory_z1"]==.25
        @test observable.provenance["angular_definition"]=="mean_of_angular_ratio"
        @test observable.provenance["angular_averaging_measure"]==angular.averaging_measure
        @test observable.provenance["qT_lower_limit_source"]==window.lower_limit_source
        @test all(observable.provenance[k]===false for k in
            ("production_ready","physical_extraction","G2PERF_pass","historical_fit_ready","experimental_data_admitted","joint_acceptance_recovered"))
        for field in ("input_momentum_bounds_GeV","qT_lower_limit_source","production_ready")
            altered=deepcopy(observable)
            altered.provenance[field]=field=="production_ready" ? true : "changed after preparation"
            @test_throws ErrorException predict_observable(altered,theta)
        end
        altered=deepcopy(observable)
        altered.momentum.weights[1]=nextfloat(altered.momentum.weights[1])
        @test_throws ErrorException predict_observable(altered,theta)
    end
    for bounds in ((-.1,1.),(0.,0.),(1.,.5),(0.,Inf),(NaN,1.),(false,1.))
        @test_throws ErrorException SIAQTWindow(bounds...;source="test",lower_limit_source="explicit test")
    end
    @test_throws UndefKeywordError SIAQTWindow(0.,1.;source="missing lower provenance")
    @test_throws UndefKeywordError SIAQTWindow(0.,1.;lower_limit_source="missing source")
    @test_throws ErrorException SIAQTWindow(0.,1.;source=" ",lower_limit_source="explicit")
    @test_throws ErrorException SIAQTWindow(0.,1.;source="explicit",lower_limit_source=" ")
    @test_throws ErrorException SIAQTWindow(0.,0.,"raw invalid constructor","explicit")
    @test_throws ErrorException ForwardDiff.derivative(t->SIAQTWindow(0.,1. + t;source="fixed",lower_limit_source="fixed").upper_GeV,0.)
    @test_throws ErrorException ForwardDiff.derivative(t->SIAQTWindow(.2+t,1.;source="fixed",lower_limit_source="fixed").lower_GeV,0.)
    for order in (true,8.0,0,-1)
        @test_throws ErrorException sia_qt_observable(context,window;ratio=:UL,angular=angular,order=order)
    end
    @test_throws ErrorException sia_qt_observable(context,window;ratio=:A12,angular=angular)
    @test_throws ErrorException sia_qt_observable(observable_fixture(:SIDIS),window;ratio=:UL,angular=angular)
    underflow=SIAQTWindow(nextfloat(0.),1.;source="synthetic underflow",lower_limit_source="explicit positive subnormal")
    @test_throws ErrorException sia_qt_observable(context,underflow;ratio=:UL,angular=angular)
    @test_throws ErrorException ObservableReference.reference_weights(observable_fixture(:SIDIS),(0.,1.);coordinate=:qT)
    @test_throws ErrorException ObservableReference.reference_weights(context,nothing;coordinate=:qT)
    @test_throws ErrorException ObservableReference.reference_weights(context,(0.,1.);coordinate=:unknown)
end

@testset "fourth-order derivative reference against closed smooth functions" begin
    contract=JSON.parsefile(QT_FD_CONTRACT);policy=contract["policy"]
    setprecision(policy["reference_bits"]) do
        x=BigFloat(.01)
        for (f,closed) in ((t->t^3,3x^2),(t->exp(7t),7exp(7x)))
            steps=BigFloat.(policy["central_steps"])
            raw=[(f(x+h)-f(x-h))/(2h) for h in steps]
            fourth=[qt_richardson(raw[i],raw[i+1]) for i in 1:2]
            @test all(abs(value-closed)/max(abs(closed),policy["relative_floor"])<=policy["derivative_max_scaled"] for value in fourth)
            @test abs(fourth[1]-fourth[2])/max(abs(fourth[2]),policy["relative_floor"])<=policy["derivative_max_scaled"]
        end
        h=BigFloat(first(policy["central_steps"]))
        raw=((x+h)^3-(x-h)^3)/(2h)
        @test abs(raw-3x^2)/abs(3x^2)>policy["derivative_max_scaled"]
        @test Float64(raw-3x^2)==Float64(h^2)
    end
end
