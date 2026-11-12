@testset "explicit SIDIS source phase and normalization" begin
    source="explicitly synthetic SIDIS convention fixture"
    removed=SIDISAsymmetryConvention(phase=:sin_phi_sum_minus_pi,normalization=:DNN_removed,source=source)
    @test removed.y===nothing && removed.source==source
    for phase in (:invalid,"sin_phi_sum",nothing)
        @test_throws ErrorException SIDISAsymmetryConvention(;phase,normalization=:DNN_removed,source=source)
    end
    for normalization in (:invalid,"DNN_removed",nothing)
        @test_throws ErrorException SIDISAsymmetryConvention(;phase=:sin_phi_sum,normalization,source=source)
    end
    for y in (nothing,0.,-.1,1.01,true,NaN,Inf)
        @test_throws ErrorException SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:fixed_y_DNN_included,y=y,source=source)
    end
    @test_throws ErrorException SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:DNN_removed,y=.4,source=source)
    @test_throws ErrorException SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:DNN_removed,source=" ")
    @test_throws ErrorException SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:DNN_removed,source=nothing)
    @test_throws UndefKeywordError SIDISAsymmetryConvention(phase=:sin_phi_sum,source=source)
    @test_throws MethodError SIDISAsymmetryConvention(:sin_phi_sum,:DNN_removed,nothing,source,-1.)
    @test_throws ErrorException ForwardDiff.derivative(t->SIDISAsymmetryConvention(
        phase=:sin_phi_sum,normalization=:fixed_y_DNN_included,y=t,source=source).y,.4)
    context=observable_fixture(:SIDIS);theta=collect(kpsy_published_parameters().values)
    @test_throws ErrorException sidis_observable(observable_fixture(:SIA),removed)
    for interval in (nothing,(.05,.35))
        momentum=interval===nothing ? nothing : MomentumRule(interval...;order=16)
        canonical=sidis_observable(context;definition=:structure_function_ratio,momentum=momentum,source=source)
        raw=sidis_observable(context;definition=:sin_phi_sum_amplitude,y=2/5,momentum=momentum,source=source)
        value=predict_observable(canonical,theta)
        gradient=ForwardDiff.gradient(t->predict_observable(canonical,t),theta)
        hessian=ForwardDiff.hessian(t->predict_observable(canonical,t),theta)
        @test isfinite(value) && value!=0 && any(!iszero,gradient) && any(!iszero,hessian)
        for phase in (:sin_phi_sum,:sin_phi_sum_minus_pi),normalization in (:DNN_removed,:fixed_y_DNN_included)
            y=normalization==:DNN_removed ? nothing : 2/5
            convention=SIDISAsymmetryConvention(;phase,normalization,y,source)
            observable=sidis_observable(context,convention;momentum)
            # Independent rational factors: at y=2/5, DNN=15/17.
            factor=(phase==:sin_phi_sum ? 1. : -1.)*(normalization==:DNN_removed ? 1. : 15/17)
            @test predict_observable(observable,theta) ≈ factor*value rtol=1e-12 atol=1e-16
            @test ForwardDiff.gradient(t->predict_observable(observable,t),theta) ≈ factor.*gradient rtol=1e-12 atol=1e-16
            @test ForwardDiff.hessian(t->predict_observable(observable,t),theta) ≈ factor.*hessian rtol=1e-12 atol=1e-16
            @test observable.provenance["sidis_convention_contract_id"]=="KPSY-SIDIS-CONVENTION-001"
            @test observable.provenance["normalization"]==string(normalization)
            @test observable.provenance["production_ready"]===false && observable.provenance["physical_extraction"]===false
            if phase==:sin_phi_sum
                previous=normalization==:DNN_removed ? canonical : raw
                @test predict_observable(observable,theta)==predict_observable(previous,theta)
            end
            for key in ("phase_definition","normalization","source","DNN_evaluation","sidis_convention_contract_id")
                old=observable.provenance[key];observable.provenance[key]="mutated source convention"
                @test_throws ErrorException predict_observable(observable,theta)
                observable.provenance[key]=old
            end
            @test predict_observable(observable,theta) ≈ factor*value rtol=1e-12 atol=1e-16
        end
        shifted=sidis_observable(context,removed;momentum)
        @test predict_observable(shifted,theta)!=value
        @test predict_observable(shifted,theta)!=-(15/17)*value
        endpoint=SIDISAsymmetryConvention(phase=:sin_phi_sum_minus_pi,normalization=:fixed_y_DNN_included,y=1.,source=source)
        @test iszero(predict_observable(sidis_observable(context,endpoint;momentum),theta))
    end
end
