using SHA
module HistoricalSIAIntegralProducerFixture
include("../bin/predict_historical_sia_integrals.jl")
end
module HistoricalSIAIntegralAuditFixture
include("../bin/audit_historical_sia_integrals.jl")
end

@testset "inward fourth-order controls retain boundary-valid samples and all13 columns" begin
    reference=HistoricalSIAIntegralProducerFixture.HIR
    policy=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))["representative_derivative_policy"]
    controls=reference.mathematical_controls(policy)
    @test controls["scope"]=="explicitly_synthetic_inward_formula_control"
    @test controls["passed"] && controls["fourth_order_halving_verified"]
    @test length(controls["records"])==14
    theta=collect(kpsy_published_parameters().values)
    candidate=t->[sum(t[7:13].^2),t[7]*t[8]+exp(t[13])]
    independent=t->[sum(t[j]*t[j] for j in 7:13),exp(t[13])+t[8]*t[7]]
    d=reference.derivative_reference(candidate,independent,theta,policy)
    @test d["first_six_gradient_columns_zero"] && d["first_six_HVP_columns_zero"]
    J=reference.matrix_rows(d["gradient"])
    @test size(J)==(2,13)
    @test J[1,7:13] ≈ 2theta[7:13]
    @test d["independent_gradient_error"]<=policy["derivative_max_scaled"]
    @test d["independent_HVP_error"]<=policy["HVP_max_scaled"]
    @test d["adjacent_gradient_error"]<=policy["derivative_max_scaled"]
    @test d["adjacent_HVP_error"]<=policy["HVP_max_scaled"]
    arithmetic=HistoricalSIAIntegralAuditFixture.HistoricalSIAIntegralAudit
    verified=arithmetic.inward_checks(d,policy)
    @test length(verified)==9
    @test all(error<=limit for (_,error,limit) in verified)
    for mutation in ("sample","gradient","zero_bit")
        broken=deepcopy(d)
        if mutation=="sample"
            broken["finite_differences"][1]["gradient_samples_decimal"][7][1][1]="1.0"
        elseif mutation=="gradient"
            broken["gradient"][1][7]+=1
        else
            broken["first_six_gradient_columns_zero"]=false
        end
        @test_throws ErrorException arithmetic.inward_checks(broken,policy)
    end
    for r in d["finite_differences"]
        @test r["gradient_error"]<=policy["derivative_max_scaled"]
        @test r["HVP_error"]<=policy["HVP_max_scaled"]
        @test length(r["gradient_samples_decimal"])==13
        @test all(length(c)==5 for c in r["gradient_samples_decimal"])
        @test length(r["HVP_jacobian_samples_decimal"])==5
    end
end

# Direct typed fixtures exercise model composition only; they never pass the
# source factory, claim experimental data, or invoke a retained native solver.
function synthetic_sia_integral_model_request(directory,name,base;measured=.0125,Pmean=.3)
    contract=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))
    points=[Dict("id"=>"explicit_synthetic_integral_$(ratio)","ratio"=>ratio,"publication"=>"synthetic",
        "kinematics"=>Dict("z1"=>base.leptonic.z1,"z2"=>base.leptonic.z2,"Q2_GeV2"=>base.leptonic.Q2),
        "reported_pT_mean_GeV"=>Pmean,"reported_pT_role"=>"source_metadata_not_model_momentum_or_bound",
        "measured"=>Dict("central_value"=>measured),
        "source_angular"=>Dict("value"=>.6,"averaging_measure"=>"explicitly synthetic fixed mathematical measure")) for ratio in ("UL","UC")]
    source=joinpath(directory,"synthetic-source.txt")
    isfile(source) || write(source,"Explicitly synthetic integral model fixture; no experimental evidence.\n")
    data=Dict("scope"=>"historical_reconstructed_sia_integral_prediction_diagnostic","fixture_scope"=>"explicitly_synthetic_model_only",
        "production_ready"=>false,"physical_extraction"=>false,"fit_performed"=>false,"points"=>points,
        "momentum_window"=>CollinsFit.historical_sia_integral_window(contract),
        "source_inputs"=>Dict(source=>bytes2hex(sha256(read(source)))))
    path=joinpath(directory,name);writejson(path,data);owned=JSON.parsefile(path)
    HistoricalSIAIntegralRequest(owned,path,bytes2hex(sha256(read(path))),CollinsFit.numerical_snapshot(owned))
end

@testset "integral model keeps mean P out of geometry and integrates absolute yields" begin
    mktempdir() do directory
        base=with_momentum(observable_fixture(:SIA),0.)
        request=synthetic_sia_integral_model_request(directory,"first.json",base)
        other=synthetic_sia_integral_model_request(directory,"changed-observations.json",base;measured=-100.,Pmean=.4)
        key=historical_sia_integral_geometry_key(first(request.data["points"]));contexts=Dict(key=>base)
        @test length(key)==3
        model=historical_sia_integral_model(request,contexts;order=64);theta=collect(kpsy_published_parameters().values)
        actual=predict_joint(model,theta)
        @test length(actual)==2 && all(isfinite,actual)
        @test actual==predict_joint(historical_sia_integral_model(other,contexts;order=64),theta)
        reference=HistoricalSIAIntegralProducerFixture.HSP.ObservableReference.prepare_generic_reference(base;
            momentum_interval=(0.,3.5),momentum_coordinate=:qT,angular=.6)
        @test actual ≈ reference.predict(theta) rtol=1e-10 atol=1e-12
        for observation in model.observations
            obs=observation.plan
            @test obs.momentum!==nothing && obs.momentum.lower==0. && obs.momentum.upper==base.leptonic.z1*3.5
            @test obs.context.first_operator===base.first_operator
            @test obs.provenance["joint_acceptance_recovered"]===false
            pairs=HistoricalSIAIntegralProducerFixture.pair_integrals(obs,theta)
            pol=reference.polarized_integrals(theta)
            for (i,name) in enumerate(("pp","pm","mp","mm"))
                @test pairs[name]["UU"] ≈ reference.denominators[i] rtol=1e-10 atol=1e-12
                @test pairs[name]["Collins"] ≈ pol[i] rtol=1e-10 atol=1e-12
            end
        end
        for order in (8,true)
            @test_throws ErrorException historical_sia_integral_model(request,contexts;order=order)
        end
        @test_throws ErrorException historical_sia_integral_model(request,Dict(key=>with_momentum(base,.1)))
        @test_throws ErrorException JointGaussianLikelihood(model,model.point_ids,zeros(2),Matrix{Float64}(I,2,2);
            scope=request.data["scope"],source="explicitly synthetic integral admission rejection")
        request.data["momentum_window"]["qT_bounds_GeV"][1]=.1
        @test_throws ErrorException historical_sia_integral_request_valid(request)
    end
end

@testset "whole integral CLIs preserve explicit failure and refuse prior output overwrite" begin
    mktempdir() do directory
        producer=HistoricalSIAIntegralProducerFixture
        audit=HistoricalSIAIntegralAuditFixture
        pred=joinpath(directory,"failed-prediction.json");out=joinpath(directory,"failed-audit.json")
        p=producer.main([joinpath(directory,"absent-request"),"unused-native","unused-input",pred])
        @test p["status"]=="FAIL" && p["physical_extraction"]===false
        @test haskey(p,"error") && isfile(pred)
        original=read(pred)
        @test_throws ErrorException producer.main(["missing","unused","unused",pred])
        @test read(pred)==original
        a=audit.main(["missing",pred,out])
        @test a["status"]=="FAIL" && a["accepted_prediction_diagnostic"]===false && a["accepted_fitting_gate"]===false
        @test haskey(a,"error") && isfile(out)
        original=read(out)
        @test_throws ErrorException audit.main(["missing",pred,out])
        @test read(out)==original
    end
end
