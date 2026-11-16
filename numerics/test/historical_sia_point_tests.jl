using SHA

module HistoricalSIASourceFixture
include("../bin/predict_historical_sia_points.jl")
end

@testset "synthetic RF0 linearization projects both ratios and rejects angular counterfactuals" begin
    policy=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-point-contract.json"))["numerical_policy"]
    points=[Dict("id"=>"explicit_synthetic_$(ratio)_$(f)","ratio"=>ratio,"source_angular"=>Dict("value"=>f))
        for f in (.5,.6,.75) for ratio in ("UL","UC")]
    rows=HistoricalSIASourceFixture.source_projection(points,policy)
    @test length(rows)==12
    for row in rows
        expected=row["angular_factor"]*row["sign"]/(row["ratio"]=="UL" ? 8 : 24)
        @test row["projected"] ≈ expected rtol=1e-12
        @test row["implementation"] ≈ expected rtol=1e-12
        @test row["scope"]=="explicitly_synthetic_source_definition_fixture"
        @test row["positive_counts_proved_by_bound"] && row["all_counterfactuals_differ_from_correct"]
        @test abs(row["counterfactuals"]["cos_phi0"])<=1e-20
        @test abs(row["counterfactuals"]["sin_2phi0"])<=1e-20
        @test row["error"]<=policy["source_projection_max_scaled"]
    end
end

# This direct typed fixture tests model composition only. It never passes the
# source-loader factory, replays raw data or invokes the real numerical job.
function synthetic_sia_model_request(directory,name,base;measured=.0125)
    points=Dict{String,Any}[]
    for ratio in ("UL","UC"),i in 1:9
        push!(points,Dict("id"=>"explicit_synthetic_$(ratio)_$(i)","ratio"=>ratio,"publication"=>"synthetic",
            "kinematics"=>Dict("z1"=>base.leptonic.z1,"z2"=>base.leptonic.z2,"Q2_GeV2"=>base.leptonic.Q2,"pT_GeV"=>.05+.1(i-1)),
            "measured"=>Dict("central_value"=>measured),
            "source_angular"=>Dict("value"=>.6,"averaging_measure"=>"explicitly synthetic fixed mathematical measure")))
    end
    source=joinpath(directory,"synthetic-source.txt")
    isfile(source) || write(source,"Explicitly synthetic model fixture; no experimental evidence.\n")
    data=Dict("scope"=>"historical_reconstructed_sia_point_prediction_diagnostic","fixture_scope"=>"explicitly_synthetic_model_only",
        "production_ready"=>false,"physical_extraction"=>false,"fit_performed"=>false,"points"=>points,
        "source_inputs"=>Dict(source=>bytes2hex(sha256(read(source)))))
    file=joinpath(directory,name);writejson(file,data);owned=JSON.parsefile(file)
    HistoricalSIAPointRequest(owned,file,bytes2hex(sha256(read(file))),CollinsFit.numerical_snapshot(owned))
end

@testset "SIA historical composition shares rows and keeps observations out of prediction" begin
    mktempdir() do directory
        base=observable_fixture(:SIA)
        request=synthetic_sia_model_request(directory,"first.json",base)
        other=synthetic_sia_model_request(directory,"changed-observations.json",base;measured=-100.)
        contexts=Dict(historical_sia_geometry_key(p)=>with_momentum(base,p["kinematics"]["pT_GeV"]) for p in request.data["points"])
        @test length(contexts)==9
        model=historical_sia_point_model(request,contexts);theta=collect(kpsy_published_parameters().values)
        actual=predict_joint(model,theta)
        @test length(actual)==18 && all(isfinite,actual)
        @test actual==predict_joint(historical_sia_point_model(other,contexts),theta)
        for (point,value,observation) in zip(request.data["points"],actual,model.observations)
            context=contexts[historical_sia_geometry_key(point)]
            pp,pm,mp,mm=[polarized_structure(with_pion_charges(context,c...),theta) for c in ((1,1),(1,-1),(-1,1),(-1,-1))]
            u=pm.Collins+mp.Collins;l=pp.Collins+mm.Collins;du=pm.UU+mp.UU;dl=pp.UU+mm.UU
            expected=.6*(u/du-(point["ratio"]=="UL" ? l/dl : (u+l)/(du+dl)))
            @test value ≈ expected rtol=1e-12
            @test observation.plan.context.first_operator===context.first_operator
        end
        changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
        @test predict_joint(model,changed)==actual
        flipped=copy(theta);flipped[[7,8]].*=-1
        @test predict_joint(model,flipped)==actual
        @test_throws ErrorException historical_sia_point_model(request,Dict())
        @test_throws ErrorException JointGaussianLikelihood(model,model.point_ids,zeros(18),Matrix{Float64}(I,18,18);
            scope=request.data["scope"],source="explicitly synthetic diagnostic-admission rejection")
        request.data["points"][1]["source_angular"]["value"]*=2
        @test_throws ErrorException historical_sia_request_valid(request)
    end
    mktempdir() do directory
        file=joinpath(directory,"wrong-scope.json")
        for scope in ("production","synthetic_closure","experimental_candidate")
            writejson(file,Dict("schema_version"=>1,"scope"=>scope))
            @test_throws ErrorException load_historical_sia_point_request(file)
        end
    end
end

module HistoricalSIAAuditCLIFixture
include("../bin/audit_historical_sia_points.jl")
struct SyntheticPassToken end
struct SyntheticFailureToken end
HistoricalSIAPointAudit.audit(::SyntheticPassToken,::SyntheticPassToken)=Dict{String,Any}(
    "status"=>"PASS","scope"=>"explicitly_synthetic_SIA_CLI_result_roundtrip",
    "accepted_prediction_diagnostic"=>false,"accepted_fitting_gate"=>false,
    "production_ready"=>false,"physical_extraction"=>false,"fit_performed"=>false)
HistoricalSIAPointAudit.audit(::SyntheticFailureToken,::SyntheticFailureToken)=error("explicitly synthetic SIA CLI exception")
end

@testset "SIA audit CLI preserves returned and failed state without admission" begin
    fixture=HistoricalSIAAuditCLIFixture
    mktempdir() do directory
        output=joinpath(directory,"synthetic-pass.json");pass=fixture.SyntheticPassToken()
        result=fixture.main([pass,pass,output])
        @test result["status"]=="PASS" && result["scope"]=="explicitly_synthetic_SIA_CLI_result_roundtrip"
        @test result["accepted_prediction_diagnostic"]===false && result["accepted_fitting_gate"]===false
        @test JSON.parsefile(output)==result
        before=read(output)
        @test_throws ErrorException fixture.main([pass,pass,output])
        @test read(output)==before
        failed=fixture.SyntheticFailureToken();errorfile=joinpath(directory,"synthetic-error.json")
        failure=fixture.main([failed,failed,errorfile])
        @test failure["status"]=="FAIL" && occursin("explicitly synthetic SIA CLI exception",failure["error"])
        @test failure["accepted_prediction_diagnostic"]===false && failure["physical_extraction"]===false
        @test JSON.parsefile(errorfile)==failure
    end
end
