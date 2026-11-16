@testset "prospective SIDIS resource policy and complete source binding" begin
    policy=historical_sidis_pilot_policy()
    @test all(policy[k]===false for k in CollinsFit.HISTORICAL_SIDIS_PILOT_FALSE_FLAGS)
    @test policy["selection"]["partial_coverage_only"]===true
    @test policy["numerics"]["derivative_evaluations"]["derivative_accuracy_accepted"]===false
    # Independent Python hashlib/json golden for sorted ASCII path/hash bindings.
    fixture=Dict("sub/z.json"=>"f"^64,"a.py"=>"0"^64)
    @test sidis_pilot_source_fingerprint(fixture)=="6a19a9b651c89389e19af6f2e67cda870133e2fccd5b61fef5093c4fa803dc5b"
    root=normpath(joinpath(@__DIR__,"../.."));sources=sidis_pilot_current_sources(root)
    @test sources["config/historical-sidis-point-contract.json"]=="831d4f3435a1eea61b0ca664b7b61b1fa56b19f930d0e516d7e23bcb7c61ba45"
    @test all(haskey(sources,p) for p in ("src/collinsflow/historical_sidis_pilot.py",
        "numerics/bin/pilot_historical_sidis_point.jl","numerics/bin/audit_historical_sidis_pilot.jl",
        "workflow/historical-sidis-pilot.smk","config/historical-sidis-pilot-policy.json"))
    @test all(!startswith(p,"runs/") && !startswith(p,"data/raw/") && !occursin("__pycache__",p) for p in keys(sources))
end

@testset "synthetic pilot scalar arithmetic cannot certify a wrong source factor" begin
    policy=historical_sidis_pilot_policy()
    for (publication,factor) in (("1006.4221v1",.8),("1408.4405",-1.),("0802.2160",-1.))
        point=Dict("publication"=>publication,"kinematics"=>Dict("y"=>.5))
        @test sidis_pilot_source_factor(point)==factor
        correct=Dict("UU"=>2.,"FUT"=>.25,"amplitude"=>factor/8)
        checks=sidis_pilot_value_checks(correct,copy(correct),point,policy)
        @test length(checks)==4 && all(c["passed"] for c in checks)
        wrong=copy(correct);wrong["amplitude"]=-wrong["amplitude"]
        @test !all(c["passed"] for c in sidis_pilot_value_checks(wrong,correct,point,policy))
        negative=copy(correct);negative["UU"]=-2.
        @test_throws ErrorException sidis_pilot_value_checks(negative,negative,point,policy)
    end
end

@testset "explicitly synthetic closed13 evaluation gradient and HVP timing samples" begin
    theta=collect(1:13).*0.01;direction=collect(1:13).*0.001;weights=Float64.(1:13)
    value(t)=sum(weights.*t.^3)+sum(t)^2/2
    samples=sidis_pilot_evaluations(value,theta,direction,3)
    @test samples.amplitude==value(theta)
    expected_gradient=3 .* weights .* theta.^2 .+ sum(theta)
    expected_hvp=6 .* weights .* theta .* direction .+ sum(direction)
    @test isapprox(samples.samples["gradient"],expected_gradient;rtol=1e-13)
    @test isapprox(samples.samples["HVP"],expected_hvp;rtol=1e-13)
    @test samples.finite && samples.repeatable
    @test samples.samples["derivative_accuracy_accepted"]===false
    @test samples.samples["independent_derivative_reference_executed"]===false
    @test length(samples.times["warm_repetitions"])==3
    @test all(isfinite(t) && t>=0 for r in samples.times["warm_repetitions"] for t in values(r["seconds"]))
    @test_throws ErrorException sidis_pilot_evaluations(value,theta[1:12],direction,3)
    @test_throws ErrorException sidis_pilot_evaluations(value,theta,direction,true)
end

# Loading either executable must define its jobs without starting a provider,
# source replay, native preparation or prediction during package/unit loading.
module SIDISPilotProducerParseProbe
include(joinpath(@__DIR__,"../bin/pilot_historical_sidis_point.jl"))
end
module SIDISPilotAuditParseProbe
include(joinpath(@__DIR__,"../bin/audit_historical_sidis_pilot.jl"))
end
@testset "pilot executables bind the current trusted source root" begin
    root=normpath(joinpath(@__DIR__,"../.."))
    @test SIDISPilotProducerParseProbe.PILOT_ROOT==root
    @test SIDISPilotAuditParseProbe.PILOT_AUDIT_ROOT==root
end
