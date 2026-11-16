# Explicitly synthetic entry-point controls. Typed test-only tokens cannot be
# supplied by command-line String arguments and approve no numerical result.
module HistoricalAuditCLIFixture
include("../bin/audit_historical_points.jl")
struct SyntheticPassToken end
struct SyntheticFailureToken end
HistoricalPointAudit.audit(::SyntheticPassToken,::SyntheticPassToken)=Dict{String,Any}(
    "status"=>"PASS","scope"=>"explicitly_synthetic_CLI_result_roundtrip",
    "accepted_prediction_diagnostic"=>false,"accepted_fitting_gate"=>false,
    "production_ready"=>false,"physical_extraction"=>false,"fit_performed"=>false)
HistoricalPointAudit.audit(::SyntheticFailureToken,::SyntheticFailureToken)=
    error("explicitly synthetic CLI exception control")
end

@testset "historical audit CLI preserves library success and error results" begin
    fixture=HistoricalAuditCLIFixture
    mktempdir() do directory
        output=joinpath(directory,"synthetic-pass.json")
        pass=fixture.SyntheticPassToken()
        result=fixture.main([pass,pass,output])
        @test result["status"]=="PASS"
        @test result["scope"]=="explicitly_synthetic_CLI_result_roundtrip"
        @test result["accepted_prediction_diagnostic"]===false && result["accepted_fitting_gate"]===false
        @test JSON.parsefile(output)==result
        before=read(output)
        @test_throws ErrorException fixture.main([pass,pass,output])
        @test read(output)==before
        failed=fixture.SyntheticFailureToken();errorfile=joinpath(directory,"synthetic-error.json")
        failure=fixture.main([failed,failed,errorfile])
        @test failure["status"]=="FAIL"
        @test occursin("explicitly synthetic CLI exception control",failure["error"])
        @test failure["accepted_prediction_diagnostic"]===false && failure["physical_extraction"]===false
        @test JSON.parsefile(errorfile)==failure
    end
end
