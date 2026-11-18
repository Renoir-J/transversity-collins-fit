using Test, JSON

# This namespace deliberately does not import Dates, JSON or now. The actual
# report runner must finalize its own successful and failed reports in isolation.
module IsolatedRetainedEntryPoint
include("../bin/SIDISRetainedPairReview.jl")
module SyntheticBindings
using JSON
require_new_output(m,path)=!ispath(path) || error("retain occupied synthetic output")
provenance(m)=Dict("scope"=>"explicitly_synthetic_reader_harness_not_a_source_manifest")
origin(run,phase,index)=run
stable(m)=true
write_record(path,record)=write(path,JSON.json(record))
end
module SyntheticSourceAPI
using SHA
false_flags()=Dict(key=>false for key in ("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted",
    "fit_performed","accepted_prediction_diagnostic","accepted_fitting_gate","G2PERF_pass","interval_coverage_verified","all132_coverage_verified"))
bindings(request,native,parton)=(fp="EXPLICITLY_SYNTHETIC_NOT_NUMERICAL_SOURCE",request=(sha256="SYNTHETIC_REQUEST_ONLY",),fail=request=="SYNTHETIC_FAIL")
point_folder(run,phase,index)=joinpath(run,phase,"p"*lpad(string(index),4,'0'))
function audit_point(ctx,producer,phase,index;complete_invocation=false)
    complete_invocation===true || error("synthetic harness requires complete-pair reader call")
    ctx.fail && index==2 && error("synthetic saved-record rejection")
    Dict("scope"=>"explicitly_synthetic_non_numerical_packet")
end
stable(ctx)=true
digest(path)=bytes2hex(sha256(read(path)))
end
transport_check(run,phase,index,packet;complete_invocation=false)=complete_invocation===true || error("missing complete transport check")
end

@testset "Isolated retained-reader entrypoint writes success, preserves failure and refuses overwrite" begin
    fixture=IsolatedRetainedEntryPoint
    @test !isdefined(fixture,:now) && !isdefined(fixture,:Dates) && !isdefined(fixture,:JSON)
    mktempdir() do temp
        rows=Any[]
        for index in 1:3
            folder=fixture.SyntheticSourceAPI.point_folder(temp,"native123",index);mkpath(folder)
            write(joinpath(folder,"invocation.json"),"EXPLICITLY_SYNTHETIC_NOT_A_COMPLETE_NUMERICAL_JOB")
            push!(rows,Dict("phase"=>"native123","index"=>index,"folder"=>folder))
        end
        data=Dict{String,Any}("request"=>"SYNTHETIC_SUCCESS","native"=>"UNUSED_SYNTHETIC_NATIVE","parton"=>"UNUSED_SYNTHETIC_PARTON",
            "retained_jobs"=>rows,"theory_scope"=>"synthetic_closure")
        manifest=(data=data,run=temp)
        output=joinpath(temp,"success.json")
        result=fixture.SIDISRetainedPairReview.run(manifest,output,fixture.SyntheticSourceAPI,fixture.SyntheticBindings,fixture.transport_check)
        @test result==0
        report=JSON.parsefile(output)
        @test report["status"]=="RETAINED_POINT_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE"
        @test report["theory_scope"]=="synthetic_closure" && report["source_fingerprint"]=="EXPLICITLY_SYNTHETIC_NOT_NUMERICAL_SOURCE"
        @test length(report["jobs_reviewed"])==length(report["retained_invocations"])==3
        @test all(report[key]===false for key in keys(fixture.SyntheticSourceAPI.false_flags()))
        @test !isempty(report["completed_at"])
        before=read(output)
        @test_throws ErrorException fixture.SIDISRetainedPairReview.run(manifest,output,fixture.SyntheticSourceAPI,fixture.SyntheticBindings,fixture.transport_check)
        @test read(output)==before
        data["request"]="SYNTHETIC_FAIL";failure=joinpath(temp,"failure.json")
        result=fixture.SIDISRetainedPairReview.run(manifest,failure,fixture.SyntheticSourceAPI,fixture.SyntheticBindings,fixture.transport_check)
        report=JSON.parsefile(failure)
        @test result==2 && report["status"]=="FAIL"
        @test occursin("synthetic saved-record rejection",report["error"])
        @test !isempty(report["completed_at"]) && length(report["jobs_reviewed"])==1
        @test all(report[key]===false for key in keys(fixture.SyntheticSourceAPI.false_flags()))
    end
end
