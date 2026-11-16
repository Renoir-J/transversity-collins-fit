using Test, CollinsFit, JSON
include(joinpath(@__DIR__,"../bin/HistoricalSIDISJobRecords.jl"))
include(joinpath(@__DIR__,"../bin/HistoricalSIDISSourceMath.jl"))

@testset "Prospective full SIDIS job policy and exact indices" begin
    p=historical_sidis_job_policy()
    @test p["execution"]["source_points"]==132
    @test p["reporting"]["scientific_tolerances_or_parameter_values_changed"]===false
    @test p["resources"]["native4"]["whole_Julia_jobs_time_budget_seconds"]==3600
    @test sidis_job_point_key(1)=="p0001" && sidis_job_point_key(132)=="p0132"
    for bad in (0,133,true,1.0,"1")
        @test_throws ErrorException sidis_job_point_key(bad)
    end
    @test_throws ErrorException sidis_job_five_point([1.,2.,3.,4.],.1)
    @test_throws ErrorException sidis_job_five_point(ones(5),0.)
end

@testset "Explicitly synthetic full132 native decision; no physical loader" begin
    ids=["SYNTHETIC_CLOSED_LAYOUT:"*sidis_job_point_key(i) for i in 1:132]
    values(level)=Dict("UU"=>1.0+level*1e-7,"FUT"=>.2+level*1e-8,"amplitude"=>.02+level*1e-9)
    baseline=[Dict{String,Any}("phase"=>"native123","point"=>Dict("id"=>id),
        "levels"=>[Dict("level"=>level,"bound_inverse_GeV"=>24,"values"=>values(level)) for level in 1:3]) for id in ids]
    result=sidis_job_native_decision(ids,baseline)
    @test result["selected_level"]==3 && result["status"]=="COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE"
    @test result["production_ready"]===result["all132_coverage_verified"]===result["accepted_prediction_diagnostic"]===false
    bad=deepcopy(baseline);bad[end]["levels"][3]["values"]["amplitude"]+=.001
    coarse=sidis_job_native_decision(ids,bad)
    @test coarse["status"]=="NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" && coarse["selected_level"]==0
    fourth=[Dict{String,Any}("phase"=>"native4","point"=>Dict("id"=>id),
        "levels"=>[Dict("level"=>4,"bound_inverse_GeV"=>24,"values"=>deepcopy(row["levels"][3]["values"]))]) for (id,row) in zip(ids,bad)]
    @test sidis_job_native_decision(ids,bad;fourth)["selected_level"]==4
    fourth[end]["levels"][1]["values"]["amplitude"]+=.002
    @test sidis_job_native_decision(ids,bad;fourth)["status"]=="FAIL"
    @test_throws ErrorException sidis_job_native_decision(ids,baseline;fourth)
    @test_throws ErrorException sidis_job_native_decision(ids[1:end-1],baseline[1:end-1])
    @test_throws ErrorException sidis_job_native_decision(ids,reverse(baseline))
    duplicate=copy(ids);duplicate[end]=ids[1]
    @test_throws ErrorException sidis_job_native_decision(duplicate,baseline)
    negative=deepcopy(baseline);negative[1]["levels"][1]["values"]["UU"]=-1.
    @test_throws ErrorException sidis_job_native_decision(ids,negative)
    wronglevel=deepcopy(baseline);wronglevel[1]["levels"][1]["level"]=true
    @test_throws ErrorException sidis_job_native_decision(ids,wronglevel)
end

@testset "Independent closed polynomial stencil oracle at256 bits" begin
    setprecision(256) do
        x=big"0.01"
        for h in (big"0.0001",big"0.00001",big"0.000001",big"0.0000001"),degree in 0:5
            derivative=degree==0 ? zero(x) : degree*x^(degree-1)
            expected=derivative-(degree==5 ? 24h^4 : zero(x))
            actual=sidis_job_five_point([(x+k*h)^degree for k in 0:4],h)
            @test abs(actual-expected)<big"1e-60"
        end
    end
end

include(joinpath(@__DIR__,"historical_sidis_job_fixtures.jl"))

@testset "All13 Julia directional HVP matches independent coupled polynomial" begin
    t=collect(.1:.1:1.3);d=collect(.01:.01:.13)
    value(x)=sum(x.^3)+x[1]*x[2]
    expected_g=3 .*t.^2;expected_g[1]+=t[2];expected_g[2]+=t[1]
    expected_h=6 .*t.*d;expected_h[1]+=d[2];expected_h[2]+=d[1]
    actual=sidis_job_derivatives(value,t,d)
    @test isapprox(actual.gradient,expected_g;rtol=1e-13)
    @test isapprox(actual.HVP,expected_h;rtol=1e-13)
    @test_throws ErrorException sidis_job_derivatives(value,t[1:12],d[1:12])
end

@testset "Independent decimal reader rejects all13 sample corruption" begin
    reference=synthetic_sidis_job_reference_records(7)
    policy=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT)
    ctx=(contract=policy,)
    packet=Dict("candidate_derivatives"=>Dict("gradient"=>Float64.(reference.J),"HVP"=>Float64.(reference.H),
        "direction"=>Float64.(reference.direction)),"centre_errors_decimal"=>reference.errors,
        "reference_centre"=>Dict("gradient"=>CollinsFit.sidis_job_decimal_vector(reference.J),
            "HVP"=>CollinsFit.sidis_job_decimal_vector(reference.H),"value"=>string(reference.centre_value)))
    previous=nothing
    for row in reference.rows
        @test HistoricalSIDISJobRecords.recompute_reference_step(row,packet,ctx;previous)
        previous=row
    end
    changed=deepcopy(reference.rows[1]);changed["coordinate_value_samples"][13][5]="100.0"
    @test_throws ErrorException HistoricalSIDISJobRecords.recompute_reference_step(changed,packet,ctx)
    changed=deepcopy(reference.rows[1]);changed["directional_Jacobian_samples"][5][13]="100.0"
    @test_throws ErrorException HistoricalSIDISJobRecords.recompute_reference_step(changed,packet,ctx)
    shifted=deepcopy(reference.rows[1])
    setprecision(256) do
        shifted["coordinate_value_samples"][13]=string.(parse.(BigFloat,shifted["coordinate_value_samples"][13]).+10)
    end
    @test_throws ErrorException HistoricalSIDISJobRecords.recompute_reference_step(shifted,packet,ctx)
    for name in ("run_historical_sidis_point_job.jl","HistoricalSIDISJobRecords.jl")
        @test Meta.parseall(read(joinpath(@__DIR__,"../bin",name),String)) isa Expr
    end
end

@testset "Explicitly synthetic all132/all13 global derivative refinement" begin
    ordinary=synthetic_sidis_job_reference_records(7)
    adjacency_limited=synthetic_sidis_job_reference_records(3000)
    stiff=synthetic_sidis_job_reference_records(30000)
    unresolved=synthetic_sidis_job_reference_records(300000)
    @test length(ordinary.rows[1]["coordinate_value_samples"])==13
    @test all(length(row)==5 for row in ordinary.rows[1]["coordinate_value_samples"])
    @test length(ordinary.rows[1]["directional_Jacobian_samples"])==5
    @test ordinary.rows[1]["adjacent_level"]===nothing
    @test ordinary.rows[3]["adjacent_level"]==2
    ids=["SYNTHETIC_CLOSED_LAYOUT:"*sidis_job_point_key(i) for i in 1:132]
    packet(id,phase,reference,levels)=Dict{String,Any}("phase"=>phase,"point"=>Dict("id"=>id),
        "native_level"=>3,"reference_steps"=>deepcopy(reference.rows[levels]),"centre_errors_decimal"=>deepcopy(reference.errors))
    initial=[packet(id,"tailref123",ordinary,1:3) for id in ids]
    @test sidis_job_derivative_decision(ids,initial,3)["selected_reference_level"]==3
    @test sidis_job_derivative_decision(ids,initial,3)["production_ready"]===false
    initial[end]=packet(ids[end],"tailref123",stiff,1:3)
    @test sidis_job_derivative_decision(ids,initial,3)["status"]=="DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"
    fourth=[packet(id,"reference4",ordinary,4:4) for id in ids]
    fourth[end]=packet(ids[end],"reference4",stiff,4:4)
    @test sidis_job_derivative_decision(ids,initial,3;fourth)["status"]=="FAIL"
    initial[end]=packet(ids[end],"tailref123",adjacency_limited,1:3)
    fourth[end]=packet(ids[end],"reference4",adjacency_limited,4:4)
    @test sidis_job_derivative_decision(ids,initial,3)["status"]=="DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"
    @test parse(BigFloat,adjacency_limited.rows[3]["errors_decimal"]["candidate_gradient_vs_inward"])<big"1e-7"
    @test parse(BigFloat,adjacency_limited.rows[3]["errors_decimal"]["adjacent_gradient"])>big"1e-7"
    @test sidis_job_derivative_decision(ids,initial,3;fourth)["selected_reference_level"]==4
    initial[end]=packet(ids[end],"tailref123",unresolved,1:3)
    fourth[end]=packet(ids[end],"reference4",unresolved,4:4)
    @test sidis_job_derivative_decision(ids,initial,3;fourth)["status"]=="FAIL"
    @test_throws ErrorException sidis_job_derivative_decision(ids,initial,3.0)
    @test_throws ErrorException sidis_job_derivative_decision(ids,initial[1:end-1],3)
    mixed=deepcopy(initial);mixed[end]["native_level"]=4
    @test_throws ErrorException sidis_job_derivative_decision(ids,mixed,3)
    missing=deepcopy(initial);pop!(missing[end]["reference_steps"][3]["inward_gradient"])
    @test_throws ErrorException sidis_job_derivative_decision(ids,missing,3)
    @test_throws ErrorException CollinsFit.sidis_job_parse_decimal("NaN")
end

@testset "Explicitly synthetic132 source-angle and target controls; no physical loader" begin
    contract=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT)
    points=[Dict("id"=>"SYNTHETIC_COUNT_LAYOUT:"*sidis_job_point_key(i),
        "dataset_id"=>"SYNTHETIC_COUNT_VIEW","publication"=>i<=40 ? "1006.4221v1" : i<=86 ? "1408.4405" : "0802.2160",
        "kinematics"=>Dict("y"=>.2+i*.001)) for i in 1:132]
    ctx=(contract=contract,request=(data=Dict("points"=>points),))
    report=Dict{String,Any}()
    HistoricalSIDISSourceMath.populate!(report,ctx)
    review=HistoricalSIDISSourceMath.audit(report,ctx)
    @test review["source_factor_errors_recomputed"]==264
    @test review["closed_functions_reevaluated"]==14
    @test review["target_sum_records_recomputed"]==2
    @test report["physical_predictions_executed"]===false
    @test last(report["synthetic_width_refinements"])["accepted"]===true
    changed=deepcopy(report);changed["source_factor_records"][132]["implementation"]*=2
    @test_throws ErrorException HistoricalSIDISSourceMath.audit(changed,ctx)
    changed=deepcopy(report);changed["synthetic_target_records"][2]["ratio"]=(.125-.0625)/2
    @test_throws ErrorException HistoricalSIDISSourceMath.audit(changed,ctx)
end
