using Test, JSON
include("../bin/SIASourceControls.jl")
const C=SIASourceControls
module InheritedSyntheticSourceControls
include("../bin/predict_historical_sia_integrals.jl")
end
const L=InheritedSyntheticSourceControls

function fixture_source_data()
    points=Any[]
    for view in (1002,1003),i in 1:36
        push!(points,Dict("id"=>"KPSY15_SIA_$view:r$(i+1)","ratio"=>view==1002 ? "UC" : "UL",
            "kinematics"=>Dict("z1"=>.2+.1*((i-1)÷6),"z2"=>.2+.1*((i-1)%6),"Q2_GeV2"=>112.36),
            "source_angular"=>Dict("value"=>(view==1002 ? .55 : .65)+.001*i)))
    end
    merge(C.S.false_flags(),Dict{String,Any}("scope"=>"historical_reconstructed_sia_integral_prediction_diagnostic",
        "theory_scope"=>"synthetic_closure","fixture_scope"=>"explicitly_synthetic_arithmetic_fixture_not_primary_data",
        "source_fingerprint"=>"01be55cc1449679f6dc26bd43ea9788604938c1fecbb25d276acaa35bc4bae47","points"=>points))
end

@testset "Synthetic SIA source controls retain scope and independent saved arithmetic" begin
    parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))
    data=fixture_source_data();before=deepcopy(data)
    binding=Dict("fixture_id"=>"explicitly_synthetic_source_control_unit_test_not_primary_data")
    report=C.produce(data,parent,binding,L.HSP.source_projection,L.HIR.mathematical_controls)
    @test report["status"]==C.COMPUTED
    @test report["theory_scope"]=="synthetic_closure"
    @test length(report["source_projection"])==144
    @test length(report["stencil_formula_control"]["records"])==14
    @test length(report["checks"])==433
    @test all(row["passed"] for row in report["checks"])
    review=C.audit(data,parent,binding,JSON.parse(JSON.json(report)))
    @test review["status"]==C.REVIEWED
    @test review["stored_arithmetic_checks_recomputed"]==433
    @test review["Fourier_quadrature_reexecuted_by_reader"]===false
    @test review["physical_SIA_evolution_or_integrals_or_AD_executed"]===false
    @test data==before
    for flag in C.FLAGS
        @test report[flag]===review[flag]===false
        changed=deepcopy(report);changed[flag]=true
        @test_throws ErrorException C.audit(data,parent,binding,changed)
    end
    for mutation in ("scope","schema","context","source-point","omit-angle","sign","angle","projection",
        "counterfactual","positivity","distinctness","saved-error","omit-check","check-decision","stencil-sample",
        "stencil-bias","stencil-aggregate","boolean-constant")
        changed=deepcopy(report);row=changed["source_projection"][end]
        if mutation=="scope";changed["theory_scope"]="reconstructed_theory"
        elseif mutation=="schema";changed["schema_version"]=true
        elseif mutation=="context";changed["source_context"]["fixture_id"]="substituted"
        elseif mutation=="source-point";changed["request_points"][end]["id"]="substituted"
        elseif mutation=="omit-angle";pop!(changed["source_projection"])
        elseif mutation=="sign";row["sign"]=-1
        elseif mutation=="angle";row["angular_factor"]*=.9
        elseif mutation=="projection";row["projected"]+=.01
        elseif mutation=="counterfactual";row["counterfactuals"]["cos_phi0"]+=.01
        elseif mutation=="positivity";row["positive_counts_proved_by_bound"]=false
        elseif mutation=="distinctness";row["all_counterfactuals_differ_from_correct"]=false
        elseif mutation=="saved-error";row["error"]+=.01
        elseif mutation=="omit-check";pop!(changed["checks"])
        elseif mutation=="check-decision";changed["checks"][1]["passed"]=false
        elseif mutation=="stencil-sample";changed["stencil_formula_control"]["records"][1]["samples"][1]="0.0"
        elseif mutation=="stencil-bias";changed["stencil_formula_control"]["records"][11]["known_stencil_bias"]="0.0"
        elseif mutation=="stencil-aggregate";changed["stencil_formula_control"]["passed"]=false
        elseif mutation=="boolean-constant";row["D_L"]=true
        end
        @test_throws ErrorException C.audit(data,parent,binding,changed)
    end
    changed=deepcopy(parent);changed["numerical_policy"]["source_projection_max_scaled"]=1.
    @test_throws ErrorException C.produce(data,changed,binding,L.HSP.source_projection,L.HIR.mathematical_controls)
    @test_throws ErrorException C.audit(data,changed,binding,report)
    changeddata=deepcopy(data);pop!(changeddata["points"])
    @test_throws ErrorException C.produce(changeddata,parent,binding,L.HSP.source_projection,L.HIR.mathematical_controls)
end
