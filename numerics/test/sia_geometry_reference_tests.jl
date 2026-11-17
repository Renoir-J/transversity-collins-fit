using Test, CollinsFit, JSON, SHA
include("../bin/SIAGeometryReferenceJobs.jl")
const R=SIAGeometryReferenceJobs
const M=R.M
module InheritedReferenceGeometryComparator
include("../bin/predict_historical_sia_integrals.jl")
include("../bin/HistoricalSIAIntegralAudit.jl")
end
const P=InheritedReferenceGeometryComparator
include("sia_geometry_fixtures.jl")

@testset "Synthetic geometry tails, scalar operands and representative all13 derivatives" begin
    mktempdir() do directory
        full=geometry_fixture_request(directory);before=deepcopy(full.data)
        parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))
        groups,keys=M.S.geometry_groups(full.data)
        contexts=Dict(keys[i]=>geometry_fixture_context(keys[i][1],keys[i][2]) for i in (1,15,36))
        binding=Dict("fixture_scope"=>"explicitly_synthetic_reference_geometry_jobs_not_physical_data","full_request_sha256"=>full.sha256)
        function prepare(view,parent,level,B)
            key=historical_sia_integral_geometry_key(first(view.data["points"]));context=contexts[key]
            context,Dict("geometry"=>collect(key),"bound_inverse_GeV"=>B,"radial_nodes"=>length(context.leptonic.radial.b),
                "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes,
                "cold_seconds_informational"=>0.,"context_provenance"=>deepcopy(context.provenance))
        end
        tails=Any[];derivatives=Any[]
        inward=P.HistoricalSIAIntegralAudit.inward_checks
        for index in (1,15,36)
            view=M.SelectedSIAGeometry(full,index);level=index==36 ? 4 : 3
            spec=R.tail_definition(full.data,parent,index,level)
            report=R.compute_tail(view,parent,spec,binding,P,prepare)
            audit=R.audit_tail(full.data,parent,spec,binding,JSON.parse(JSON.json(report)))
            @test report["status"]==R.TAIL_COMPUTED && audit["status"]==R.TAIL_REVIEWED
            @test audit["saved_value_points_recomputed"]==20
            @test audit["saved_scalar_reference_comparisons_recomputed"]==10
            @test audit["fixed_value_control_comparisons_recomputed"]==10
            @test audit["infinite_tail_proved"]===false && audit["job_resource_and_complete_provenance_ledger_verified"]===false
            @test [value["bound_inverse_GeV"] for value in report["tail_values"]]==[12,32]
            @test all([value["order"] for value in reference["order_comparisons"]]==[16,32,64,128,256] for reference in report["independent_references"])
            push!(tails,report)
            derivative_spec=R.derivative_definition(full.data,parent,index,level,64)
            derivative=R.compute_derivative(view,parent,derivative_spec,binding,P,prepare)
            derivative_audit=R.audit_derivative(full.data,parent,derivative_spec,binding,JSON.parse(JSON.json(derivative)),inward)
            @test derivative["status"]==R.DERIVATIVE_COMPUTED && derivative_audit["status"]==R.DERIVATIVE_REVIEWED
            @test derivative_audit["stored_derivative_comparisons_recomputed"]==9
            @test derivative_audit["all72_production_derivative_certificate"]===false
            @test derivative_audit["independent_cold_AD_reexecuted_by_reader"]===false
            @test derivative["derivative_record"]["point_ids"]==[full.data["points"][groups[index][2]]["id"],full.data["points"][groups[index][1]]["id"]]
            push!(derivatives,derivative)
        end
        @test full.data==before
        tail=last(tails);spec=tail["definition"]
        for flag in R.N.FLAGS
            @test tail[flag]===false
            bad=deepcopy(tail);bad[flag]=true
            @test_throws ErrorException R.audit_tail(full.data,parent,spec,binding,bad)
        end
        for mutation in ("missing-bound","B24-substitution","wrong-resource-bound","source-angular","reference-scope","evolution-claim",
            "omit-reference","reference-order","reference-error","reference-pair","reference-category","reference-algebra","reference-bool","reference-nan",
            "omit-fixed-order","fixed-theta","fixed-output","fixed-error")
            bad=deepcopy(tail)
            if mutation=="missing-bound";pop!(bad["tail_values"])
            elseif mutation=="B24-substitution";bad["tail_values"][1]["bound_inverse_GeV"]=24
            elseif mutation=="wrong-resource-bound";bad["tail_values"][2]["resources"][1]["bound_inverse_GeV"]=24
            elseif mutation=="source-angular";bad["independent_references"][1]["angular_value"]*=.9
            elseif mutation=="reference-scope";bad["independent_references"][1]["scope"]="independent_evolution"
            elseif mutation=="evolution-claim";bad["independent_references"][1]["independent_evolution"]=true
            elseif mutation=="omit-reference";pop!(bad["independent_references"])
            elseif mutation=="reference-order";reverse!(bad["independent_references"][1]["order_comparisons"])
            elseif mutation=="reference-error";bad["independent_references"][1]["order_comparisons"][1]["max_scaled"]+=1.
            elseif mutation=="reference-pair";bad["independent_references"][1]["reference"]["pairs"]["pp"]["UU"]+=1.
            elseif mutation=="reference-category";bad["independent_references"][1]["reference"]["categories"]["U"]["UU"]+=1.
            elseif mutation=="reference-algebra";bad["independent_references"][1]["reference"]["amplitude"]+=.01
            elseif mutation=="reference-bool"
                pair=Dict{String,Any}(bad["independent_references"][1]["reference"]["pairs"]["pp"]);pair["UU"]=true
                pairs=Dict{String,Any}(bad["independent_references"][1]["reference"]["pairs"]);pairs["pp"]=pair
                bad["independent_references"][1]["reference"]["pairs"]=pairs
            elseif mutation=="reference-nan";bad["independent_references"][1]["reference"]["amplitude"]=NaN
            elseif mutation=="omit-fixed-order";pop!(bad["fixed_value_controls"])
            elseif mutation=="fixed-theta";bad["fixed_value_controls"][1]["flipped_theta"][7]*=.9
            elseif mutation=="fixed-output";bad["fixed_value_controls"][1]["flipped"][1]+=.01
            elseif mutation=="fixed-error";bad["fixed_value_controls"][1]["flipped_max_scaled"]+=.01
            end
            @test_throws ErrorException R.audit_tail(full.data,parent,spec,binding,bad)
        end
        derivative=last(derivatives);spec=derivative["definition"]
        for mutation in ("source-id","direction","theta","scope","zero-column","stencil-sample","saved-error","input-kind","resource-bound","value-angle")
            bad=deepcopy(derivative);record=bad["derivative_record"]
            if mutation=="source-id";record["geometry_source_id"]="KPSY15_SIA_1002:r2"
            elseif mutation=="direction";record["direction"][13]*=.9
            elseif mutation=="theta";record["theta"][13]*=.9
            elseif mutation=="scope";record["scope"]="all72_production_certificate"
            elseif mutation=="zero-column";record["first_six_gradient_columns_zero"]=false
            elseif mutation=="stencil-sample";record["finite_differences"][1]["gradient_samples_decimal"][7][1][1]="1.0"
            elseif mutation=="saved-error";record["independent_HVP_error"]+=1.
            elseif mutation=="input-kind";bad["initial_source_kind"]="pinned_inputs"
            elseif mutation=="resource-bound";bad["resources"][1]["bound_inverse_GeV"]=24
            elseif mutation=="value-angle";bad["value_points"][1]["angular_factor"]*=.9
            end
            @test_throws ErrorException R.audit_derivative(full.data,parent,spec,binding,bad,inward)
        end
        for (index,level,order) in ((2,3,64),(1,2,64),(1,3,32),(1,3,true))
            @test_throws ErrorException R.derivative_definition(full.data,parent,index,level,order)
        end
        partial=Ref{Any}(nothing);view=M.SelectedSIAGeometry(full,36);spec=tail["definition"]
        progress=stage->get(stage,"bound_inverse_GeV",nothing)===32 && get(stage,"momentum_order",nothing)===128 ?
            error("explicit synthetic finite-tail interruption") : nothing
        @test_throws ErrorException R.compute_tail(view,parent,spec,binding,P,prepare;progress=progress,checkpoint=r->(partial[]=deepcopy(r)))
        @test partial[]["status"]=="IN_PROGRESS_NOT_A_GATE" && length(partial[]["tail_values"])==2
        @test length(partial[]["tail_values"][1]["momentum_levels"])==5 && length(partial[]["tail_values"][2]["momentum_levels"])==3
        @test partial[]["accepted_fitting_gate"]===false
    end
end
