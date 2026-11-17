using Test, CollinsFit, JSON, SHA
include("../bin/SIAGeometryNativeJobs.jl")
const N=SIAGeometryNativeJobs
const M=N.M
module InheritedNativeGeometryComparator
include("../bin/predict_historical_sia_integrals.jl")
end
const P=InheritedNativeGeometryComparator
include("sia_geometry_fixtures.jl")

@testset "Synthetic native geometry computations preserve all36 and all momentum orders" begin
    mktempdir() do directory
        full=geometry_fixture_request(directory);before=deepcopy(full.data)
        parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))
        groups,keys=M.S.geometry_groups(full.data)
        contexts=Dict(key=>geometry_fixture_context(key[1],key[2]) for key in keys)
        binding=Dict("fixture_scope"=>"explicitly_synthetic_native_geometry_jobs_not_physical_data",
            "full_request_sha256"=>full.sha256)
        calls=Any[]
        function prepare(view,parent,level)
            push!(calls,(view.index,level))
            key=historical_sia_integral_geometry_key(first(view.data["points"]));context=contexts[key]
            context,Dict("geometry"=>collect(key),"radial_nodes"=>length(context.leptonic.radial.b),
                "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes,
                "cold_seconds_informational"=>0.,"context_provenance"=>deepcopy(context.provenance))
        end
        reports=Any[]
        for index in 1:36
            view=M.SelectedSIAGeometry(full,index);spec=N.definition(full.data,parent,"native123",index)
            report=N.compute(view,parent,spec,binding,P,prepare)
            audit=N.audit(full.data,parent,spec,binding,JSON.parse(JSON.json(report)))
            @test report["status"]==N.COMPUTED && audit["status"]==N.REVIEWED
            @test audit["saved_value_points_recomputed"]==30
            @test audit["global_common_order_or_native_selection_performed"]===false
            @test audit["job_resource_and_complete_provenance_ledger_verified"]===false
            @test [value["native_level"] for value in report["native_values"]]==[1,2,3]
            @test all([row["order"] for row in value["momentum_levels"]]==[16,32,64,128,256] for value in report["native_values"])
            push!(reports,report)
        end
        @test calls==[(index,level) for index in 1:36 for level in 1:3]
        native=[[report["native_values"][level] for report in reports] for level in 1:3]
        global_selection=N.S.select_native(full.data,native,parent)
        @test global_selection["selected_level"]==3
        @test all(result["selected_order"]==64 for result in global_selection["levels"])
        @test length(last(global_selection["levels"])["points"])==72
        @test_throws ErrorException N.S.select_native(full.data,[level[1:35] for level in native],parent)
        for index in (1,15,36)
            view=M.SelectedSIAGeometry(full,index);spec=N.definition(full.data,parent,"native4",index)
            report=N.compute(view,parent,spec,binding,P,prepare);audit=N.audit(full.data,parent,spec,binding,report)
            @test audit["saved_value_points_recomputed"]==10
            @test length(report["native_values"])==1 && only(report["native_values"])["native_level"]===4
        end
        @test full.data==before
        spec=N.definition(full.data,parent,"native123",36);original=last(reports)
        for flag in N.FLAGS
            @test original[flag]===false
            changed=deepcopy(original);changed[flag]=true
            @test_throws ErrorException N.audit(full.data,parent,spec,binding,changed)
        end
        for mutation in ("phase","binding","scope","parameters","point-source","native-omit","native-order","native-setting-type",
            "momentum-omit","momentum-type","point-omit","id","angular","window","category","negative-yield",
            "source-kind","identity-attestation","resource-missing","resource-type","resource-nan")
            changed=deepcopy(original);value=last(changed["native_values"]);row=last(value["momentum_levels"])["points"][2]
            if mutation=="phase";changed["definition"]["phase"]="native4"
            elseif mutation=="binding";changed["source_context"]["full_request_sha256"]="substituted"
            elseif mutation=="scope";changed["theory_scope"]="reconstructed_theory"
            elseif mutation=="parameters";changed["parameters"][7]*=.9
            elseif mutation=="point-source";changed["request_points"][2]["source_angular"]["value"]*=.9
            elseif mutation=="native-omit";pop!(changed["native_values"])
            elseif mutation=="native-order";reverse!(changed["native_values"])
            elseif mutation=="native-setting-type";value["setting"]["radial_order"]=true
            elseif mutation=="momentum-omit";pop!(value["momentum_levels"])
            elseif mutation=="momentum-type";value["momentum_levels"][1]["order"]=16.
            elseif mutation=="point-omit";pop!(last(value["momentum_levels"])["points"])
            elseif mutation=="id";row["id"]=full.data["points"][1]["id"]
            elseif mutation=="angular";row["angular_factor"]*=.9
            elseif mutation=="window";row["mapped_P_bounds_GeV"][2]+=1.
            elseif mutation=="category";row["categories"]["U"]["UU"]+=1.
            elseif mutation=="negative-yield";row["pairs"]["pp"]["UU"]=-1.
            elseif mutation=="source-kind";value["initial_source_kind"]="pinned_inputs"
            elseif mutation=="identity-attestation";value["same_native_rows_across_ratios_producer_attestation"]=false
            elseif mutation=="resource-missing";empty!(value["resources"])
            elseif mutation=="resource-type";value["resources"][1]["radial_nodes"]=true
            elseif mutation=="resource-nan";value["resources"][1]["cold_seconds_informational"]=NaN
            end
            @test_throws ErrorException N.audit(full.data,parent,spec,binding,changed)
        end
        # A common factor cancels in a ratio. The reader must still preserve
        # absolute yields and reject a different pair normalization for UL.
        changed=deepcopy(original);row=last(last(changed["native_values"])["momentum_levels"])["points"][2]
        for collection in (row["pairs"],row["categories"]),value in values(collection),key in ("UU","Collins")
            value[key]*=2
        end
        @test N.S.validate_point(row,full.data["points"][72],256,parent["numerical_policy"])===row
        @test_throws ErrorException N.audit(full.data,parent,spec,binding,changed)
        for (phase,index) in (("tailref",1),("native123",true),("native4",37))
            @test_throws ErrorException N.definition(full.data,parent,phase,index)
        end
        wrong=deepcopy(spec);wrong["geometry_index"]=35
        @test_throws ErrorException N.compute(M.SelectedSIAGeometry(full,36),parent,wrong,binding,P,prepare)
        badtheta=deepcopy(full.data);badtheta["parameters"]["values"]=Any[badtheta["parameters"]["values"]...]
        badtheta["parameters"]["values"][1]=true
        @test_throws ErrorException N.literal_theta(badtheta)
        partial=Ref{Any}(nothing)
        progress=stage->get(stage,"momentum_order",nothing)===128 ? error("explicit synthetic interruption") : nothing
        @test_throws ErrorException N.compute(M.SelectedSIAGeometry(full,36),parent,spec,binding,P,prepare;
            progress=progress,checkpoint=record->(partial[]=deepcopy(record)))
        @test partial[]["status"]=="IN_PROGRESS_NOT_A_GATE"
        @test length(partial[]["native_values"])==1 && length(only(partial[]["native_values"])["momentum_levels"])==3
        @test partial[]["accepted_fitting_gate"]===false
    end
end
