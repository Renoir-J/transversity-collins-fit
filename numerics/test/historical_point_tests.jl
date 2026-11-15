using SHA

# Representation fixtures only. These invented scalars never replay a primary
# source, invoke the real producer, or confer an experimental/numerical gate.
function historical_layout_fixture(directory)
    root=normpath(joinpath(@__DIR__,"../.."))
    contractfile=joinpath(root,"config/historical-point-contract.json")
    initialfile=joinpath(root,"config/initial-contract.json")
    contract=JSON.parsefile(contractfile);initial=JSON.parsefile(initialfile)
    views=Dict{String,Any}[];points=Dict{String,Any}[]
    for selected in contract["selection"]["source_views"]
        rows=Dict{String,Any}[]
        for i in 1:4
            means=Dict("x"=>.15,"y"=>.4+.01(i-1),"z"=>.3,"Q2_GeV2"=>4.,"pT_GeV"=>.2)
            row=Dict{String,Any}("id"=>selected["id"]*":r"*string(i+1),"averages"=>means,
                "central_value"=>(-1)^i*.0125,"statistical_error"=>.05,
                "systematic_components"=>Dict("experimental"=>.002,"model"=>.003),
                "locator"=>Dict("sheet"=>"synthetic layout only","excel_row"=>i+1),
                "source_remap"=>Dict("kind"=>"reviewed_legacy_column_remap",
                    "primary_tables"=>["table:kine","table:results"],"primary_lines"=>[10+i,30+i]))
            push!(rows,row)
            push!(points,Dict("id"=>row["id"],"dataset_id"=>selected["id"],
                "publication"=>contract["selection"]["publication"],"source_sha256"=>selected["sha256"],
                "target"=>"neutron","target_material"=>"3He","target_result"=>contract["selection"]["target_result"],
                "hadron"=>selected["hadron"],"charge"=>selected["charge"],"kinematics"=>deepcopy(means),
                "measured"=>Dict(key=>deepcopy(row[key]) for key in ("central_value","statistical_error","systematic_components")),
                "source_locator"=>deepcopy(row["locator"]),"source_remap"=>deepcopy(row["source_remap"]),
                "original_depolarization_removed"=>nothing,
                "source_blockers"=>["synthetic layout; no primary replay or experimental admission"]))
        end
        push!(views,Dict("id"=>selected["id"],"source_id"=>selected["id"],"process"=>"SIDIS","projection"=>"x",
            "publication"=>contract["selection"]["publication"],"target"=>"neutron","target_material"=>"3He",
            "target_result"=>contract["selection"]["target_result"],"source_sha256"=>selected["sha256"],
            "hadron"=>selected["hadron"],"source_observable"=>"AUTcollins","output_units"=>"fraction",
            "production_ready"=>false,"historical_fit_ready"=>false,"depolarization_removed"=>nothing,
            "blockers"=>copy(last(points)["source_blockers"]),"points"=>rows,"excluded"=>[]))
    end
    fingerprint="explicitly_synthetic_layout_no_source_replay"
    candidate=Dict("scope"=>"historical_benchmark_candidate","production_ready"=>false,
        "physical_extraction"=>false,"historical_fit_ready"=>false,"source_fingerprint"=>fingerprint,"views"=>views)
    report=Dict("status"=>"CANDIDATES_VALIDATED","production_ready"=>false,"physical_extraction"=>false,
        "historical_fit_ready"=>false,"source_fingerprint"=>fingerprint)
    candidatefile=joinpath(directory,"candidates.json");reportfile=joinpath(directory,"report.json")
    sourcefile=joinpath(directory,"synthetic-layout.txt")
    write(sourcefile,"invented source layout only, never experimental evidence\n")
    writejson(candidatefile,candidate);writejson(reportfile,report)
    digest(path)=bytes2hex(sha256(read(path)))
    data=Dict{String,Any}("schema_version"=>1,"scope"=>contract["scope"],"theory_scope"=>"reconstructed_theory",
        "observation_scope"=>"historical_benchmark_candidate","production_ready"=>false,"physical_extraction"=>false,
        "historical_fit_ready"=>false,"experimental_data_admitted"=>false,"fit_performed"=>false,
        "contract_id"=>contract["id"],"contract_sha256"=>digest(contractfile),"initial_contract_sha256"=>digest(initialfile),
        "source_fingerprint"=>fingerprint,"parameters"=>Dict("order"=>initial["parameter_order"],"values"=>initial["published_parameters"]),
        "model_convention"=>Dict(key=>contract["prediction"][key] for key in ("phase","normalization","model_y","source_mapping_status")),
        "source_preparation"=>Dict("status"=>"CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED",
            "scope"=>"historical_benchmark_candidate","candidates"=>directory,
            "verification_kind"=>"Python_source_preparation_attestation_not_Julia_raw_table_replay",
            "candidate_sha256"=>digest(candidatefile),"report_sha256"=>digest(reportfile)),
        "points"=>points,"source_inputs"=>Dict(path=>digest(path) for path in
            (candidatefile,reportfile,sourcefile,contractfile,initialfile)))
    requestfile=joinpath(directory,"request.json");writejson(requestfile,data)
    (;requestfile,sourcefile,data,candidate,candidatefile)
end

@testset "synthetic historical request layout preserves source identity and rejects admission" begin
    mktempdir() do directory
        f=historical_layout_fixture(directory);request=load_historical_point_request(f.requestfile)
        @test length(request.data["points"])==8
        @test length(unique(historical_geometry_key.(request.data["points"])))==4
        @test [p["charge"] for p in request.data["points"]]==vcat(fill(1,4),fill(-1,4))
        @test request.data["points"][1]["measured"]["systematic_components"]==Dict("experimental"=>.002,"model"=>.003)
        @test request.data["points"][1]["measured"]["central_value"]==-.0125
        @test all(p["original_depolarization_removed"]===nothing for p in request.data["points"])
        @test historical_request_valid(request)===request
        request.data["points"][1]["kinematics"]["x"]+=.001
        @test_throws ErrorException historical_request_valid(request)
    end
    mutations=(d->d["scope"]="production",d->d["scope"]="synthetic_closure",
        d->delete!(d,"scope"),d->d["schema_version"]=true,d->d["production_ready"]=true,
        d->d["physical_extraction"]=true,d->d["historical_fit_ready"]=true,
        d->d["experimental_data_admitted"]=true,d->d["fit_performed"]=true,
        d->d["contract_sha256"]="0"^64,d->d["parameters"]["values"][1]+=.01,
        d->reverse!(d["parameters"]["order"]),d->d["model_convention"]["phase"]="minus_sin_phi_sum",
        d->d["model_convention"]["normalization"]="DNN_removed",
        d->d["source_preparation"]["verification_kind"]="Julia_primary_replay_claim",
        d->pop!(d["points"]),d->reverse!(d["points"]),d->d["points"][1]["charge"]=-1,
        d->d["points"][1]["target"]="proton",d->d["points"][1]["original_depolarization_removed"]=false,
        d->d["points"][1]["kinematics"]["x"]+=.001,
        d->d["points"][1]["measured"]["central_value"]*=2,
        d->delete!(d["points"][1]["measured"]["systematic_components"],"model"),
        d->empty!(d["points"][1]["source_blockers"]),d->d["source_inputs"]=Dict())
    for mutate in mutations
        mktempdir() do directory
            f=historical_layout_fixture(directory);mutate(f.data);writejson(f.requestfile,f.data)
            @test_throws Exception load_historical_point_request(f.requestfile)
        end
    end
    for mutation in (:source_bytes,:missing_source,:request_bytes)
        mktempdir() do directory
            f=historical_layout_fixture(directory);request=load_historical_point_request(f.requestfile)
            mutation==:source_bytes && write(f.sourcefile,"changed synthetic source bytes")
            mutation==:missing_source && rm(f.sourcefile)
            mutation==:request_bytes && write(f.requestfile,read(f.requestfile,String)*"\n")
            @test_throws ErrorException historical_request_valid(request)
        end
    end
end

@testset "historical adapter uses existing predictions without a likelihood on synthetic rows" begin
    mktempdir() do directory
        f=historical_layout_fixture(directory);request=load_historical_point_request(f.requestfile)
        base=observable_fixture(:SIDIS)
        neutron=polarized_leptonic_plan(base.leptonic,base.unpolarized_first,base.unpolarized_second,
            base.first_operator,base.first_samples;second_operator=base.second_operator,
            second_samples=base.second_samples,target=NucleonTarget(:neutron),charges=(1,0))
        contexts=Dict(historical_geometry_key(point)=>neutron for point in request.data["points"])
        model=historical_point_model(request,contexts);theta=collect(kpsy_published_parameters().values)
        actual=predict_joint(model,theta)
        @test length(actual)==8 && all(isfinite,actual)
        for (point,value,observation) in zip(request.data["points"],actual,model.observations)
            y=point["kinematics"]["y"];row=polarized_structure(with_pion_charges(neutron,point["charge"]),theta)
            @test value ≈ 2(1-y)/(1+(1-y)^2)*row.Collins/row.UU rtol=1e-13
            @test observation.plan.context.first_operator===neutron.first_operator
        end
        @test_throws ErrorException historical_point_model(request,Dict())
        wrong=Dict(key=>base for key in keys(contexts))
        @test_throws ErrorException historical_point_model(request,wrong)
        # The general optimizer still rejects this real-observation diagnostic scope.
        @test_throws ErrorException JointGaussianLikelihood(model,model.point_ids,ones(8),Matrix{Float64}(I,8,8);
            scope=request.data["scope"],source="synthetic representation rejection control")
    end
end
