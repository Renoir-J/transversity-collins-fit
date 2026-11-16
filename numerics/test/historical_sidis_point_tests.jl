using SHA

# Invented layout and analytic rows only. The raw-source loader cannot admit
# these fixtures: they have no pinned primary PDFs/source archives or replay.
function historical_sidis_layout_fixture(directory)
    root=normpath(joinpath(@__DIR__,"../.."));contract=JSON.parsefile(joinpath(root,"config/historical-sidis-point-contract.json"))
    views=Dict{String,Any}[];points=Dict{String,Any}[];index=0
    for selected in contract["selection"]["source_views"]
        rows=Dict{String,Any}[];mode=contract["model_conventions"][selected["publication"]]
        for i in 1:selected["points"]
            index+=1;k=Dict("x"=>.15,"y"=>.3+.001index,"z"=>.3,"Q2_GeV2"=>4.,"pT_GeV"=>.2+index*1e-6)
            row=Dict{String,Any}("id"=>selected["id"]*":r"*string(i+1),"averages"=>k,"bins"=>Dict(),
                "central_value"=>(-1)^i*.0125,"statistical_error"=>.05,
                "systematic_components"=>(mode["id"]=="HERMES2010" ? Dict("experimental"=>.002) : Dict()),
                "locator"=>Dict("sheet"=>"explicitly synthetic layout only","excel_row"=>i+1))
            if mode["id"]=="COMPASS2014"
                row["statistical_error"]=nothing;row["reported_error_u"]=.06
                row["uncertainty_status"]="reported_error_u_composition_unresolved"
            end
            push!(rows,row)
            push!(points,Dict{String,Any}("id"=>row["id"],"dataset_id"=>selected["id"],
                (key=>selected[key] for key in ("publication","target","hadron","projection","source_sha256"))...,
                "charge"=>(selected["hadron"]=="pi+" ? 1 : -1),"kinematics"=>deepcopy(k),"source_row"=>deepcopy(row),
                "measured"=>Dict(key=>deepcopy(row[key]) for key in ("central_value","statistical_error","systematic_components")),
                "original_depolarization_removed"=>selected["original_depolarization_removed"],
                "model_convention"=>deepcopy(mode),"source_blockers"=>["explicitly synthetic layout only"]))
        end
        n=selected["projection"]=="pt" ? (mode["id"]=="HERMES2010" ? 3 : 2) : 0
        push!(views,Dict{String,Any}("id"=>selected["id"],"source_id"=>selected["id"],"process"=>"SIDIS",
            (key=>selected[key] for key in ("publication","target","hadron","projection","source_sha256"))...,
            "source_observable"=>"AUTcollins","output_units"=>"fraction","production_ready"=>false,
            "historical_fit_ready"=>false,"depolarization_removed"=>selected["original_depolarization_removed"],
            "blockers"=>["explicitly synthetic layout only"],"points"=>rows,
            "excluded"=>[Dict("reason"=>"synthetic exclusion layout","excel_row"=>20+i) for i in 1:n]))
    end
    candidate=Dict("views"=>views)
    data=Dict{String,Any}("scope"=>contract["scope"],"source_fingerprint"=>"explicitly_synthetic_layout_not_raw_replay",
        (key=>false for key in CollinsFit.HISTORICAL_SIDIS_FALSE_FLAGS)...,
        "accepted_fitting_gate"=>false,"likelihood_constructed"=>false,
        "points"=>points,"target_weights"=>deepcopy(contract["prediction"]["target_weights"]),
        "view_evidence"=>Dict(v["id"]=>Dict(k=>deepcopy(x) for (k,x) in v if k!="points") for v in views))
    source=joinpath(directory,"synthetic-layout.txt");write(source,"invented layout, never primary evidence\n")
    data["source_inputs"]=Dict(source=>bytes2hex(sha256(read(source))))
    path=joinpath(directory,"request.json");writejson(path,data)
    request=HistoricalSIDISPointRequest(deepcopy(data),path,bytes2hex(sha256(read(path))),CollinsFit.numerical_snapshot(data))
    (;data,candidate,contract,request,path,source)
end

@testset "synthetic full132 SIDIS layouts retain exact means, exclusions and unknown errors" begin
    mktempdir() do directory
        f=historical_sidis_layout_fixture(directory)
        @test CollinsFit._validate_historical_sidis_layout(f.data,f.contract,f.candidate)===nothing
        @test length(unique(historical_sidis_geometry_key.(f.data["points"])))==132
        @test count(p->p["measured"]["statistical_error"]===nothing,f.data["points"])==46
        @test historical_sidis_request_valid(f.request)===f.request
        # The direct toy representation above is never a valid source-loaded request.
        @test_throws Exception load_historical_sidis_point_request(f.path)
        write(f.source,"changed invented layout bytes")
        @test_throws ErrorException historical_sidis_request_valid(f.request)
    end
    for mutate in (d->pop!(d["points"]),d->reverse!(d["points"]),
        d->d["points"][1]["model_convention"]["phase"]="sin_phi_sum_minus_pi",
        d->d["points"][41]["measured"]["statistical_error"]=.06,
        d->d["points"][1]["original_depolarization_removed"]=true,
        d->d["points"][1]["kinematics"]["pT_GeV"]+=1e-5,
        d->empty!(d["view_evidence"]["KPSY15_SIDIS_3026"]["excluded"]))
        mktempdir() do directory
            f=historical_sidis_layout_fixture(directory);mutate(f.data)
            @test_throws Exception CollinsFit._validate_historical_sidis_layout(f.data,f.contract,f.candidate)
        end
    end
end

@testset "synthetic SIDIS composition preserves source factors and p/n yield sums" begin
    mktempdir() do directory
        f=historical_sidis_layout_fixture(directory);base=observable_fixture(:SIDIS)
        targets=Dict(name=>polarized_leptonic_plan(base.leptonic,base.unpolarized_first,base.unpolarized_second,
            base.first_operator,base.first_samples;second_operator=base.second_operator,second_samples=base.second_samples,
            target=historical_sidis_target(name),charges=(1,0)) for name in ("proton","deuteron"))
        neutron=polarized_leptonic_plan(base.leptonic,base.unpolarized_first,base.unpolarized_second,
            base.first_operator,base.first_samples;second_operator=base.second_operator,second_samples=base.second_samples,
            target=NucleonTarget(:neutron),charges=(1,0))
        contexts=Dict(historical_sidis_geometry_key(p)=>with_momentum(targets[p["target"]],p["kinematics"]["pT_GeV"])
            for p in f.data["points"])
        model=historical_sidis_point_model(f.request,contexts);theta=collect(kpsy_published_parameters().values)
        actual=predict_joint(model,theta)
        function source_formula(point,t)
            P=point["kinematics"]["pT_GeV"];charge=point["charge"]
            proton=polarized_structure(with_pion_charges(with_momentum(targets["proton"],P),charge),t)
            if point["target"]=="deuteron"
                n=polarized_structure(with_pion_charges(with_momentum(neutron,P),charge),t)
                ratio=(proton.Collins+n.Collins)/(proton.UU+n.UU)
            else
                ratio=proton.Collins/proton.UU
            end
            y=point["kinematics"]["y"]
            (point["publication"]=="1006.4221v1" ? 2(1-y)/(1+(1-y)^2) : -1)*ratio
        end
        expected=[source_formula(p,theta) for p in f.data["points"]]
        @test length(actual)==132 && all(isfinite,actual)
        @test actual ≈ expected rtol=1e-12
        # This independently builds p/n sums instead of applying the composed
        # deuteron context, and covers all132 source factors and both charges.
        jacobian=ForwardDiff.jacobian(t->predict_joint(model,t),theta)
        reference=ForwardDiff.jacobian(t->[source_formula(p,t) for p in f.data["points"]],theta)
        @test size(jacobian)==(132,13)
        @test jacobian ≈ reference rtol=1e-11
        point=f.data["points"][87];P=point["kinematics"]["pT_GeV"]
        p=polarized_structure(with_pion_charges(with_momentum(targets["proton"],P),point["charge"]),theta)
        n=polarized_structure(with_pion_charges(with_momentum(neutron,P),point["charge"]),theta)
        @test !isapprox((p.Collins+n.Collins)/(p.UU+n.UU),(p.Collins/p.UU+n.Collins/n.UU)/2;rtol=1e-8)
        @test historical_sidis_point_observable(f.request,point["id"],contexts[historical_sidis_geometry_key(point)]).angular==-1
        @test_throws ErrorException historical_sidis_point_observable(f.request,"not-a-source-row",base)
        @test_throws ErrorException historical_sidis_point_model(f.request,Dict())
        wrong=copy(contexts);wrong[historical_sidis_geometry_key(point)]=with_momentum(targets["proton"],P)
        @test_throws ErrorException historical_sidis_point_model(f.request,wrong)
        @test_throws ErrorException JointGaussianLikelihood(model,model.point_ids,zeros(132),Matrix{Float64}(I,132,132);
            scope=f.data["scope"],source="explicitly synthetic diagnostic rejection control")
    end
end
