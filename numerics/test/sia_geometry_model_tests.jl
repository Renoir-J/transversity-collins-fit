using Test, CollinsFit, JSON, SHA, LinearAlgebra, ForwardDiff
include("../bin/SIAGeometryModel.jl")
const M=SIAGeometryModel
module InheritedGeometryModelComparator
include("../bin/predict_historical_sia_integrals.jl")
end
const P=InheritedGeometryModelComparator

# The inherited synthetic-row construction is used only for model composition.
# There is no retained native evolution, raw experimental intake or physical job.
function geometry_fixture_context(z1,z2)
    model=KPSYCSS();radial=RadialRule([0.,.25,1.,4.,8.];order=4)
    plan=leptonic_plan(:SIA,model,radial,112.36,0.,z1;z2=z2)
    grid=sort!(unique(vcat(.05,.1,[z1,z2],.75,1.)))
    q2=css_mu_b2.(Ref(model),radial.b);tau=evolution_time.(Ref(model.coupling),q2,model.Q02)
    alpha=alpha_s.(Ref(model.coupling),q2);blocks=CollinsFit.HoppetPointBlock[]
    for query in (z1,z2)
        rows=zeros(length(q2),length(grid));rows[:,findfirst(==(query),grid)].=1 .+ .01 .*radial.b
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),rows))
    end
    payload=(grid,[z1,z2],q2,model.Q02,tau,:KPSY_Hhat3_homogeneous,5,32,:synthetic_fixture,
        blocks,sum(sizeof(b.rows) for b in blocks),0.)
    op=HoppetMatchingPlan(payload...,:SIA,alpha,CollinsFit.hoppet_point_snapshot(payload...;process=:SIA,alpha=alpha))
    samples=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op.grid,
        [(1+.3i)*x*(1-x)^3 for i in 1:4,x in op.grid];input_kind=:analytic)
    first=[(1+.07i)*(1+.03b) for i in 1:10,b in radial.b]
    second=[(1+.11i)*(1+.05b) for i in 1:10,b in radial.b]
    polarized_leptonic_plan(plan,first,second,op,samples;charges=(1,-1))
end

function geometry_fixture_request(directory;name="full72.json",measured=.02,Pmean=.3)
    parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))
    points=Any[]
    for view in (1002,1003),i in 1:36
        push!(points,Dict("id"=>"KPSY15_SIA_$view:r$(i+1)","ratio"=>view==1002 ? "UC" : "UL",
            "publication"=>"explicitly_synthetic_fixture_not_primary_data",
            "kinematics"=>Dict("z1"=>.2+.1*((i-1)÷6),"z2"=>.2+.1*((i-1)%6),"Q2_GeV2"=>112.36),
            "measured"=>Dict("central_value"=>measured),"reported_pT_mean_GeV"=>Pmean,
            "source_angular"=>Dict("value"=>(view==1002 ? .55 : .65)+.001*i,
                "averaging_measure"=>"explicitly synthetic fixed mathematical measure")))
    end
    data=merge(M.S.false_flags(),Dict{String,Any}("scope"=>"historical_reconstructed_sia_integral_prediction_diagnostic",
        "theory_scope"=>"synthetic_closure","fixture_scope"=>"explicitly_synthetic_arithmetic_fixture_not_primary_data",
        "source_fingerprint"=>"01be55cc1449679f6dc26bd43ea9788604938c1fecbb25d276acaa35bc4bae47",
        "points"=>points,"momentum_window"=>CollinsFit.historical_sia_integral_window(parent),"source_inputs"=>Dict{String,String}()))
    path=joinpath(directory,name);writejson(path,data);owned=JSON.parsefile(path)
    HistoricalSIAIntegralRequest(owned,path,bytes2hex(sha256(read(path))),CollinsFit.numerical_snapshot(owned))
end

@testset "Synthetic two-row geometry views match the unchanged full72 model" begin
    mktempdir() do directory
        full=geometry_fixture_request(directory);original=deepcopy(full.data)
        groups,keys=M.S.geometry_groups(full.data)
        contexts=Dict(key=>geometry_fixture_context(key[1],key[2]) for key in keys)
        theta=collect(kpsy_published_parameters().values)
        for order in (16,32,64,128,256)
            whole=P.predictions(full,contexts,theta,order)
            for index in 1:36
                view=M.SelectedSIAGeometry(full,index);selected=Dict(keys[index]=>contexts[keys[index]])
                result=P.predictions(view,selected,theta,order)
                @test length(view.full.data["points"])==72 && length(view.data["points"])==2
                @test !(view isa HistoricalSIAIntegralRequest) && view.data["scope"]==M.VIEW_SCOPE
                @test [r["id"] for r in result.points]==[full.data["points"][i]["id"] for i in groups[index]]
                @test all(M.S.same_json(result.points[j],whole.points[i]) for (j,i) in enumerate(groups[index]))
                @test all(o.plan.context.first_operator===contexts[keys[index]].first_operator for o in result.model.observations)
            end
        end
        @test full.data==original
        other=geometry_fixture_request(directory;name="changed-measurements.json",measured=-999.,Pmean=999.)
        whole_model=historical_sia_integral_model(full,contexts;order=64)
        direction=[.02,.03,.01,.02,.01,.01,.02,-.02,.01,.02,.03,.03,.001]
        for index in (1,15,36)
            selected=Dict(keys[index]=>contexts[keys[index]])
            view=M.SelectedSIAGeometry(full,index);model=historical_sia_integral_model(view,selected;order=64)
            altered=historical_sia_integral_model(M.SelectedSIAGeometry(other,index),selected;order=64)
            @test predict_joint(model,theta)==predict_joint(altered,theta)
            f=t->predict_joint(model,t);g=t->predict_joint(whole_model,t)[groups[index]]
            J=ForwardDiff.jacobian(f,theta);RJ=ForwardDiff.jacobian(g,theta)
            H=ForwardDiff.derivative(a->vec(ForwardDiff.jacobian(f,theta.+a.*direction)),0.)
            RH=ForwardDiff.derivative(a->vec(ForwardDiff.jacobian(g,theta.+a.*direction)),0.)
            @test J==RJ && H==RH
            @test all(iszero,J[:,1:6]) && all(iszero,reshape(H,2,13)[:,1:6])
            @test_throws ErrorException historical_sia_integral_model(view,contexts;order=64)
            for order in (true,8,64.)
                @test_throws ErrorException historical_sia_integral_model(view,selected;order=order)
            end
            @test_throws ErrorException historical_sia_integral_model(view,Dict(keys[index]=>with_momentum(contexts[keys[index]],.1)))
        end
        for index in (true,0,37,1.)
            @test_throws ErrorException M.SelectedSIAGeometry(full,index)
        end
        for field in ("points","momentum_window","geometry_index","full_request_sha256")
            view=M.SelectedSIAGeometry(full,1)
            if field=="points";view.data[field][1]["source_angular"]["value"]*=.9
            elseif field=="momentum_window";view.data[field]["qT_bounds_GeV"][2]=3.4
            elseif field=="geometry_index";view.data[field]=2
            else;view.data[field]="substituted"
            end
            @test_throws ErrorException M.valid(view)
        end
        view=M.SelectedSIAGeometry(full,1);full.data["points"][1]["measured"]["central_value"]+=1
        @test_throws ErrorException M.valid(view)
        @test_throws ErrorException M.SelectedSIAGeometry(full,1)
    end
end
