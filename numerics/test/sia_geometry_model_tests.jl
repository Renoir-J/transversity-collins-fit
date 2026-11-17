using Test, CollinsFit, JSON, SHA, LinearAlgebra, ForwardDiff
include("../bin/SIAGeometryModel.jl")
const M=SIAGeometryModel
module InheritedGeometryModelComparator
include("../bin/predict_historical_sia_integrals.jl")
end
const P=InheritedGeometryModelComparator

include("sia_geometry_fixtures.jl")

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
