include("../bin/ObservableChecks.jl")

@testset "independent scalar/adaptive observable reference on synthetic rows" begin
    theta=[.5,-.4,1.,1.,2.,2.,.3,-.2,1.,1.,0.,0.,.01]
    for process in (:SIDIS,:SIA),interval in (nothing,(.1,.3)),
            angular_kind in (process==:SIA ? (:ratio_of_angular_means,:mean_of_angular_ratio) : (:ratio_of_angular_means,))
        context=observable_fixture(process)
        pair=observable_pair(context;interval=interval,angular_kind=angular_kind)
        prediction(t)=observable_values(pair,t)
        reference=prepare_generic_reference(context;momentum_interval=interval,y=OBS_Y,angular=OBS_ANGULAR).predict
        @test prediction(theta) ≈ reference(theta) rtol=1e-12 atol=1e-14
        @test ForwardDiff.jacobian(prediction,theta) ≈ ForwardDiff.jacobian(reference,theta) rtol=1e-11 atol=1e-14
        for i in 1:2
            @test ForwardDiff.hessian(t->prediction(t)[i],theta) ≈ ForwardDiff.hessian(t->reference(t)[i],theta) rtol=1e-10 atol=1e-13
        end
        @test Float64.(reference(BigFloat.(theta))) ≈ prediction(theta) rtol=1e-12 atol=1e-14
    end
end
