include("../bin/PaperMoments.jl")
using .PaperMoments

@testset "synthetic partial-moment integration and source-sensitive checks" begin
    # An affine analytic oracle catches an omitted dx=x*dlog(x) Jacobian and
    # partial first/last cell mistakes without copying the integration algorithm.
    grid=[.0065,.02,.08,.3,.6,1.]
    lower,upper=.011,.35
    weights=moment_weights(grid,lower,upper)
    values=2 .- 3grid
    exact=2(upper-lower)-1.5(upper^2-lower^2)
    @test all(>=(0),weights)
    @test sum(weights) ≈ upper-lower atol=2e-16
    @test dot(weights,values) ≈ exact atol=2e-16
    @test dot(moment_weights(grid,first(grid),last(grid)),ones(length(grid))) ≈ 1-first(grid)
    quadrature=reference_integral(grid,values,lower,upper)
    @test quadrature.value ≈ exact atol=3e-15
    @test quadrature.estimated_error < 1e-12
    @test_throws ErrorException moment_weights(grid,.001,.35)
    @test_throws ErrorException moment_weights(grid,.01,1.01)
    @test_throws ErrorException moment_weights(grid,.35,.01)
    @test_throws ErrorException moment_weights(grid,NaN,.35)
    @test_throws ErrorException reference_integral(grid,[1.],lower,upper)

    theta=collect(kpsy_published_parameters().values)
    soffer=[2. 1. 0.;1. .5 0.]
    x=[.02,.2,1.]
    reference=direct_initial(x,soffer,theta)
    direct=[initial_transversity(KPSYParameters(theta),pid,v,soffer[f,j]) for (j,v) in enumerate(x),(f,pid) in enumerate((2,1))]
    @test reference ≈ direct rtol=1e-14
    @test reference[end,:]==[0.,0.]
    @test all(>=(0),reference[:,1]) && all(<=(0),reference[:,2])

    contract=JSON.parsefile(joinpath(@__DIR__,"../../config/kpsy-paper-moment-contract.json"))
    bounds=contract["printed_value_policy"]["parameter_rounding"]["bounds"]
    vertices=rounding_vertices(theta,bounds)
    @test length(vertices)==64
    @test length(unique(Tuple(t[1:6]) for t in vertices))==64
    @test all(t->t[7:13]==theta[7:13],vertices)
    @test extrema(t[2] for t in vertices)==(-1.,-.95)
    @test all(t->KPSYParameters(t) isa KPSYParameters,vertices)
    @test_throws ErrorException rounding_vertices(theta,[[NaN,1.] for _ in 1:6])

    @test compare_check("equal",[.3],[.3],0.)["passed"]
    failure=compare_check("paper mismatch",[.28],[.30],.005;metric=:absolute)
    @test !failure["passed"]
    @test failure["actual"]==[.28] && failure["reference"]==[.30]
    @test !compare_check("missing finite result",[NaN],[.3],.005;metric=:absolute)["passed"]
    @test_throws ErrorException compare_check("missing",Float64[],Float64[],1.)
    @test_throws ErrorException compare_check("bad tolerance",[1.],[1.],-1.)
    @test_throws ErrorException compare_check("bad metric",[1.],[1.],1.;metric=:unknown)
end
