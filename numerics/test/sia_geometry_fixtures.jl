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
        "parameters"=>Dict("values"=>collect(kpsy_published_parameters().values)),
        "points"=>points,"momentum_window"=>CollinsFit.historical_sia_integral_window(parent),"source_inputs"=>Dict{String,String}()))
    path=joinpath(directory,name);writejson(path,data);owned=JSON.parsefile(path)
    HistoricalSIAIntegralRequest(owned,path,bytes2hex(sha256(read(path))),CollinsFit.numerical_snapshot(owned))
end
