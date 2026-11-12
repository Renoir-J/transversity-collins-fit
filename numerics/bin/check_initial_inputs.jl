# Pinned-input initial/evolved distribution pilot; not a fit or production receipt.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/initial-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"config/parton-contract.json"),joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function condition_grid(n,minimum,transition)
    nlo=n÷4+1;nhi=n-nlo+1
    low=exp.(range(log(minimum),log(transition);length=nlo));low[1]=minimum;low[end]=transition
    high=1 .- range(sqrt(1-transition),0.;length=nhi).^2
    vcat(low,high[2:end])
end

function interpolate(x,v,target)
    first(x)<=target<=last(x) || error("no input-grid extrapolation")
    i=clamp(searchsortedlast(x,target),1,length(x)-1)
    w=(target-x[i])/(x[i+1]-x[i])
    (1-w)*v[i]+w*v[i+1]
end

function initial_input_checks(buildfile)
    sources=source_hashes();build=CollinsFit.verify_parton_build(buildfile)
    build_hash=bytes2hex(sha256(read(buildfile)))
    policy=JSON.parsefile(CONTRACT)["validation_policy"]
    scaled(a,b)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    pdf=CT10Input(buildfile)
    hel=load_dssv08(build["input_root"]);ff=load_dss14(build["input_root"])
    parameters=kpsy_published_parameters()
    xpoints=[.0065,.01,.03,.1,.2,.35,.6,.8]
    zpoints=[.05,.1,.2,.3,.5,.7,.9]
    refinements=[]
    try
        for n in policy["physical_grids"]
            grids=initial_source_grids(pdf,hel,ff;background_nodes=n,x_probes=xpoints,z_probes=zpoints)
            x,z=grids.x,grids.z
            inputs=prepare_initial_conditions(pdf,hel,ff,x,z)
            initial=initial_distributions(inputs,parameters)
            hop=Cop=nothing
            preparation=@elapsed begin
                hop=collinear_lo_operator(x,10.,2.4,AsymptoticTwoLoopCoupling())
                Cop=collinear_lo_operator(z,10.,2.4,AsymptoticTwoLoopCoupling();quantity=:KPSY_Hhat3_homogeneous)
            end
            evolved=evolve_initial_conditions(inputs,parameters,hop,Cop)
            hvalues=[interpolate(x,evolved.h1[row,:],u) for row in (7,6),u in xpoints]
            Hvalues=[interpolate(z,evolved.Hhat3[row,:],u) for row in (7,6,8),u in zpoints]
            excess=maximum(abs(initial.h1[pid==2 ? 7 : 6,j])-inputs.soffer[pid==2 ? 1 : 2,j]
                for pid in (2,1),j in eachindex(x))
            check("$n source Soffer envelope",max(0.,excess),policy["positivity_max_relative_excess"])
            check("$n exact evolved endpoints",max(maximum(abs,evolved.h1[:,end]),maximum(abs,evolved.Hhat3[:,end])),policy["endpoint_max_abs"])
            finite=all(isfinite,evolved.h1) && all(isfinite,evolved.Hhat3)
            check("$n finite initial/evolved physical-input distributions",finite ? 0. : 1.,0.)
            affine_errors=Float64[]
            for grid in (x,z)
                u=grid[1:end-1]
                reference=(4/3).*(2u.*log.(u).+(1.5.+2log1p.(-u)).*(1 .-u))
                push!(affine_errors,maximum(abs.((collinear_lo_generator(grid)*(1 .-grid))[1:end-1].-reference)))
            end
            check("$n aligned-grid analytic affine identity",maximum(affine_errors),policy["affine_generator_max_abs"])
            push!(refinements,Dict("background_nodes_each"=>n,"x_nodes"=>length(x),"z_nodes"=>length(z),
                "grid_definition"=>grids.definition,"Q2_GeV2"=>10.,"h1_values"=>[collect(r) for r in eachrow(hvalues)],
                "Hhat3_values"=>[collect(r) for r in eachrow(Hvalues)],
                "dense_operator_bytes"=>sizeof(hop.matrix)+sizeof(Cop.matrix),
                "operator_preparation_seconds_informational"=>preparation,
                "input_provenance"=>inputs.provenance))
            # Release the dense reference before constructing the next grid.
            hop=Cop=nothing
            GC.gc()
        end
        for field in ("h1_values","Hhat3_values")
            a=reduce(vcat,last(refinements)[field]);b=reduce(vcat,refinements[end-1][field])
            check("$field fine-grid refinement",scaled(a,b),policy["evolved_refinement_max_scaled"])
        end
        x=condition_grid(81,1e-5,.02);z=condition_grid(81,.05,.2)
        inputs=prepare_initial_conditions(pdf,hel,ff,x,z)
        plus=initial_distributions(inputs,parameters;charge=1)
        minus=initial_distributions(inputs,parameters;charge=-1)
        neutral=initial_distributions(inputs,parameters;charge=0)
        check("charged/neutral initial isospin relation",maximum(abs.(neutral.Hhat3-(plus.Hhat3+minus.Hhat3)/2)),policy["formula_max_scaled"])
        jacobian=ForwardDiff.jacobian(t -> vec(initial_distributions(inputs,KPSYParameters(t)).Hhat3),collect(parameters.values))
        check("pinned-input beta=0 initial Jacobian is finite",all(isfinite,jacobian) ? 0. : 1.,0.)
        check("source/build hashes stable during evolved input pilot",source_hashes()==sources &&
            bytes2hex(sha256(read(buildfile)))==build_hash && CollinsFit.verify_parton_build(buildfile)==build ? 0. : 1.,0.)
        Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
            "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
            "build_sha256"=>build_hash,"build"=>build,"contract_id"=>"KPSY-INITIAL-001",
            "printed_parameters"=>collect(parameters.values),"x_probes"=>xpoints,"z_probes"=>zpoints,
            "refinements"=>refinements,"julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),
            "dense_many_scale_estimate"=>Dict("assumed_radial_nodes"=>2336,
                "bytes_if_both_dense_operators_cached_per_node"=>last(refinements)["dense_operator_bytes"]*2336,
                "qualification"=>"resource estimate only; no speedup, fitting throughput or production performance claim"),
            "input_scope"=>"strict initial-scale pinned inputs, evolved rounded published parameters; no experiment or fit",
            "limitations"=>JSON.parsefile(CONTRACT)["limitations"])
    finally
        close(pdf)
    end
end

buildfile=isempty(ARGS) ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[1]
output=length(ARGS)<2 ? joinpath(ROOT,"runs/initial-inputs/validation.json") : ARGS[2]
try
    report=initial_input_checks(buildfile);writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
