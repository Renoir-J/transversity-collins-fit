# Physical-input backend pilot. Every failed accuracy case remains explicit.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, Statistics
include("CollinearReference.jl")
using .CollinearReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/hoppet-input-pilot.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/initial-contract.json"),joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"config/parton-contract.json"),joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"numerics/bin/CollinearReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function input_checks(hoppet_build,parton_build)
    sources=source_hashes();native=CollinsFit.verified_hoppet_build(hoppet_build)
    parton=CollinsFit.verify_parton_build(parton_build)
    build_hashes=Dict(p=>bytes2hex(sha256(read(p))) for p in (hoppet_build,parton_build))
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    q2=Float64.(contract["Q2_GeV2"]);coupling=AsymptoticTwoLoopCoupling()
    xp,zp=Float64.(contract["x_probes"]),Float64.(contract["z_probes"])
    parameters=kpsy_published_parameters();theta=collect(parameters.values)
    pdf=CT10Input(parton_build);hel=load_dssv08(parton["input_root"]);ff=load_dss14(parton["input_root"])
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    scaled(a,b)=maximum(abs.(a-b)./max.(abs.(b),policy["relative_floor"]))
    refinements=Dict{Symbol,Any}(:transversity=>[],:KPSY_Hhat3_homogeneous=>[])
    timings=[];dense=[];dense_previous=nothing
    try
        for quantity in (:transversity,:KPSY_Hhat3_homogeneous)
            ispdf=quantity==:transversity
            probes=ispdf ? xp : zp
            minimum=contract[ispdf ? "x_min" : "z_min"]
            rows=ispdf ? [7,6] : [7,6,8]
            for dy in contract["dy_refinements"]
                plan=nothing
                preparation=@elapsed plan=HoppetPlan(hoppet_build,q2,2.4,coupling;quantity=quantity,
                    xmin=minimum,dy=dy,du=contract["du"])
                try
                    inputs=ispdf ? prepare_initial_conditions(pdf,hel,ff,plan.grid,[.05,.5,1.]) :
                                   prepare_initial_conditions(pdf,hel,ff,[1e-5,.1,1.],plan.grid)
                    function initial(t)
                        family=initial_distributions(inputs,KPSYParameters(t))
                        ispdf ? permutedims(family.h1[rows,:]) : permutedims(family.Hhat3[rows,:])
                    end
                    initial_values=initial(theta)
                    evolved=nothing;warm=@elapsed evolved=apply_operator(plan,initial_values)
                    values=hoppet_evaluate(plan,evolved,probes)
                    check("$quantity dy=$dy finite",all(isfinite,evolved) ? 0. : 1.,0.)
                    check("$quantity dy=$dy exact endpoints",maximum(abs,evolved[end,:,:]),policy["endpoint_max_abs"])
                    check("$quantity dy=$dy Q0 identity",maximum(abs.(evolved[:,:,2]-initial_values)),0.)
                    push!(refinements[quantity],Dict("dy"=>dy,"values"=>values,"unique_nodes"=>length(plan.grid),
                        "native_nodes"=>plan.native_nodes,"coefficient_payload_bytes"=>plan.coefficient_payload_bytes,
                        "preparation_seconds_informational"=>preparation,"apply_seconds_informational"=>warm,
                        "xmin"=>first(plan.grid),"initial_provenance"=>inputs.provenance))
                    if dy==last(contract["dy_refinements"])
                        prediction(t)=vec(hoppet_evaluate(plan,apply_operator(plan,initial(t)),probes))
                        direction=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,0.]
                        J=ForwardDiff.jacobian(prediction,theta)
                        check("$quantity finite physical Jacobian at published boundaries",all(isfinite,J) ? 0. : 1.,0.)
                        for step in policy["finite_difference_steps"]
                            fd=(-3prediction(theta)+4prediction(theta+step*direction)-prediction(theta+2step*direction))/(2step)
                            check("$quantity inward JVP step $step",scaled(fd,J*direction),policy["derivative_max_scaled"])
                        end
                        loss(t)=sum(abs2,prediction(t))
                        gradient(t)=ForwardDiff.gradient(loss,t)
                        hv=ForwardDiff.jacobian(s -> gradient(theta+s[1]*direction),[0.])[:,1]
                        for step in policy["finite_difference_steps"]
                            fd=(-3gradient(theta)+4gradient(theta+step*direction)-gradient(theta+2step*direction))/(2step)
                            check("$quantity inward HVP step $step",scaled(hv,fd),policy["derivative_max_scaled"])
                        end
                    end
                finally
                    close(plan)
                end
                println(JSON.json(Dict("completed_quantity"=>string(quantity),"dy"=>dy)));flush(stdout)
            end
            fine=refinements[quantity][end]["values"]
            previous=refinements[quantity][end-1]["values"]
            for k in eachindex(q2)
                check("$quantity refinement Q2=$(q2[k])",scaled(fine[:,:,k],previous[:,:,k]),policy["value_refinement_max_scaled"])
            end
        end
        # A fresh dense route is reconstructed from strict source/probe-aligned
        # initial grids. Its own refinement is checked before treating it as a
        # comparison; no old receipt or unreviewed legacy object is loaded.
        for n in contract["dense_background_refinements"]
            grids=initial_source_grids(pdf,hel,ff;background_nodes=n,x_probes=xp,z_probes=zp)
            inputs=prepare_initial_conditions(pdf,hel,ff,grids.x,grids.z)
            samples=Dict(:transversity=>zeros(length(xp),2,length(q2)),
                :KPSY_Hhat3_homogeneous=>zeros(length(zp),3,length(q2)))
            for k in eachindex(q2)
                op=Cop=nothing
                cold=@elapsed begin
                    op=collinear_lo_operator(grids.x,q2[k],2.4,coupling)
                    Cop=collinear_lo_operator(grids.z,q2[k],2.4,coupling;quantity=:KPSY_Hhat3_homogeneous)
                end
                evolved=evolve_initial_conditions(inputs,parameters,op,Cop)
                samples[:transversity][:,:,k]=[interpolate_grid(grids.x,evolved.h1[row,:],v) for v in xp,row in (7,6)]
                samples[:KPSY_Hhat3_homogeneous][:,:,k]=[interpolate_grid(grids.z,evolved.Hhat3[row,:],v) for v in zp,row in (7,6,8)]
                push!(timings,Dict("background_nodes"=>n,"Q2_GeV2"=>q2[k],
                    "dense_operator_bytes"=>sizeof(op.matrix)+sizeof(Cop.matrix),"preparation_seconds_informational"=>cold))
                op=Cop=nothing;GC.gc()
            end
            push!(dense,Dict("background_nodes"=>n,"x_nodes"=>length(grids.x),"z_nodes"=>length(grids.z),
                "h1_values"=>samples[:transversity],"Hhat3_values"=>samples[:KPSY_Hhat3_homogeneous]))
            if dense_previous!==nothing
                for quantity in keys(samples),k in eachindex(q2)
                    check("$quantity dense-reference refinement Q2=$(q2[k])",
                        scaled(samples[quantity][:,:,k],dense_previous[quantity][:,:,k]),policy["value_refinement_max_scaled"])
                    check("$quantity retained/dense comparison Q2=$(q2[k])",
                        scaled(refinements[quantity][end]["values"][:,:,k],samples[quantity][:,:,k]),policy["cross_backend_max_scaled"])
                end
            end
            dense_previous=samples
        end
        check("source/build hashes stable",source_hashes()==sources && all(bytes2hex(sha256(read(p)))==d for (p,d) in build_hashes) &&
            CollinsFit.verified_hoppet_build(hoppet_build)==native && CollinsFit.verify_parton_build(parton_build)==parton ? 0. : 1.,0.)
        Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
            "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
            "contract_id"=>contract["id"],"build_sha256"=>build_hashes,"native_build"=>native,"parton_build"=>parton,
            "printed_parameters"=>theta,"Q2_GeV2"=>q2,"x_probes"=>xp,"z_probes"=>zp,
            "refinements"=>Dict(string(k)=>v for (k,v) in refinements),"dense_refinements"=>dense,
            "dense_timings_informational"=>timings,"limitations"=>contract["limitations"],
            "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
    finally
        close(pdf)
    end
end

native=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
parton=length(ARGS)<2 ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(ROOT,"runs/hoppet-inputs/validation.json") : ARGS[3]
try
    report=input_checks(native,parton);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
