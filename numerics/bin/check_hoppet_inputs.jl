# Physical-input backend pilot. Every failed accuracy case remains explicit.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, Statistics
include("CollinearReference.jl")
using .CollinearReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/hoppet-input-pilot.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-point-contract.json"),joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/initial-contract.json"),joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"config/parton-contract.json"),joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"numerics/bin/CollinearReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function input_checks(hoppet_build,parton_build;points=false)
    sources=source_hashes();native=CollinsFit.verified_hoppet_build(hoppet_build)
    parton=CollinsFit.verify_parton_build(parton_build)
    build_hashes=Dict(p=>bytes2hex(sha256(read(p))) for p in (hoppet_build,parton_build))
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    point_contract=points ? JSON.parsefile(joinpath(ROOT,"config/hoppet-point-contract.json")) : nothing
    q2=Float64.(contract["Q2_GeV2"]);coupling=AsymptoticTwoLoopCoupling()
    xp,zp=Float64.(contract["x_probes"]),Float64.(contract["z_probes"])
    parameters=kpsy_published_parameters();theta=collect(parameters.values)
    pdf=CT10Input(parton_build);hel=load_dssv08(parton["input_root"]);ff=load_dss14(parton["input_root"])
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    scaled(a,b)=maximum(abs.(a-b)./max.(abs.(b),policy["relative_floor"]))
    refinements=Dict{Symbol,Any}(:transversity=>[],:KPSY_Hhat3_homogeneous=>[])
    timings=[];dense=[];dense_previous=nothing;dense_coarse=nothing
    try
        for quantity in (:transversity,:KPSY_Hhat3_homogeneous)
            ispdf=quantity==:transversity
            probes=ispdf ? xp : zp
            minimum=contract[ispdf ? "x_min" : "z_min"]
            rows=ispdf ? [7,6] : [7,6,8]
            resolutions=points ? point_contract["validation_policy"]["interval_refinements"] : contract["dy_refinements"]
            for resolution in resolutions
                plan=nothing
                preparation=@elapsed plan=points ? HoppetPointPlan(hoppet_build,vcat(probes,1.),q2,2.4,coupling;
                    quantity=quantity,intervals=resolution,du=contract["du"]) :
                    HoppetPlan(hoppet_build,q2,2.4,coupling;quantity=quantity,xmin=minimum,dy=resolution,du=contract["du"])
                try
                    inputs=ispdf ? prepare_initial_conditions(pdf,hel,ff,plan.grid,[.05,.5,1.]) :
                                   prepare_initial_conditions(pdf,hel,ff,[1e-5,.1,1.],plan.grid)
                    function initial(t)
                        p=KPSYParameters(t)
                        if ispdf
                            [initial_transversity(p,pid,x,inputs.soffer[f,j]) for (j,x) in enumerate(inputs.x),(f,pid) in enumerate((2,1))]
                        else
                            [initial_collins(p,pid,1,z,view(inputs.ff,:,j)) for (j,z) in enumerate(inputs.z),pid in (2,1,3)]
                        end
                    end
                    initial_values=initial(theta)
                    evolved=nothing;warm=@elapsed evolved=apply_operator(plan,initial_values)
                    values=points ? evolved[1:end-1,:,:] : hoppet_evaluate(plan,evolved,probes)
                    target_initial=points ? initial_values[searchsortedfirst.(Ref(plan.grid),plan.queries),:] : initial_values
                    check("$quantity resolution=$resolution finite",all(isfinite,evolved) ? 0. : 1.,0.)
                    check("$quantity resolution=$resolution exact endpoints",maximum(abs,evolved[end,:,:]),policy["endpoint_max_abs"])
                    check("$quantity resolution=$resolution Q0 identity",maximum(abs.(evolved[:,:,2]-target_initial)),0.)
                    push!(refinements[quantity],Dict("resolution"=>resolution,"kind"=>points ? "query_bound_rows" : "global_grid",
                        "values"=>values,"unique_nodes"=>length(plan.grid),
                        "native_nodes"=>points ? length(probes)*(resolution+1) : plan.native_nodes,"coefficient_payload_bytes"=>plan.coefficient_payload_bytes,
                        "preparation_seconds_informational"=>preparation,"apply_seconds_informational"=>warm,
                        "xmin"=>first(plan.grid),"initial_provenance"=>inputs.provenance))
                    if resolution==last(resolutions)
                        prediction(t)=vec(points ? apply_operator(plan,initial(t))[1:end-1,:,:] : hoppet_evaluate(plan,apply_operator(plan,initial(t)),probes))
                        direction=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,0.]
                        J=ForwardDiff.jacobian(prediction,theta)
                        check("$quantity finite physical Jacobian at published boundaries",all(isfinite,J) ? 0. : 1.,0.)
                        for step in policy["finite_difference_steps"]
                            fd=if points
                                setprecision(point_contract["validation_policy"]["finite_difference_reference_bits"]) do
                                    t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
                                    Float64.((-3prediction(t)+4prediction(t+h*d)-prediction(t+2h*d))/(2h))
                                end
                            else
                                (-3prediction(theta)+4prediction(theta+step*direction)-prediction(theta+2step*direction))/(2step)
                            end
                            check("$quantity inward JVP step $step",scaled(fd,J*direction),policy["derivative_max_scaled"])
                        end
                        loss(t)=sum(abs2,prediction(t))
                        gradient(t)=ForwardDiff.gradient(loss,t)
                        hv=ForwardDiff.jacobian(s -> gradient(theta+s[1]*direction),[0.])[:,1]
                        for step in policy["finite_difference_steps"]
                            fd=if points
                                setprecision(point_contract["validation_policy"]["finite_difference_reference_bits"]) do
                                    t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
                                    Float64.((-3gradient(t)+4gradient(t+h*d)-gradient(t+2h*d))/(2h))
                                end
                            else
                                (-3gradient(theta)+4gradient(theta+step*direction)-gradient(theta+2step*direction))/(2step)
                            end
                            check("$quantity inward HVP step $step",scaled(hv,fd),policy["derivative_max_scaled"])
                        end
                    end
                finally
                    !points && close(plan)
                end
                println(JSON.json(Dict("completed_quantity"=>string(quantity),"resolution"=>resolution)));flush(stdout)
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
                op=Cop=nothing;evolved=nothing;bytes=0
                cold=@elapsed begin
                    if q2[k]==2.4
                        evolved=initial_distributions(inputs,parameters)
                    else
                        op=collinear_lo_operator(grids.x,q2[k],2.4,coupling)
                        Cop=collinear_lo_operator(grids.z,q2[k],2.4,coupling;quantity=:KPSY_Hhat3_homogeneous)
                        evolved=evolve_initial_conditions(inputs,parameters,op,Cop)
                        bytes=sizeof(op.matrix)+sizeof(Cop.matrix)
                    end
                end
                samples[:transversity][:,:,k]=[interpolate_grid(grids.x,evolved.h1[row,:],v) for v in xp,row in (7,6)]
                samples[:KPSY_Hhat3_homogeneous][:,:,k]=[interpolate_grid(grids.z,evolved.Hhat3[row,:],v) for v in zp,row in (7,6,8)]
                push!(timings,Dict("background_nodes"=>n,"Q2_GeV2"=>q2[k],
                    "dense_operator_bytes"=>bytes,"preparation_seconds_informational"=>cold))
                op=Cop=nothing;GC.gc()
            end
            push!(dense,Dict("background_nodes"=>n,"x_nodes"=>length(grids.x),"z_nodes"=>length(grids.z),
                "h1_values"=>samples[:transversity],"Hhat3_values"=>samples[:KPSY_Hhat3_homogeneous]))
            dense_coarse=dense_previous
            dense_previous=samples
        end
        reference_comparisons=[];extra_inputs=nothing;extra_grids=nothing
        for quantity in keys(dense_previous),k in eachindex(q2)
            coarse=dense_coarse[quantity][:,:,k]
            fine=dense_previous[quantity][:,:,k]
            refinement=scaled(fine,coarse);last_n=last(contract["dense_background_refinements"])
            ladder=Any[Dict("nodes"=>last_n,"refinement_error"=>refinement)]
            if points && refinement>policy["value_refinement_max_scaled"]
                for n in point_contract["validation_policy"]["dense_reference_background_budget"]
                    n<=last_n && continue
                    if extra_inputs===nothing
                        extra_grids=initial_source_grids(pdf,hel,ff;background_nodes=n,x_probes=xp,z_probes=zp)
                        extra_inputs=prepare_initial_conditions(pdf,hel,ff,extra_grids.x,extra_grids.z)
                    end
                    initial=initial_distributions(extra_inputs,parameters)
                    ispdf=quantity==:transversity
                    grid=ispdf ? extra_grids.x : extra_grids.z
                    rows=ispdf ? [7,6] : [7,6,8]
                    probes=ispdf ? xp : zp
                    op=nothing;elapsed=@elapsed op=collinear_lo_operator(grid,q2[k],2.4,coupling;quantity=quantity)
                    evolved=apply_operator(op,permutedims((ispdf ? initial.h1 : initial.Hhat3)[rows,:]))
                    finer=[interpolate_grid(grid,evolved[:,j],v) for v in probes,j in eachindex(rows)]
                    refinement=scaled(finer,fine);coarse=fine;fine=finer;last_n=n
                    push!(ladder,Dict("nodes"=>n,"refinement_error"=>refinement,
                        "dense_operator_bytes"=>sizeof(op.matrix),"preparation_seconds_informational"=>elapsed))
                    op=nothing;GC.gc()
                    refinement<=policy["value_refinement_max_scaled"] && break
                end
            end
            cross=scaled(refinements[quantity][end]["values"][:,:,k],fine)
            check("$quantity dense-reference refinement Q2=$(q2[k])",refinement,policy["value_refinement_max_scaled"])
            check("$quantity retained/dense comparison Q2=$(q2[k])",cross,policy["cross_backend_max_scaled"])
            push!(reference_comparisons,Dict("quantity"=>string(quantity),"Q2_GeV2"=>q2[k],
                "ladder"=>ladder,"final_cross_error"=>cross,"reference_values"=>fine))
        end
        check("source/build hashes stable",source_hashes()==sources && all(bytes2hex(sha256(read(p)))==d for (p,d) in build_hashes) &&
            CollinsFit.verified_hoppet_build(hoppet_build)==native && CollinsFit.verify_parton_build(parton_build)==parton ? 0. : 1.,0.)
        Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
            "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
            "contract_id"=>points ? point_contract["id"] : contract["id"],
            "backend"=>points ? "query_bound_rows" : "global_grid",
            "finite_difference_reference_bits"=>points ? point_contract["validation_policy"]["finite_difference_reference_bits"] : 53,
            "build_sha256"=>build_hashes,"native_build"=>native,"parton_build"=>parton,
            "printed_parameters"=>theta,"Q2_GeV2"=>q2,"x_probes"=>xp,"z_probes"=>zp,
            "refinements"=>Dict(string(k)=>v for (k,v) in refinements),"dense_refinements"=>dense,
            "dense_timings_informational"=>timings,"reference_comparisons"=>reference_comparisons,
            "resource_pilot"=>Dict("status"=>"SEPARATE_JOB","runner"=>"numerics/bin/check_hoppet_points.jl",
                "scope"=>"analytic_resource_pilot","production_ready"=>false),
            "limitations"=>points ? vcat(contract["limitations"],point_contract["limits"]) : contract["limitations"],
            "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
    finally
        close(pdf)
    end
end

native=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
parton=length(ARGS)<2 ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(ROOT,"runs/hoppet-inputs/validation.json") : ARGS[3]
try
    mode=length(ARGS)<4 ? "global_grid" : ARGS[4]
    mode in ("global_grid","query_rows") || error("choose global_grid or query_rows")
    report=input_checks(native,parton;points=mode=="query_rows");writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
