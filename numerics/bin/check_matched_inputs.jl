# Strict pinned initial inputs, matched evolved legs. Never a physical extraction.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/matched-rows-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/hoppet-point-contract.json"),joinpath(ROOT,"config/initial-contract.json"),
    joinpath(ROOT,"config/leptonic-contract.json"),joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"config/parton-contract.json"),joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

"Independent regular matching of piecewise-linear dense evolved values."
function dense_regular_rows(grid,queries,order)
    nodes,weights=CollinsFit.gauss_rule(order)
    rows=zeros(length(queries),length(grid))
    for (j,x) in enumerate(queries)
        first=searchsortedfirst(grid,x)
        grid[first]==x || error("dense reference query is not an exact grid knot")
        for i in first:length(grid)-1
            low,high=grid[i:i+1];half=(log(high)-log(low))/2;center=(log(high)+log(low))/2
            for k in eachindex(nodes)
                y=exp(center+half*nodes[k]);v=x/y
                L=v==1 ? -1. : log(v)/(1-v)
                w=half*weights[k]*(8/3)*v*L
                fraction=(y-low)/(high-low)
                rows[j,i]+=w*(1-fraction);rows[j,i+1]+=w*fraction
            end
        end
    end
    rows
end

function matched_input_checks(nativefile,partonfile)
    sources=source_hashes();native=CollinsFit.verified_hoppet_build(nativefile);parton=CollinsFit.verify_parton_build(partonfile)
    builds=Dict(p=>bytes2hex(sha256(read(p))) for p in (nativefile,partonfile))
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    q2=Float64.(policy["Q2_GeV2"]);q0=policy["Q02_GeV2"];coupling=AsymptoticTwoLoopCoupling()
    xp,zp=Float64.(policy["x_probes"]),Float64.(policy["z_probes"])
    p=kpsy_published_parameters();theta=collect(p.values)
    pdf=CT10Input(partonfile);hel=load_dssv08(parton["input_root"]);ff=load_dss14(parton["input_root"])
    cases=((:SIDIS,:transversity),(:SIDIS,:KPSY_Hhat3_homogeneous),(:SIA,:KPSY_Hhat3_homogeneous))
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    scaled(a,b)=maximum(abs.(a-b)./max.(abs.(b),policy["physical_relative_floor"]))
    refinements=Dict(case=>Any[] for case in cases);native_values=Dict{Tuple{Symbol,Symbol},Any}()
    dense_timings=[];dense_grids=[];reference_comparisons=[];q0_references=Dict{Int,Any}()
    try
        for (process,quantity) in cases
            ispdf=quantity==:transversity;probes=ispdf ? xp : zp;label="$process $quantity"
            for resolution in policy["physical_interval_refinements"]
                plan=nothing;cold=@elapsed plan=HoppetMatchingPlan(nativefile,vcat(probes,1.),q2,q0,coupling;
                    process=process,quantity=quantity,intervals=resolution)
                inputs=ispdf ? prepare_initial_conditions(pdf,hel,ff,plan.grid,[.05,.5,1.]) :
                    prepare_initial_conditions(pdf,hel,ff,[1e-5,.1,1.],plan.grid)
                function initial(t)
                    parameters=KPSYParameters(t)
                    ispdf ? [initial_transversity(parameters,pid,x,inputs.soffer[f,j]) for (j,x) in enumerate(inputs.x),(f,pid) in enumerate((2,1))] :
                        [initial_collins(parameters,pid,1,z,view(inputs.ff,:,j)) for (j,z) in enumerate(inputs.z),pid in (2,1,3)]
                end
                values=nothing;warm=@elapsed values=apply_operator(plan,initial(theta))
                check("$label $resolution finite",all(isfinite,values) ? 0. : 1.,0.)
                check("$label $resolution endpoint",maximum(abs,values[end,:,:]),0.)
                check("$label $resolution native row comparison",plan.row_native_max_scaled,policy["row_native_max_scaled"])
                push!(refinements[(process,quantity)],Dict("intervals"=>resolution,"values"=>values[1:end-1,:,:],
                    "unique_nodes"=>length(plan.grid),"coefficient_payload_bytes"=>plan.coefficient_payload_bytes,
                    "preparation_seconds_informational"=>cold,"apply_seconds_informational"=>warm,"initial_provenance"=>inputs.provenance))
                levels=refinements[(process,quantity)]
                error=length(levels)>1 ? maximum(scaled(levels[end]["values"][:,:,k],levels[end-1]["values"][:,:,k]) for k in eachindex(q2)) : nothing
                converged=resolution>=policy["minimum_physical_intervals"] && error!==nothing && error<=policy["physical_value_max_scaled"]
                levels[end]["previous_level_max_error"]=error
                levels[end]["refinement_accepted"]=converged
                finish=converged || resolution==last(policy["physical_interval_refinements"])
                if finish
                    native_values[(process,quantity)]=values[1:end-1,:,:]
                    prediction(t)=vec(apply_operator(plan,initial(t))[1:end-1,:,:])
                    direction=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,0.]
                    J=ForwardDiff.jacobian(prediction,theta)
                    loss(t)=sum(abs2,prediction(t));gradient(t)=ForwardDiff.gradient(loss,t)
                    hv=ForwardDiff.jacobian(s -> gradient(theta+s[1]*direction),[0.])[:,1]
                    check("$label finite physical boundary Jacobian/HVP",all(isfinite,J) && all(isfinite,hv) ? 0. : 1.,0.)
                    for step in policy["finite_difference_steps"]
                        fd,fdh=setprecision(policy["finite_difference_reference_bits"]) do
                            t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
                            (Float64.((-3prediction(t)+4prediction(t+h*d)-prediction(t+2h*d))/(2h)),
                             Float64.((-3gradient(t)+4gradient(t+h*d)-gradient(t+2h*d))/(2h)))
                        end
                        check("$label inward JVP step $step",scaled(fd,J*direction),policy["derivative_max_scaled"])
                        check("$label inward HVP step $step",scaled(fdh,hv),policy["derivative_max_scaled"])
                    end
                end
                println(JSON.json(Dict("completed_process"=>string(process),"quantity"=>string(quantity),"intervals"=>resolution)));flush(stdout)
                finish && break
            end
            fine,coarse=refinements[(process,quantity)][end]["values"],refinements[(process,quantity)][end-1]["values"]
            for k in eachindex(q2)
                check("$label native refinement Q2=$(q2[k])",scaled(fine[:,:,k],coarse[:,:,k]),policy["physical_value_max_scaled"])
            end
        end
        # Direct strict input matching at Q0, split at DSS14 source knots.
        direct_inputs=prepare_initial_conditions(pdf,hel,ff,vcat(xp,1.),[.05,.5,1.])
        h0=[initial_transversity(p,pid,x,direct_inputs.soffer[f,j]) for (j,x) in enumerate(xp),(f,pid) in enumerate((2,1))]
        alpha0=alpha_s(coupling,q0)
        for order in (16,32,64)
            direct=Dict(case=>zeros(case[2]==:transversity ? (length(xp),2) : (length(zp),3)) for case in cases)
            direct[(:SIDIS,:transversity)]=(1-(8/3)*alpha0/pi).*h0
            for (j,z) in enumerate(zp)
                rule=matching_rule(:SIDIS,:collins,z,alpha0;order=order,knots=ff.x)
                bases=[pion_ff(ff,[2,1,3,-3],1,y,q0) for y in rule.coordinates]
                values=[initial_collins(p,pid,1,y,bases[k]) for (k,y) in enumerate(rule.coordinates),pid in (2,1,3)]
                for process in (:SIDIS,:SIA)
                    other=process==:SIDIS ? rule : matching_rule(process,:collins,z,alpha0;order=order,knots=ff.x)
                    other.coordinates==rule.coordinates || error("reference matching coordinates differ")
                    direct[(process,:KPSY_Hhat3_homogeneous)][j,:]=vec(transpose(other.quark_weights)*values)
                end
            end
            q0_references[order]=direct
        end
        k0=findfirst(==(q0),q2)
        for case in cases
            check("$case Q0 direct quadrature refinement",scaled(q0_references[64][case],q0_references[32][case]),policy["physical_value_max_scaled"])
            check("$case Q0 independent strict-input matching",scaled(native_values[case][:,:,k0],q0_references[64][case]),policy["physical_value_max_scaled"])
        end
        # A separately reconstructed dense evolution and independent cell-wise
        # matching quadrature. Cache only small probe outputs, not dense matrices.
        function dense_values(quantity,n,ks)
            ispdf=quantity==:transversity;probes=ispdf ? xp : zp;rows=ispdf ? [7,6] : [7,6,8]
            grids=initial_source_grids(pdf,hel,ff;background_nodes=n,x_probes=xp,z_probes=zp)
            inputs=prepare_initial_conditions(pdf,hel,ff,grids.x,grids.z)
            initial=initial_distributions(inputs,p);grid=ispdf ? grids.x : grids.z
            base=permutedims((ispdf ? initial.h1 : initial.Hhat3)[rows,:])
            indices=searchsortedfirst.(Ref(grid),probes);grid[indices]==probes || error("dense reference probes not aligned")
            W16=ispdf ? nothing : dense_regular_rows(grid,probes,16)
            W32=ispdf ? nothing : dense_regular_rows(grid,probes,32)
            processes=ispdf ? (:SIDIS,) : (:SIDIS,:SIA)
            values=Dict(process=>zeros(length(probes),length(rows),length(q2)) for process in processes)
            push!(dense_grids,Dict("quantity"=>string(quantity),"background_nodes"=>n,"nodes"=>length(grid)))
            for k in ks
                op=nothing;evolved=nothing
                cold=@elapsed if q2[k]==q0
                    evolved=base
                else
                    op=collinear_lo_operator(grid,q2[k],q0,coupling;quantity=quantity)
                    evolved=apply_operator(op,base)
                end
                regular=ispdf ? zeros(length(probes),length(rows)) : W32*evolved
                if !ispdf
                    check("$quantity dense n=$n quadrature Q2=$(q2[k])",scaled(W16*evolved,regular),policy["physical_value_max_scaled"])
                end
                a=alpha_s(coupling,q2[k])/pi
                for process in processes
                    delta=process==:SIDIS ? -8/3 : (pi^2-8)/3
                    values[process][:,:,k]=(1+a*delta).*evolved[indices,:]+a.*regular
                end
                push!(dense_timings,Dict("quantity"=>string(quantity),"background_nodes"=>n,"Q2_GeV2"=>q2[k],
                    "operator_bytes"=>isnothing(op) ? 0 : sizeof(op.matrix),"preparation_seconds_informational"=>cold))
                op=evolved=nothing;GC.gc()
            end
            println(JSON.json(Dict("completed_dense"=>string(quantity),"background_nodes"=>n,"scales"=>collect(ks))));flush(stdout)
            values
        end
        dense=Dict{Tuple{Symbol,Int},Any}();budget=Int.(policy["dense_reference_background_budget"])
        for n in budget[1:2],quantity in (:transversity,:KPSY_Hhat3_homogeneous)
            dense[(quantity,n)]=dense_values(quantity,n,eachindex(q2))
        end
        extra=Dict{Tuple{Symbol,Int,Int},Any}()
        for case in cases,k in eachindex(q2)
            process,quantity=case
            coarse=dense[(quantity,budget[1])][process][:,:,k];fine=dense[(quantity,budget[2])][process][:,:,k]
            error=scaled(fine,coarse);ladder=Any[Dict("background_nodes"=>budget[2],"refinement_error"=>error)]
            if error>policy["physical_value_max_scaled"]
                for n in budget[3:end]
                    key=(quantity,n,k)
                    finer=get!(extra,key) do
                        dense_values(quantity,n,[k])
                    end
                    next=finer[process][:,:,k];error=scaled(next,fine);coarse=fine;fine=next
                    push!(ladder,Dict("background_nodes"=>n,"refinement_error"=>error))
                    error<=policy["physical_value_max_scaled"] && break
                end
            end
            cross=scaled(native_values[case][:,:,k],fine)
            check("$case dense own refinement Q2=$(q2[k])",error,policy["physical_value_max_scaled"])
            check("$case matched rows/dense Q2=$(q2[k])",cross,policy["physical_value_max_scaled"])
            push!(reference_comparisons,Dict("process"=>string(process),"quantity"=>string(quantity),"Q2_GeV2"=>q2[k],
                "refinement_ladder"=>ladder,"row_dense_error"=>cross))
        end
        check("source/build stability",source_hashes()==sources && all(bytes2hex(sha256(read(p)))==d for (p,d) in builds) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(partonfile)==parton ? 0. : 1.,0.)
        Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
            "production_ready"=>false,"physical_extraction"=>false,"contract_id"=>contract["id"],"sources"=>sources,
            "checks"=>checks,"build_sha256"=>builds,"native_build"=>native,"parton_build"=>parton,
            "finite_difference_reference_bits"=>policy["finite_difference_reference_bits"],"printed_parameters"=>theta,
            "Q2_GeV2"=>q2,"refinements"=>Dict(string(k)=>v for (k,v) in refinements),
            "Q0_direct_refinements"=>Dict(string(order)=>Dict(string(k)=>v for (k,v) in values) for (order,values) in q0_references),
            "dense_grids"=>dense_grids,"dense_timings_informational"=>dense_timings,"reference_comparisons"=>reference_comparisons,
            "limitations"=>contract["limits"],"julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
    finally
        close(pdf)
    end
end

native=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
parton=length(ARGS)<2 ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(ROOT,"runs/matched-inputs/validation.json") : ARGS[3]
try
    report=matched_input_checks(native,parton);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"error"=>sprint(showerror,exception)))
    rethrow()
end
