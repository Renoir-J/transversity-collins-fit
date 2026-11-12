# Retained construction validation. Analytic inputs and resource measurements are
# not an experimental fit or a production G2PERF certificate.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, Statistics
include("MatchedReference.jl")
using .MatchedReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/hoppet-batch-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/hoppet-point-contract.json"),joinpath(ROOT,"config/matched-rows-contract.json"),
    joinpath(ROOT,"config/leptonic-contract.json"),joinpath(ROOT,"config/parton-contract.json"),
    joinpath(ROOT,"config/parton-sources.json"),joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"src/collinsflow/hoppet.py"),
    joinpath(ROOT,"scripts/prepare_hoppet.py"),joinpath(ROOT,"numerics/bin/MatchedReference.jl"),
    joinpath(ROOT,"numerics/bin/CollinearReference.jl"),joinpath(ROOT,"numerics/Project.toml"),
    joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)
scaled(a,b)=maximum(abs.(a.-b)./max.(abs.(b),.01maximum(abs,b),1e-15))
rejected(f)=try f();false catch;true end
initial_columns(x)=hcat(x.*(1 .-x).^3,x.^2 .* (1 .-x).^2,x.^3 .* (1 .-x).^3)

function batch_checks(buildfile,inputfile)
    sources=source_hashes();build=CollinsFit.verified_hoppet_build(buildfile)
    inputs=CollinsFit.verify_parton_build(inputfile)
    hashes=(bytes2hex(sha256(read(buildfile))),bytes2hex(sha256(read(inputfile))))
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error)&&error<=tolerance))
    coupling=AsymptoticTwoLoopCoupling();mixed=Float64.(policy["mixed_Q2_GeV2"])
    for layout in (:single,:nested)
        saved=nothing
        for construction in (:independent,:batched_probes)
            plan=HoppetPlan(buildfile,mixed,2.4,coupling;xmin=.0065,dy=.1,layout=layout,construction=construction)
            try
                h=initial_columns(plan.grid);values=apply_operator(plan,h)
                rows=CollinsFit.hoppet_rows(plan,1)
                g=reshape(sin.(1:length(values)),size(values))
                check("$layout $construction native adjoint",abs(dot(values,g)-dot(h,operator_adjoint(plan,g)))/max(1.,abs(dot(values,g))),policy["adjoint_max_scaled"])
                check("$layout $construction identity",maximum(abs.(values[:,:,3]-h)),policy["identity_max_abs"])
                check("$layout $construction repeated identity",maximum(abs.(values[:,:,3]-values[:,:,7])),0.)
                check("$layout $construction repeated negative scale",maximum(abs.(values[:,:,2]-values[:,:,5])),0.)
                check("$layout $construction constrained endpoint",maximum(abs,values[end,:,:]),0.)
                if saved===nothing
                    saved=(values=values,rows=rows,grid=copy(plan.grid))
                else
                    check("$layout grid unchanged",plan.grid==saved.grid ? 0. : 1.,0.)
                    check("$layout mixed-scale values versus independent",scaled(values,saved.values),policy["value_max_scaled"])
                    check("$layout raw rows versus independent",scaled(rows,saved.rows),policy["operator_row_max_scaled"])
                end
            finally
                close(plan)
            end
        end
    end
    queries=Float64.(policy["queries"])
    snapshots=[]
    for process in (:SIDIS,:SIA)
        plans=[HoppetMatchingPlan(buildfile,vcat(queries,1.,queries[2]),mixed,2.4,coupling;
            process=process,intervals=1024,construction=construction) for construction in (:independent,:batched_probes)]
        values=[apply_operator(p,initial_columns(p.grid)) for p in plans]
        check("$process matched mixed-scale values versus independent",scaled(values[2],values[1]),policy["value_max_scaled"])
        p=plans[2];h=initial_columns(p.grid)
        reference(a,b)=[x==1 ? 0. : matched_polynomial(x,t,alpha,process,:KPSY_Hhat3_homogeneous;a=a,b=b)
            for x in p.queries,(t,alpha) in zip(p.tau,p.alpha)]
        R=hcat([vec(reference(a,b)) for (a,b) in ((1,3),(2,2),(3,3))]...)
        for (j,(a,b)) in enumerate(((1,3),(2,2),(3,3)))
            check("$process independent Mellin column $j",scaled(values[2][:,j,:],reference(a,b)),policy["value_max_scaled"])
        end
        coeff(t)=[exp(t[1]),sin(t[2]),t[2]*t[3]]
        prediction(t)=vec(apply_operator(p,h*coeff(t)))
        theta=[.2,-.3,.1];direction=[.3,-.2,.4]
        J=R*[exp(theta[1]) 0. 0.;0. cos(theta[2]) 0.;0. theta[3] theta[2]]
        check("$process batched tangent Jacobian versus Mellin",scaled(ForwardDiff.jacobian(prediction,theta),J),policy["derivative_max_scaled"])
        loss(t)=sum(abs2,prediction(t))
        H=2J'*J;v=R*coeff(theta)
        H[1,1]+=2dot(v,R[:,1]*exp(theta[1]));H[2,2]-=2dot(v,R[:,2]*sin(theta[2]))
        H[2,3]+=2dot(v,R[:,3]);H[3,2]+=2dot(v,R[:,3])
        check("$process nested batched Hessian versus Mellin",scaled(ForwardDiff.hessian(loss,theta),H),policy["derivative_max_scaled"])
        check("$process nested batched HVP versus Mellin",scaled(ForwardDiff.jacobian(s->ForwardDiff.gradient(loss,theta+s[1]*direction),[0.])[:,1],H*direction),policy["derivative_max_scaled"])
        big=BigFloat.(theta)
        check("$process Float64 tangent batch versus generic BigFloat",scaled(ForwardDiff.hessian(loss,theta),Float64.(ForwardDiff.hessian(loss,big))),policy["derivative_max_scaled"])
        check("$process repeated query",maximum(abs.(values[2][2,:,:]-values[2][end,:,:])),0.)
        push!(snapshots,(plan=p,initial=h,values=values[2]))
    end
    for snapshot in snapshots
        check("$(snapshot.plan.process) detached rows after later contexts",maximum(abs.(apply_operator(snapshot.plan,snapshot.initial)-snapshot.values)),0.)
    end
    check("unknown construction rejected",rejected(()->HoppetPlan(buildfile,[10.],2.4,coupling;construction=:interpolated)) ? 0. : 1.,0.)
    # Identical 64-scale workloads, with the same inputs/settings/row payload.
    endpoints=Float64.(policy["resource_Q2_GeV2"])
    scales=exp.(range(log(endpoints[1]),log(endpoints[2]);length=policy["resource_scales"]))
    resource=[];resource_values=[]
    for construction in (:independent,:batched_probes)
        plan=nothing
        cold=@elapsed plan=HoppetMatchingPlan(buildfile,queries,scales,2.4,coupling;process=:SIDIS,intervals=1024,construction=construction)
        h=initial_columns(plan.grid);value=apply_operator(plan,h)
        warm=[@elapsed apply_operator(plan,h) for _ in 1:5]
        push!(resource_values,value)
        push!(resource,Dict("construction"=>string(construction),"preparation_seconds_informational"=>cold,
            "warm_apply_seconds"=>warm,"median_warm_apply_seconds"=>median(warm),
            "row_coefficient_payload_bytes"=>plan.coefficient_payload_bytes,"row_native_max_scaled"=>plan.row_native_max_scaled))
    end
    check("identical 64-scale matched workloads",scaled(resource_values[2],resource_values[1]),policy["value_max_scaled"])
    check("identical 64-scale row payload",resource[1]["row_coefficient_payload_bytes"]==resource[2]["row_coefficient_payload_bytes"] ? 0. : 1.,0.)
    # Actual source/scale-knot-split Fourier nodes. All physical scales are used;
    # analytic initial inputs isolate evolution/matching accuracy and cost.
    model=KPSYCSS();ff=load_dss14(inputs["input_root"])
    radial=RadialRule(parton_radial_edges(model,ff,inputs["input_root"],24.);order=32)
    radial_scales=css_mu_b2.(Ref(model),radial.b)
    full=nothing
    cold=@elapsed full=HoppetMatchingPlan(buildfile,queries,radial_scales,2.4,coupling;process=:SIDIS,intervals=1024)
    h=initial_columns(full.grid);value=apply_operator(full,h)
    warm=[@elapsed apply_operator(full,h) for _ in 1:5]
    selected=unique(round.(Int,range(1,length(radial.b);length=36)))
    independent=HoppetMatchingPlan(buildfile,queries,radial_scales[selected],2.4,coupling;process=:SIDIS,intervals=1024,construction=:independent)
    comparison=apply_operator(independent,initial_columns(independent.grid))
    check("all actual radial scales use the supplied mu_b squared",full.Q2==radial_scales && length(unique(full.Q2))==length(radial.b) ? 0. : 1.,0.)
    check("actual radial workload selected independent scales",scaled(value[:,:,selected],comparison),policy["value_max_scaled"])
    truth=[matched_polynomial(x,full.tau[k],full.alpha[k],:SIDIS,:KPSY_Hhat3_homogeneous) for x in queries,k in selected]
    check("actual radial workload independent Mellin values",scaled(value[:,1,selected],truth),policy["value_max_scaled"])
    check("actual radial workload row versus native",full.row_native_max_scaled,1e-10)
    check("source/build stability",source_hashes()==sources &&
        (bytes2hex(sha256(read(buildfile))),bytes2hex(sha256(read(inputfile))))==hashes &&
        CollinsFit.verified_hoppet_build(buildfile)==build && CollinsFit.verify_parton_build(inputfile)==inputs ? 0. : 1.,0.)
    Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"contract_id"=>contract["id"],"checks"=>checks,
        "sources"=>sources,"build_sha256"=>hashes[1],"build"=>build,"input_build_sha256"=>hashes[2],
        "completed_at"=>string(now(UTC)),"julia_version"=>string(VERSION),"limits"=>contract["limits"],
        "resource_pilot"=>Dict("scope"=>"identical_analytic_64_scale_comparison","G2PERF_pass"=>false,"production_ready"=>false,
            "scales"=>length(scales),"queries"=>queries,"intervals"=>1024,"initial_columns"=>3,"measurements"=>resource),
        "radial_workload"=>Dict("scope"=>"actual_radial_scales_with_analytic_initial_inputs","production_ready"=>false,"G2PERF_pass"=>false,
            "radial_nodes"=>length(radial.b),"radial_edges"=>radial.edges,"radial_order"=>32,"radial_bound_inverse_GeV"=>24.,
            "Q2_GeV2"=>radial_scales,"selected_independent_indices"=>selected,"queries"=>queries,"intervals"=>1024,
            "preparation_seconds_informational"=>cold,"warm_apply_seconds"=>warm,"median_warm_apply_seconds"=>median(warm),
            "row_coefficient_payload_bytes"=>full.coefficient_payload_bytes,"output_bytes"=>sizeof(value)),
        "environment"=>Dict("cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,"julia_threads"=>Threads.nthreads(),
            "BLAS_threads"=>BLAS.get_num_threads(),"BLAS_configuration"=>sprint(show,BLAS.get_config())))
end

buildfile=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
inputfile=length(ARGS)<2 ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(ROOT,"runs/hoppet-batch/validation.json") : ARGS[3]
try
    report=batch_checks(buildfile,inputfile);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "failed"=>filter(c->!c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"error"=>sprint(showerror,exception)))
    rethrow()
end
