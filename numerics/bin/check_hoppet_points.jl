# Query-bound row checks and a measured, explicitly non-production resource job.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, Statistics
include("CollinearReference.jl")
using .CollinearReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/hoppet-point-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/hoppet-input-pilot.json"),joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"numerics/bin/CollinearReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function point_checks(buildfile)
    sources=source_hashes();build=CollinsFit.verified_hoppet_build(buildfile)
    build_hash=bytes2hex(sha256(read(buildfile)))
    contract=JSON.parsefile(CONTRACT);base=JSON.parsefile(joinpath(ROOT,"config/hoppet-contract.json"))
    policy=base["validation_policy"];coupling=AsymptoticTwoLoopCoupling()
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    rejected(f)=try f();false catch;true end
    queries=[.005,.017,.1,.28,.61,.85,.95,1.,.28];q2=[1.2,2.4,10.,1000.]
    plan=HoppetPointPlan(buildfile,queries,q2,2.4,coupling;intervals=512)
    x=plan.grid;h=hcat(x.*(1 .-x).^3,-2x.^2 .* (1 .-x).^2)
    values=apply_operator(plan,h)
    reference=[xx==1 ? 0. : mellin_polynomial_tau(xx,t) for xx in queries,t in plan.tau]
    check("query-bound rows versus independent Mellin",max_scaled(values[:,1,:],reference),policy["mellin_max_scaled"])
    check("rows versus full native forward maps",plan.row_native_max_scaled,contract["validation_policy"]["row_native_max_scaled"])
    check("exact query endpoint",maximum(abs,values[8,:,:]),0.)
    check("repeated query order preserved",maximum(abs.(values[4,:,:]-values[9,:,:])),0.)
    check("shared rows counted once",plan.coefficient_payload_bytes==7*4*513*8 ? 0. : 1.,0.)
    cotangent=reshape(cos.(1:length(values)),size(values))
    adjoint=operator_adjoint(plan,cotangent)
    check("query-row scatter adjoint dot product",abs(dot(values,cotangent)-dot(h,adjoint))/max(1.,abs(dot(values,cotangent))),policy["adjoint_max_scaled"])
    check("query-row endpoint cotangent removed",maximum(abs,adjoint[end,:]),0.)
    amplitude(t)=exp(t[1])+t[2]*t[3]
    prediction(t)=vec(apply_operator(plan,amplitude(t).*h[:,1]))
    theta=[.2,-.3,.1];R=vec(reference)
    J=R*permutedims([exp(theta[1]),theta[3],theta[2]])
    check("query-row generic ForwardDiff Jacobian",max_scaled(ForwardDiff.jacobian(prediction,theta),J),policy["derivative_max_scaled"])
    loss(t)=sum(abs2,prediction(t))
    gradA=[exp(theta[1]),theta[3],theta[2]]
    HA=[exp(theta[1]) 0. 0.;0. 0. 1.;0. 1. 0.]
    expected=2dot(R,R)*(gradA*gradA'+amplitude(theta)*HA)
    check("query-row nested Hessian versus closed expression",max_scaled(ForwardDiff.hessian(loss,theta),expected),policy["derivative_max_scaled"])
    # Cached rows own no native resources. A fresh context and a second quantity
    # can be constructed while both detached plans remain numerically usable.
    native=HoppetPlan(buildfile,[10.],2.4,coupling)
    close(native)
    other=HoppetPointPlan(buildfile,[.05,.3,.9,1.],q2,2.4,coupling;
        quantity=:KPSY_Hhat3_homogeneous,intervals=256)
    check("detached plan stable across later native construction",maximum(abs.(apply_operator(plan,h)-values)),0.)
    check("second detached quantity has exact endpoint",maximum(abs,apply_operator(other,other.grid.*(1 .-other.grid).^3)[end,:]),0.)
    saved=plan.queries[2];plan.queries[2]=nextfloat(saved)
    check("mutated query rejected",rejected(() -> apply_operator(plan,h)) ? 0. : 1.,0.)
    plan.queries[2]=saved
    bad=copy(h);bad[end,1]=1.
    check("nonzero initial endpoint rejected",rejected(() -> apply_operator(plan,bad)) ? 0. : 1.,0.)
    check("query setting AD rejected",rejected(() -> ForwardDiff.derivative(t -> HoppetPointPlan(buildfile,[.1+t],q2,2.4,coupling),0.)) ? 0. : 1.,0.)
    one_loop=HoppetPointPlan(buildfile,[.01,.2,.8], [4.,100.],4.,.3;nf=3,intervals=512)
    one=apply_operator(one_loop,one_loop.grid.*(1 .-one_loop.grid).^3)
    expected_one=[mellin_polynomial(xx,q,4.,.3;nf=3) for xx in one_loop.queries,q in one_loop.Q2]
    check("scalar one-loop constructor versus Mellin",max_scaled(one,expected_one),policy["mellin_max_scaled"])
    dimensions=contract["validation_policy"]["identity_dimension_stress"]
    stress=HoppetPointPlan(buildfile,[.2],fill(2.4,dimensions["scales"]),2.4,coupling;intervals=dimensions["intervals"])
    stress_initial=stress.grid.*(1 .-stress.grid).^3
    check("2336-scale row export identity",maximum(abs.(apply_operator(stress,stress_initial).-stress_initial[1])),0.)
    check("2336-scale row storage is linear",stress.coefficient_payload_bytes==dimensions["scales"]*(dimensions["intervals"]+1)*8 ? 0. : 1.,0.)
    # This resource job measures the implemented row path at checked analytic
    # accuracy. It does not certify a physical matching/likelihood workload.
    workload=JSON.parsefile(joinpath(ROOT,"config/hoppet-input-pilot.json"))["resource_pilot"]
    scales=exp.(range(log(workload["Q2_min_GeV2"]),log(workload["Q2_max_GeV2"]);length=workload["scales"]))
    resource=nothing;cold=@elapsed resource=HoppetPointPlan(buildfile,[.0065,.1,.8],scales,2.4,coupling;intervals=1024)
    base_initial=resource.grid.*(1 .-resource.grid).^3
    initial=hcat(base_initial,-2base_initial,3base_initial)
    result=apply_operator(resource,initial)
    samples=[@elapsed apply_operator(resource,initial) for _ in 1:5]
    for k in (1,length(scales))
        truth=[mellin_polynomial_tau(xx,resource.tau[k]) for xx in resource.queries]
        check("resource workload endpoint scale $(scales[k]) analytic accuracy",max_scaled(result[:,1,k],truth),policy["mellin_max_scaled"])
    end
    resource_report=Dict("scope"=>"analytic_resource_pilot","production_ready"=>false,
        "G2PERF_pass"=>false,"distinct_queries"=>3,"scales"=>length(scales),"initial_columns"=>3,
        "intervals_per_query"=>resource.intervals,"initial_unique_nodes"=>length(resource.grid),
        "preparation_seconds_informational"=>cold,"warm_apply_seconds"=>samples,
        "median_warm_apply_seconds"=>median(samples),"row_coefficient_payload_bytes"=>resource.coefficient_payload_bytes,
        "output_bytes"=>sizeof(result),"qualification"=>workload["qualification"],
        "environment"=>Dict("cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,
            "julia_threads"=>Threads.nthreads(),"BLAS_threads"=>BLAS.get_num_threads(),
            "BLAS_configuration"=>sprint(show,BLAS.get_config())),
        "comparison"=>"No ratio to a differently accurate backend or full fitting-throughput claim")
    check("source/build hashes stable",source_hashes()==sources && bytes2hex(sha256(read(buildfile)))==build_hash &&
        CollinsFit.verified_hoppet_build(buildfile)==build ? 0. : 1.,0.)
    Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
        "contract_id"=>contract["id"],"build_sha256"=>build_hash,"build"=>build,
        "queries"=>queries,"Q2_GeV2"=>q2,"row_native_error"=>plan.row_native_max_scaled,
        "identity_dimension_stress"=>Dict("scales"=>dimensions["scales"],"nodes_per_query"=>dimensions["intervals"]+1,
            "row_coefficient_bytes"=>stress.coefficient_payload_bytes,"scope"=>dimensions["scope"],"production_ready"=>false),
        "resource_pilot"=>resource_report,"limitations"=>contract["limits"],
        "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
end

buildfile=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
output=length(ARGS)<2 ? joinpath(ROOT,"runs/hoppet-points/validation.json") : ARGS[2]
try
    report=point_checks(buildfile);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
