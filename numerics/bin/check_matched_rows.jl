# Evolved/matched polarized legs only; no observable, likelihood or production fit.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, QuadGK, Statistics
include("MatchedReference.jl")
using .MatchedReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/matched-rows-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/hoppet-point-contract.json"),joinpath(ROOT,"config/leptonic-contract.json"),
    joinpath(ROOT,"config/evolution-contract.json"),joinpath(ROOT,"fortran/hoppet_bridge.f90"),
    joinpath(ROOT,"src/collinsflow/hoppet.py"),joinpath(ROOT,"scripts/prepare_hoppet.py"),
    joinpath(ROOT,"numerics/bin/MatchedReference.jl"),joinpath(ROOT,"numerics/bin/CollinearReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function matched_checks(buildfile)
    sources=source_hashes();build=CollinsFit.verified_hoppet_build(buildfile)
    build_hash=bytes2hex(sha256(read(buildfile)));contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    rejected(f)=try f();false catch;true end
    scaled(a,b)=maximum(abs.(a-b)./max.(abs.(b),.01maximum(abs,b),1e-15))
    for N in ComplexF64[.2,.5,1,2,3,2+3im,.5+10im]
        direct=quadgk(v -> v^N*log(v)/(1-v),0.,1.;rtol=1e-12,atol=1e-14)[1]
        check("independent regular moment N=$N",abs(direct+contour_trigamma(N+1)),policy["kernel_moment_max_abs"])
    end
    for z in ComplexF64[.2+.3im,2+5im,-3.2+1.2im,1.5+1000im]
        check("trigamma continuation recurrence z=$z",abs(contour_trigamma(z)-contour_trigamma(z+1)-inv(z^2)),policy["kernel_moment_max_abs"])
    end
    for process in (:SIDIS,:SIA),x in (.005,.02,.1,.3,.6,.85,.95)
        rule=matching_rule(process,:collins,x,.327;order=64)
        direct=match_leg(rule,[y*(1-y)^3 for y in rule.coordinates])
        oracle=matched_polynomial(x,0.,.327,process,:KPSY_Hhat3_homogeneous)
        check("Q0 independent quadrature $process x=$x",abs(direct-oracle)/max(abs(oracle),1e-8),policy["analytic_Q0_matching_max_scaled"])
    end
    for tau in (-.05,.05,.15),x in (.02,.3,.85)
        a=matched_polynomial(x,tau,.327,:SIDIS,:KPSY_Hhat3_homogeneous)
        b=matched_polynomial(x,tau,.327,:SIDIS,:KPSY_Hhat3_homogeneous;angle=2pi/3)
        check("independent inversion contour tau=$tau x=$x",abs(a-b)/max(abs(b),1e-8),policy["analytic_Q0_matching_max_scaled"])
    end
    coupling=AsymptoticTwoLoopCoupling();q2=Float64.(policy["Q2_GeV2"])
    queries=[.005,.02,.1,.3,.6,.85,.95,1.,.3]
    plans=HoppetMatchingPlan[];snapshots=[]
    for (process,quantity) in ((:SIDIS,:transversity),(:SIDIS,:KPSY_Hhat3_homogeneous),(:SIA,:KPSY_Hhat3_homogeneous))
        label="$process $quantity"
        plan=HoppetMatchingPlan(buildfile,queries,q2,2.4,coupling;process=process,quantity=quantity,intervals=1024)
        x=plan.grid;h=hcat(x.*(1 .-x).^3,x.^2 .* (1 .-x).^2,x.^3 .* (1 .-x).^3)
        values=apply_operator(plan,h)
        function reference(a,b)
            [xx==1 ? 0. : matched_polynomial(xx,t,alpha_s(coupling,q),process,quantity;a=a,b=b)
                for xx in queries,(t,q) in zip(plan.tau,q2)]
        end
        R=hcat([vec(reference(a,b)) for (a,b) in ((1,3),(2,2),(3,3))]...)
        check("$label independent Mellin values",scaled(reshape(values,length(queries),3,length(q2))[:,1,:],reference(1,3)),policy["analytic_value_max_scaled"])
        check("$label rows versus full native matching",plan.row_native_max_scaled,policy["row_native_max_scaled"])
        check("$label physical alpha metadata",maximum(abs.(plan.alpha-alpha_s.(Ref(coupling),q2))),0.)
        check("$label exact query endpoint",maximum(abs,values[8,:,:]),0.)
        check("$label repeated query preserved",maximum(abs.(values[4,:,:]-values[9,:,:])),0.)
        check("$label shared row storage",plan.coefficient_payload_bytes==7*4*1025*8 ? 0. : 1.,0.)
        coeff(t)=[exp(t[1]),sin(t[2]),t[2]*t[3]]
        prediction(t)=vec(apply_operator(plan,h*coeff(t)))
        theta=[.2,-.3,.1];J=R*[exp(theta[1]) 0. 0.;0. cos(theta[2]) 0.;0. theta[3] theta[2]]
        check("$label Jacobian versus independent Mellin",scaled(ForwardDiff.jacobian(prediction,theta),J),policy["derivative_max_scaled"])
        target=collect(range(-.01,.02;length=size(R,1)))
        loss(t)=sum(abs2,prediction(t)-target)
        residual=R*coeff(theta)-target
        H=2J'*J
        H[1,1]+=2dot(residual,R[:,1]*exp(theta[1]));H[2,2]-=2dot(residual,R[:,2]*sin(theta[2]))
        H[2,3]+=2dot(residual,R[:,3]);H[3,2]+=2dot(residual,R[:,3])
        check("$label gradient versus independent Mellin",scaled(ForwardDiff.gradient(loss,theta),2J'*residual),policy["derivative_max_scaled"])
        check("$label nested Hessian versus independent Mellin",scaled(ForwardDiff.hessian(loss,theta),H),policy["derivative_max_scaled"])
        direction=[.3,-.2,.4]
        hv=ForwardDiff.jacobian(s -> ForwardDiff.gradient(loss,theta+s[1]*direction),[0.])[:,1]
        check("$label HVP versus independent Mellin",scaled(hv,H*direction),policy["derivative_max_scaled"])
        for step in policy["finite_difference_steps"]
            fd=(prediction(theta+step*direction)-prediction(theta-step*direction))/(2step)
            check("$label JVP difference step $step",scaled(fd,J*direction),policy["derivative_max_scaled"])
        end
        cotangent=reshape(cos.(1:length(values)),size(values))
        adjoint=operator_adjoint(plan,cotangent)
        check("$label scatter adjoint",abs(dot(values,cotangent)-dot(h,adjoint))/max(1.,abs(dot(values,cotangent))),policy["row_native_max_scaled"])
        check("$label constrained endpoint adjoint",maximum(abs,adjoint[end,:]),0.)
        saved=plan.alpha[1];plan.alpha[1]=nextfloat(saved)
        check("$label changed physical alpha rejected",rejected(() -> apply_operator(plan,h)) ? 0. : 1.,0.)
        plan.alpha[1]=saved
        push!(plans,plan);push!(snapshots,(initial=h,values=values))
    end
    for (plan,snapshot) in zip(plans,snapshots)
        check("$(plan.process) $(plan.quantity) detached after later contexts",maximum(abs.(apply_operator(plan,snapshot.initial)-snapshot.values)),0.)
    end
    check("SIA transversity rejected",rejected(() -> HoppetMatchingPlan(buildfile,[.1],q2,2.4,coupling;process=:SIA,quantity=:transversity)) ? 0. : 1.,0.)
    check("operator scale AD rejected",rejected(() -> ForwardDiff.derivative(t -> HoppetMatchingPlan(buildfile,[.1],[10.0+t],2.4,coupling;process=:SIDIS),0.)) ? 0. : 1.,0.)
    native=HoppetPlan(buildfile,q2,2.4,coupling;xmin=.05,dy=-log(.05)/32,layout=:single)
    try
        n=length(native.grid);basis=Matrix{Float64}(I,n,n-1)
        rows=CollinsFit.hoppet_rows(native,1;collins=true)
        full=CollinsFit.hoppet_collins_regular(native,apply_operator(native,basis))
        check("regular rows complete constrained basis",scaled(rows[:,1:end-1],permutedims(full[1,:,:])),policy["row_native_max_scaled"])
        check("regular row constrained endpoint",maximum(abs,CollinsFit.hoppet_rows(native,n;collins=true)),0.)
    finally
        close(native)
    end
    one=HoppetMatchingPlan(buildfile,[.01,.2,.8],[4.,100.],4.,.3;nf=3,process=:SIDIS,intervals=512)
    expected=[matched_polynomial(x,t,alpha_lo(q,4.,.3;nf=3),:SIDIS,:KPSY_Hhat3_homogeneous)
        for x in one.queries,(t,q) in zip(one.tau,one.Q2)]
    check("scalar one-loop physical matching",scaled(apply_operator(one,one.grid.*(1 .-one.grid).^3),expected),policy["analytic_value_max_scaled"])
    dimensions=policy["identity_dimension_stress"]
    stress=HoppetMatchingPlan(buildfile,[.2],fill(2.4,dimensions["scales"]),2.4,coupling;process=:SIA,intervals=dimensions["intervals"])
    identity_values=apply_operator(stress,stress.grid.*(1 .-stress.grid).^3)
    truth=matched_polynomial(.2,0.,alpha_s(coupling,2.4),:SIA,:KPSY_Hhat3_homogeneous)
    check("2336 Q0 matched scales",maximum(abs.(identity_values.-truth))/abs(truth),policy["analytic_value_max_scaled"])
    check("2336 Q0 matching storage is linear",stress.coefficient_payload_bytes==2336*513*8 ? 0. : 1.,0.)
    workload=policy["resource_pilot"]
    scales=exp.(range(log(workload["Q2_min_GeV2"]),log(workload["Q2_max_GeV2"]);length=workload["scales"]))
    resource=nothing;cold=@elapsed resource=HoppetMatchingPlan(buildfile,Float64.(workload["queries"]),scales,2.4,coupling;process=:SIDIS,intervals=workload["intervals"])
    base=resource.grid.*(1 .-resource.grid).^3;initial=hcat(base,-2base,3base)
    result=apply_operator(resource,initial);samples=[@elapsed apply_operator(resource,initial) for _ in 1:5]
    for k in (1,length(scales))
        truth=[matched_polynomial(x,resource.tau[k],resource.alpha[k],:SIDIS,:KPSY_Hhat3_homogeneous) for x in resource.queries]
        check("resource endpoint scale $(scales[k]) Mellin",scaled(result[:,1,k],truth),policy["analytic_value_max_scaled"])
    end
    check("source/build stability",source_hashes()==sources && bytes2hex(sha256(read(buildfile)))==build_hash && CollinsFit.verified_hoppet_build(buildfile)==build ? 0. : 1.,0.)
    Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"contract_id"=>contract["id"],"checks"=>checks,
        "sources"=>sources,"build_sha256"=>build_hash,"build"=>build,"julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),
        "identity_dimension_stress"=>Dict("scales"=>2336,"row_coefficient_bytes"=>stress.coefficient_payload_bytes,"scope"=>dimensions["scope"],"production_ready"=>false),
        "resource_pilot"=>Dict("scope"=>"analytic_matched_row_resource_pilot","production_ready"=>false,"G2PERF_pass"=>false,
            "scales"=>length(scales),"distinct_queries"=>length(resource.queries),"initial_columns"=>3,
            "preparation_seconds_informational"=>cold,"warm_apply_seconds"=>samples,"median_warm_apply_seconds"=>median(samples),
            "row_coefficient_payload_bytes"=>resource.coefficient_payload_bytes,"output_bytes"=>sizeof(result),
            "environment"=>Dict("cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,"julia_threads"=>Threads.nthreads(),
                "BLAS_threads"=>BLAS.get_num_threads(),"BLAS_configuration"=>sprint(show,BLAS.get_config()))),"limitations"=>contract["limits"])
end

buildfile=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
output=length(ARGS)<2 ? joinpath(ROOT,"runs/matched-rows/validation.json") : ARGS[2]
try
    report=matched_checks(buildfile);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"error"=>sprint(showerror,exception)))
    rethrow()
end
