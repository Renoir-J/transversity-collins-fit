# Reconstructed fixed-nf LO backend, independently checked; no physical extraction.
using CollinsFit, ForwardDiff, LinearAlgebra, JSON, SHA, Dates, QuadGK, Libdl
include("CollinearReference.jl")
using .CollinearReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/hoppet-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"src/collinsflow/hoppet.py"),
    joinpath(ROOT,"scripts/prepare_hoppet.py"),joinpath(ROOT,"numerics/bin/CollinearReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function beta_mellin(x,tau,a,b)
    direction=cis(3pi/4)
    function integrand(t)
        N=1.5+t*direction
        moment=factorial(b)/prod(N+a+j for j in 0:b)
        PN=(4/3)*(1.5-2*(contour_digamma(N+1)+Float64(Base.MathConstants.eulergamma)))
        imag(direction*exp(-N*log(x)+tau*PN)*moment)/pi
    end
    quadgk(integrand,0,1,5,20,Inf;rtol=1e-11,atol=1e-13)[1]
end

function moments(plan,values;order=8)
    nodes,weights=CollinsFit.gauss_rule(order)
    queries=Float64[]; integration=Float64[]
    for i in 1:length(plan.grid)-1
        a,b=plan.grid[i:i+1]
        append!(queries,(a+b)/2 .+(b-a)/2 .*nodes)
        append!(integration,(b-a)/2 .*weights)
    end
    sampled=hoppet_evaluate(plan,values,queries)
    [dot(integration.*queries.^(n-1),sampled[:,1,k]) for n in 1:3,k in eachindex(plan.tau)]
end

function backend_checks(buildfile)
    sources=source_hashes(); build=CollinsFit.verified_hoppet_build(buildfile)
    build_sha=bytes2hex(sha256(read(buildfile)))
    contract=JSON.parsefile(CONTRACT); policy=contract["validation_policy"]
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    rejected(f)=try f(); false catch; true end
    q2=[1.2,2.4,10.,1000.]; coupling=AsymptoticTwoLoopCoupling()
    probes=Float64.(policy["point_probes"]); refinements=[]
    reference=nothing; fine=nothing; last_moments=nothing
    for dy in policy["dy_refinements"]
        plan=nothing; elapsed=@elapsed plan=HoppetPlan(buildfile,q2,2.4,coupling;dy=dy,du=.0125)
        try
            x=plan.grid; h=reshape(x.*(1 .-x).^3,:,1)
            result=apply_operator(plan,h)
            values=dropdims(hoppet_evaluate(plan,result,probes);dims=2)
            reference=[mellin_polynomial_tau(x,t) for x in probes,t in plan.tau]
            generator=hoppet_evaluate(plan,hoppet_generator(plan,1 .-x),probes)
            affine=(4/3).*(2probes.*log.(probes).+(1.5.+2log1p.(-probes)).*(1 .-probes))
            numerical=moments(plan,result)
            expected=[6/prod(n+j for j in 1:4)*exp(t*transversity_moment(n)) for n in 1:3,t in plan.tau]
            check("$dy exact Q0 identity",maximum(abs.(result[:,1,2]-h[:,1])),policy["identity_max_abs"])
            check("$dy endpoint remains zero",maximum(abs,result[end,:,:]),0.)
            push!(refinements,Dict("dy"=>dy,"unique_nodes"=>length(x),"native_nodes"=>plan.native_nodes,
                "tau"=>plan.tau,"preparation_seconds_informational"=>elapsed,
                "coefficient_payload_bytes"=>plan.coefficient_payload_bytes,
                "point_error"=>maximum(max_scaled(values[:,k],reference[:,k]) for k in eachindex(q2)),
                "generator_error"=>max_scaled(generator,affine),
                "moment_error"=>maximum(abs.(numerical-expected)./abs.(expected))))
            if dy==last(policy["dy_refinements"])
                fine=values; last_moments=numerical
                check("fine grid versus independent inverse Mellin",last(refinements)["point_error"],policy["mellin_max_scaled"])
                check("fine affine plus-prescription generator",last(refinements)["generator_error"],policy["generator_max_scaled"])
                check("fine first three Mellin moments",last(refinements)["moment_error"],policy["moment_max_relative"])
                check("moment quadrature refinement",maximum(abs.(moments(plan,result;order=12)-numerical)),1e-10)
                # Independent polynomial amplitudes make transport mistakes visible
                # in the full gradient and nested Hessian, including cross terms.
                bases=hcat(x.*(1 .-x).^3,x.^2 .* (1 .-x).^2,x.^3 .* (1 .-x).^3)
                coeff(t)=[exp(t[1]),sin(t[2]),t[2]*t[3]]
                prediction(t)=vec(hoppet_evaluate(plan,apply_operator(plan,bases*coeff(t)),probes))
                theta=[.2,-.3,.1]
                R=hcat([vec([beta_mellin(xx,t,a,b) for xx in probes,t in plan.tau]) for (a,b) in ((1,3),(2,2),(3,3))]...)
                Jcoeff=[exp(theta[1]) 0. 0.;0. cos(theta[2]) 0.;0. theta[3] theta[2]]
                J=R*Jcoeff
                check("batched ForwardDiff Jacobian versus independent Mellin amplitudes",max_scaled(ForwardDiff.jacobian(prediction,theta),J),policy["derivative_max_scaled"])
                target=collect(range(-.01,.02;length=size(R,1)))
                loss(t)=sum(abs2,prediction(t)-target)
                residual=R*coeff(theta)-target
                expected_gradient=2J'*residual
                second=hcat(R[:,1]*exp(theta[1]),-R[:,2]*sin(theta[2]),R[:,3])
                H=2J'*J
                H[1,1]+=2dot(residual,second[:,1]); H[2,2]+=2dot(residual,second[:,2])
                H[2,3]+=2dot(residual,second[:,3]); H[3,2]+=2dot(residual,second[:,3])
                check("loss gradient versus closed expression",max_scaled(ForwardDiff.gradient(loss,theta),expected_gradient),policy["derivative_max_scaled"])
                check("nested-dual Hessian versus independent expression",max_scaled(ForwardDiff.hessian(loss,theta),H),policy["derivative_max_scaled"])
                direction=[.3,-.2,.4]
                hv=ForwardDiff.jacobian(s -> ForwardDiff.gradient(loss,theta+s[1]*direction),[0.])[:,1]
                check("nested-dual HVP",max_scaled(hv,H*direction),policy["derivative_max_scaled"])
                for step in (1e-4,1e-5)
                    finite=(prediction(theta+step*direction)-prediction(theta-step*direction))/(2step)
                    check("JVP finite-difference step $step",max_scaled(finite,J*direction),policy["derivative_max_scaled"])
                end
                cot=reshape(cos.(1:length(result)),size(result))
                check("fine-grid adjoint dot product",abs(dot(result,cot)-dot(h,operator_adjoint(plan,cot)))/max(1.,abs(dot(result,cot))),policy["adjoint_max_scaled"])
                check("second simultaneous context rejected",rejected(() -> HoppetPlan(buildfile,q2,2.4,coupling)) ? 0. : 1.,0.)
                check("nonzero endpoint rejected",rejected(() -> apply_operator(plan,ones(length(x)))) ? 0. : 1.,0.)
                check("NaN input rejected",rejected(() -> apply_operator(plan,fill(NaN,length(x)))) ? 0. : 1.,0.)
                check("query below support rejected",rejected(() -> hoppet_evaluate(plan,h,[first(x)/2])) ? 0. : 1.,0.)
                check("query one ULP below declared support rejected",rejected(() -> hoppet_evaluate(plan,h,[prevfloat(first(x))])) ? 0. : 1.,0.)
                saved=plan.grid[2];plan.grid[2]=nextfloat(saved)
                check("mutated grid rejected",rejected(() -> apply_operator(plan,h)) ? 0. : 1.,0.)
                plan.grid[2]=saved
            end
        finally
            close(plan)
        end
        check("closed context rejected",rejected(() -> apply_operator(plan,zeros(length(plan.grid)))) ? 0. : 1.,0.)
    end
    boundary=HoppetPlan(buildfile,[2.4],2.4,coupling;xmin=.05)
    try
        check("declared DSS boundary retained exactly",first(boundary.grid)==.05 ? 0. : 1.,0.)
        value=hoppet_evaluate(boundary,1 .-boundary.grid,[.05,1.])
        check("exact declared boundary query",max_scaled(value,[.95,0.]),policy["identity_max_abs"])
    finally
        close(boundary)
    end
    # Step-size and requested quadrature epsilon are varied independently of mesh.
    solver=[]
    for (du,eps) in ((.05,1e-10),(.025,1e-10),(.0125,1e-8))
        plan=HoppetPlan(buildfile,q2,2.4,coupling;du=du,eps=eps)
        try
            values=hoppet_evaluate(plan,apply_operator(plan,plan.grid.*(1 .-plan.grid).^3),probes)
            push!(solver,Dict("du"=>du,"requested_eps"=>eps,"difference_from_fine"=>max_scaled(values,fine)))
        finally
            close(plan)
        end
    end
    check("RK step refinement at adjacent fine steps",solver[2]["difference_from_fine"],policy["composition_max_scaled"])
    check("coarse RK step change is effective",solver[1]["difference_from_fine"]>1e-12 ? 0. : 1.,0.)
    check("requested convolution epsilon stability",solver[3]["difference_from_fine"],policy["derivative_max_scaled"])
    # Full ordinate-basis transpose comparisons do not assume smooth inputs.
    for layout in (:single,:nested)
        plan=HoppetPlan(buildfile,[2.4,10.],2.4,coupling;dy=.2,layout=layout)
        try
            n=length(plan.grid); basis=vcat(Matrix{Float64}(I,n-1,n-1),zeros(1,n-1))
            U=apply_operator(plan,basis)
            for node in (1,n÷2,n)
                rows=CollinsFit.hoppet_rows(plan,node)
                reference=permutedims(U[node,:,:])
                check("$layout exported node $node rows versus full constrained basis",max_scaled(rows[:,1:end-1],reference),policy["adjoint_max_scaled"])
            end
            g=reshape(sin.(1:n*2),n,1,2)
            expected=sum(transpose(U[:,:,k])*g[:,1,k] for k in 1:2)
            actual=operator_adjoint(plan,g)[:,1]
            check("$layout transpose versus complete constrained forward basis",max_scaled(actual[1:end-1],expected),policy["adjoint_max_scaled"])
            check("$layout endpoint cotangent projected out",abs(actual[end]),0.)
            check("$layout invalid row node rejected",rejected(() -> CollinsFit.hoppet_rows(plan,0)) ? 0. : 1.,0.)
        finally
            close(plan)
        end
    end
    check("setting AD rejected",rejected(() -> ForwardDiff.derivative(t -> HoppetPlan(buildfile,[10.0+t],2.4,coupling),0.)) ? 0. : 1.,0.)
    check("generic Collins convention rejected",rejected(() -> HoppetPlan(buildfile,q2,2.4,coupling;quantity=:Collins)) ? 0. : 1.,0.)
    check("coupling nf mismatch rejected",rejected(() -> HoppetPlan(buildfile,q2,2.4,coupling;nf=3)) ? 0. : 1.,0.)
    check("source/build stable during validation",source_hashes()==sources && bytes2hex(sha256(read(buildfile)))==build_sha &&
        CollinsFit.verified_hoppet_build(buildfile)==build ? 0. : 1.,0.)
    Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
        "build_sha256"=>build_sha,"build"=>build,"contract_id"=>contract["id"],
        "Q2_GeV2"=>q2,"probes"=>probes,"refinements"=>refinements,"solver_refinements"=>solver,
        "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),
        "limitations"=>["The original modified HOPPET and fitter remain unavailable",
        "Requested upstream quadrature epsilon is a solver setting, not an attained error bound; upstream precision warnings are retained in job stdout and the independent epsilon-stability check is required",
        "No implicit scale interpolation, physical many-scale assembly, historical fit or production performance/coverage admission"])
end

buildfile=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
output=length(ARGS)<2 ? joinpath(ROOT,"runs/hoppet/validation.json") : ARGS[2]
try
    report=backend_checks(buildfile);writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
