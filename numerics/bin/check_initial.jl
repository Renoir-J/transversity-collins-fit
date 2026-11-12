# Independent algebraic/derivative checks of the initial family, no experimental fit.
using CollinsFit, ForwardDiff, LinearAlgebra, JSON, SHA, Dates

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/initial-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"numerics/Project.toml"),
    joinpath(ROOT,"numerics/Manifest.toml")],filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function reference_h1(theta,pid,x,S)
    !(pid in (1,2)) && return big"0"
    x==1 && return big"0"
    t=BigFloat.(theta); u=BigFloat(x)
    N,a,b = pid==2 ? t[[1,3,5]] : t[[2,4,6]]
    power(v)=iszero(v) ? one(v) : v^v
    N*u^a*(1-u)^b*power(a+b)/(power(a)*power(b))*BigFloat(S)
end

function reference_hat(theta,pid,charge,z,ff)
    z==1 && return big"0"
    abs(pid)>=4 && return big"0"
    charge==0 && return (reference_hat(theta,pid,1,z,ff)+reference_hat(theta,pid,-1,z,ff))/2
    basis = charge==1 ? Dict(2=>1,-1=>1,1=>2,-2=>2,3=>3,-3=>4) :
                        Dict(2=>2,-1=>2,1=>1,-2=>1,3=>3,-3=>4)
    index=basis[pid]
    shape=index==1 ? (7,9,11) : (8,10,12)
    t=BigFloat.(theta);u=BigFloat(z)
    t[shape[1]]*u^t[shape[2]]*(1-u)^t[shape[3]]*BigFloat(ff[index])
end

function reference_h1_derivatives(theta,pid,x,S)
    j=zeros(13);H=zeros(13,13)
    i,a,b=pid==2 ? (1,3,5) : (2,4,6)
    n,aa,bb=theta[[i,a,b]]
    value=Float64(reference_h1(theta,pid,x,S))
    basis=Float64(reference_h1([k==i ? 1. : theta[k] for k in 1:13],pid,x,S))
    A=log(x)+log((aa+bb)/aa)
    B=log1p(-x)+log((aa+bb)/bb)
    j[i],j[a],j[b]=basis,value*A,value*B
    H[i,a]=H[a,i]=basis*A;H[i,b]=H[b,i]=basis*B
    H[a,a]=value*(A^2+inv(aa+bb)-inv(aa))
    H[b,b]=value*(B^2+inv(aa+bb)-inv(bb))
    H[a,b]=H[b,a]=value*(A*B+inv(aa+bb))
    j,H
end

function initial_checks()
    sources=source_hashes()
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    scale(a,b)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    p=kpsy_published_parameters();theta=collect(p.values)
    check("published boundary values retained exactly",theta[2]==-1 && theta[11]==theta[12]==0 ? 0. : 1.,0.)
    herrors,Haterrors,derivatives,hessians=Float64[],Float64[],Float64[],Float64[]
    ff=[3.,1.,.5,.3]
    for x in (1e-5,.0065,.1,.35,.8,.99,1.),pid in (2,1,-2,-1,3,-3,4,-4,5,-5)
        S=x*(1-x)^2
        actual=initial_transversity(p,pid,x,S)
        expected=Float64(reference_h1(theta,pid,x,S))
        push!(herrors,scale([actual],[expected]))
        if pid in (1,2) && x<1
            f(t)=initial_transversity(KPSYParameters(t),pid,x,S)
            J,H=reference_h1_derivatives(theta,pid,x,S)
            push!(derivatives,scale(ForwardDiff.gradient(f,theta),J))
            push!(hessians,scale(ForwardDiff.hessian(f,theta),H))
        end
    end
    for z in (.05,.13,.4,.75,.99,1.),charge in (-1,0,1),pid in (2,1,-2,-1,3,-3,4,-4,5,-5)
        actual=initial_collins(p,pid,charge,z,ff)
        expected=Float64(reference_hat(theta,pid,charge,z,ff))
        push!(Haterrors,scale([actual],[expected]))
    end
    check("h1 versus independent 256-bit direct-power expression",maximum(herrors),policy["formula_max_scaled"])
    check("Hhat3 flavor/charge family versus independent dictionary and powers",maximum(Haterrors),policy["formula_max_scaled"])
    check("h1 gradients versus analytic logarithmic derivatives",maximum(derivatives),policy["derivative_max_scaled"])
    check("h1 Hessians versus analytic second derivatives",maximum(hessians),policy["derivative_max_scaled"])
    positivity=maximum(abs(initial_transversity(p,pid,x,1.))-1 for x in range(.00001,.99999;length=1001) for pid in (2,1))
    check("Soffer envelope on dense initial x probes",max(0.,positivity),policy["positivity_max_relative_excess"])
    boundary=[]
    for (index,pid) in ((11,2),(12,1)),z in (.05,.3,.8)
        f(t)=initial_collins(KPSYParameters(t),pid,1,z,ff)
        ad=ForwardDiff.gradient(f,theta)[index]
        for h in (1e-4,1e-5)
            direction=zeros(13);direction[index]=h
            fd=setprecision(256) do
                (-3reference_hat(theta,pid,1,z,ff)+4reference_hat(theta+direction,pid,1,z,ff)-
                    reference_hat(theta+2direction,pid,1,z,ff))/(2BigFloat(h))
            end
            push!(boundary,scale([ad],[Float64(fd)]))
        end
    end
    check("beta=0 derivatives versus high-precision inward differences",maximum(boundary),policy["derivative_max_scaled"])
    endpoint=maximum(abs,ForwardDiff.hessian(t -> initial_collins(KPSYParameters(t),2,1,1.,zeros(4)),theta))
    check("beta=0 endpoint has exactly zero first and second derivatives",endpoint,policy["endpoint_max_abs"])
    x,z=collect(range(.01,1.;length=31)),collect(range(.05,1.;length=29))
    S=[(f==1 ? 2. : 1.)*u*(1-u)^3 for f in 1:2,u in x]
    D=[(1+.2f)*u*(1-u)^2 for f in 1:4,u in z]
    inputs=CollinsFit.initial_conditions(x,z,S,D)
    op=collinear_lo_operator(x,10.,2.4,AsymptoticTwoLoopCoupling())
    Cp=collinear_lo_operator(z,10.,2.4,AsymptoticTwoLoopCoupling();quantity=:KPSY_Hhat3_homogeneous)
    # An independent initial Jacobian at each point is then acted on by the
    # separately validated fixed linear operator; no differentiation of settings.
    for (pid,row) in ((2,7),(1,6))
        f(t)=evolve_initial_conditions(inputs,KPSYParameters(t),op,Cp).h1[row,:]
        initial_j=reduce(vcat,[permutedims(u==1 ? zeros(13) : first(reference_h1_derivatives(theta,pid,u,S[pid==2 ? 1 : 2,j])))
            for (j,u) in enumerate(x)])
        check("evolved h1 flavor $pid Jacobian against independent initial derivatives",
            scale(ForwardDiff.jacobian(f,theta),op.matrix*initial_j),policy["derivative_max_scaled"])
    end
    initial=initial_distributions(inputs,p)
    evolved=evolve_initial_conditions(inputs,p,op,Cp)
    check("matrix evolution equals independent per-flavor applications",maximum(abs.(evolved.Hhat3-
        reduce(vcat,[permutedims(apply_operator(Cp,initial.Hhat3[f,:])) for f in 1:10]))),policy["formula_max_scaled"])
    check("source hashes stable during initial check",source_hashes()==sources ? 0. : 1.,0.)
    Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL",
        "checks"=>checks,"sources"=>sources,"contract_id"=>contract["id"],"printed_parameters"=>theta,
        "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),
        "input_scope"=>"analytic fixed inputs; rounded historical parameter values; no fitted data",
        "limitations"=>contract["limitations"])
end

output=isempty(ARGS) ? joinpath(ROOT,"runs/initial/validation.json") : ARGS[1]
try
    report=initial_checks();writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
