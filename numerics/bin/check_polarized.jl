# Full reconstructed structure-function validation, with analytic initial inputs.
# This producer admits no experimental data, historical fit or production result.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff
include("PolarizedReference.jl")
using .PolarizedReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/polarized-process-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(@__DIR__,"PolarizedReference.jl"),
    joinpath(@__DIR__,"LeptonicReference.jl"),joinpath(@__DIR__,"MatchedReference.jl"),joinpath(@__DIR__,"CollinearReference.jl"),
    joinpath(ROOT,"config/hoppet-contract.json"),joinpath(ROOT,"config/hoppet-point-contract.json"),
    joinpath(ROOT,"config/hoppet-batch-contract.json"),joinpath(ROOT,"config/matched-rows-contract.json"),
    joinpath(ROOT,"config/css-contract.json"),joinpath(ROOT,"config/leptonic-contract.json"),
    joinpath(ROOT,"config/initial-contract.json"),joinpath(ROOT,"fortran/hoppet_bridge.f90"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)
edges(B)=vcat([0.,1e-7,1e-6,1e-5,1e-4,.001,.005,.02,.05,.1,.2,.4,.7,1.,1.5,2.,3.,4.,6.,8.],filter(<(B),[12.,16.,20.,24.]),B)

function analytic_context(buildfile,process,model,radial;intervals=1024,target=NucleonTarget(:proton))
    q2=css_mu_b2.(Ref(model),radial.b)
    Q2=process==:SIDIS ? 4. : 110.
    plan=leptonic_plan(process,model,radial,Q2,.3,.3;z2=process==:SIA ? .5 : nothing)
    unpol1=[1+.1pid for pid in plan.flavors,b in radial.b]
    unpol2=[.8-.05pid for pid in plan.flavors,b in radial.b]
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    op=HoppetMatchingPlan(buildfile,process==:SIDIS ? [.1] : [.3,.5],q2,2.4,model.coupling;
        process=process,quantity=kind,intervals=intervals)
    bases=process==:SIDIS ? [(i==1 ? 1.3 : .8)*(1-x) for i in 1:2,x in op.grid] :
        [PolarizedReference.FF_BASES[i]*(1-x)^3 for i in 1:4,x in op.grid]
    samples=PolarizedInitialSamples(kind,op.grid,bases;input_kind=:analytic)
    if process==:SIDIS
        op2=HoppetMatchingPlan(buildfile,[.3],q2,2.4,model.coupling;process=process,intervals=intervals)
        samples2=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,
            [PolarizedReference.FF_BASES[i]*(1-x)^3 for i in 1:4,x in op2.grid];input_kind=:analytic)
        polarized_leptonic_plan(plan,unpol1,unpol2,op,samples;second_operator=op2,second_samples=samples2,target=target,charges=(1,0))
    else
        polarized_leptonic_plan(plan,unpol1,unpol2,op,samples;charges=(1,-1))
    end
end

function polarized_checks(buildfile)
    sources=source_hashes();build=CollinsFit.verified_hoppet_build(buildfile);build_hash=bytes2hex(sha256(read(buildfile)))
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"];floor=policy["relative_floor"]
    scaled(a,b)=maximum(abs.(a.-b)./max.(abs.(b),floor))
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error)&&error<=tolerance))
    theta=[.5,-.4,1.,1.,2.,2.,.3,-.2,1.,1.,0.,0.,.01];model=KPSYCSS()
    amplitude_indices=[1,2,7,8,13];a=theta[amplitude_indices]
    refinements=Dict{String,Any}();oracles=Dict{Symbol,Any}();fine=Dict{Symbol,Any}()
    for process in (:SIDIS,:SIA)
        values=[];levels=[]
        for (order,intervals) in ((8,512),(16,1024),(32,2048))
            radial=RadialRule(edges(24.);order=order)
            context=analytic_context(buildfile,process,model,radial;intervals=intervals)
            result=polarized_structure(context,theta)
            push!(values,[result.UU,result.Collins])
            push!(levels,Dict("radial_order"=>order,"intervals"=>intervals,"radial_nodes"=>length(radial.b),"values"=>last(values)))
            fine[process]=context
            println(JSON.json(Dict("completed_process"=>string(process),"radial_order"=>order,"intervals"=>intervals)));flush(stdout)
        end
        refinements[string(process)]=levels;context=fine[process]
        plan=context.leptonic
        oracle=analytic_reference(model,process,plan.Q2,plan.P,.1,plan.z1,plan.z2,context.charges...,
            context.target,a,24.;moments=true)
        oracles[process]=oracle
        check("$process complete Mellin/adaptive structure functions",scaled(values[end],[oracle.values.UU,oracle.values.Collins]),policy["physical_value_max_scaled"])
        check("$process complete radial/evolution refinement",scaled(values[end],values[end-1]),policy["physical_value_max_scaled"])
        prediction(t)=polarized_structure(context,t).Collins
        J=ForwardDiff.gradient(prediction,theta);H=ForwardDiff.hessian(prediction,theta)
        check("$process amplitude/gc derivatives versus Mellin/adaptive moments",scaled(J[amplitude_indices],oracle.gradient),policy["derivative_max_scaled"])
        check("$process amplitude/gc Hessian versus Mellin/adaptive moments",scaled(H[amplitude_indices,amplitude_indices],oracle.hessian),policy["derivative_max_scaled"])
        for c1 in (-1,0,1),c2 in (process==:SIDIS ? (0,) : (-1,0,1))
            changed=with_pion_charges(context,c1,c2)
            actual=polarized_structure(changed,theta)
            ref=analytic_reference(model,process,plan.Q2,plan.P,.1,plan.z1,plan.z2,c1,c2,context.target,a,24.)
            check("$process charge ($c1,$c2) independent flavor structure",scaled([actual.UU,actual.Collins],[ref.UU,ref.Collins]),policy["physical_value_max_scaled"])
        end
        check("$process exact zero-P Collins",abs(polarized_structure(with_momentum(context,0.),theta).Collins),0.)
        check("$process detached momentum/charge operator reuse",with_momentum(context,.5).first_operator===context.first_operator &&
            with_pion_charges(context,0).first_samples===context.first_samples ? 0. : 1.,0.)
    end
    tails=Dict{String,Any}()
    for process in (:SIDIS,:SIA)
        values=[]
        for B in policy["radial_bounds"]
            context=B==24 ? fine[process] : analytic_context(buildfile,process,model,RadialRule(edges(B);order=32);intervals=2048)
            result=polarized_structure(context,theta)
            push!(values,Dict("B_inverse_GeV"=>B,"values"=>[result.UU,result.Collins]))
        end
        tails[string(process)]=values
        check("$process full finite-bound extension 24 to 32",scaled(values[end]["values"],values[end-1]["values"]),policy["physical_value_max_scaled"])
    end
    # The physical boundary theta used here has beta=0. Inward finite differences
    # avoid leaving that parameter domain; generic BigFloat rows are independent
    # of Float64 tangent batching. Use the full coarse process, not a row alone.
    derivatives=[]
    for process in (:SIDIS,:SIA)
        context=analytic_context(buildfile,process,model,RadialRule(edges(24.);order=8);intervals=512)
        prediction(t)=polarized_structure(context,t).Collins
        gradient(t)=ForwardDiff.gradient(prediction,t)
        d=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,.001]
        J=gradient(theta);hv=ForwardDiff.jacobian(s->gradient(theta+s[1]*d),[0.])[:,1]
        for step in policy["finite_difference_steps"]
            fd,fdh=setprecision(policy["finite_difference_reference_bits"]) do
                t=BigFloat.(theta);h=BigFloat(step);direction=BigFloat.(d)
                (Float64((-3prediction(t)+4prediction(t+h*direction)-prediction(t+2h*direction))/(2h)),
                 Float64.((-3gradient(t)+4gradient(t+h*direction)-gradient(t+2h*direction))/(2h)))
            end
            e1,e2=scaled(fd,dot(J,d)),scaled(fdh,hv)
            check("$process complete 13-parameter inward JVP step $step",e1,policy["derivative_max_scaled"])
            check("$process complete 13-parameter inward HVP step $step",e2,policy["derivative_max_scaled"])
            push!(derivatives,Dict("process"=>string(process),"step"=>step,"JVP_error"=>e1,"HVP_error"=>e2,
                "radial_order"=>8,"radial_nodes"=>length(context.leptonic.radial.b),"intervals"=>512,"reference_bits"=>policy["finite_difference_reference_bits"]))
        end
    end
    check("source/build stability",source_hashes()==sources && bytes2hex(sha256(read(buildfile)))==build_hash && CollinsFit.verified_hoppet_build(buildfile)==build ? 0. : 1.,0.)
    Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"contract_id"=>contract["id"],"checks"=>checks,
        "sources"=>sources,"build_sha256"=>build_hash,"build"=>build,"refinements"=>refinements,"tails"=>tails,
        "derivative_checks"=>derivatives,"analytic_parameters"=>theta,"analytic_input_scope"=>"polynomial initial inputs and explicit analytic unpolarized pi+ bases",
        "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),"limitations"=>contract["limitations"])
end

buildfile=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
output=length(ARGS)<2 ? joinpath(ROOT,"runs/polarized/validation.json") : ARGS[2]
try
    report=polarized_checks(buildfile);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"error"=>sprint(showerror,exception)))
    rethrow()
end
