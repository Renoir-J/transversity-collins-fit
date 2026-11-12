# Shared report plumbing. Numerical work remains in complete Julia jobs.
include("PolarizedReference.jl")
include("PolarizedFixture.jl")
include("ObservableReference.jl")
using .ObservableReference

const OBS_ROOT=normpath(joinpath(@__DIR__,"../.."))
const OBS_CONTRACT=joinpath(OBS_ROOT,"config/leptonic-observable-contract.json")
const OBS_BINS=((.05,.2),(.1,.3),(.2,.4))
const OBS_Y=.4
const OBS_ANGULAR=.6
const OBS_NAMES=Dict(:SIDIS=>["structure_function_ratio","sin_phi_sum_amplitude"],:SIA=>["UL","UC"])

function observable_source_hashes(producer)
    files=vcat([producer,OBS_CONTRACT],
        [joinpath(@__DIR__,name*".jl") for name in ("ObservableChecks","ObservableReference","PolarizedFixture","PolarizedReference","LeptonicReference","MatchedReference","CollinearReference")],
        [joinpath(OBS_ROOT,"config/"*name*"-contract.json") for name in
            ("polarized-process","hoppet","hoppet-point","hoppet-batch","matched-rows","css","leptonic","initial","parton")],
        [joinpath(OBS_ROOT,"config/parton-sources.json"),joinpath(OBS_ROOT,"fortran/hoppet_bridge.f90"),
         joinpath(OBS_ROOT,"numerics/Project.toml"),joinpath(OBS_ROOT,"numerics/Manifest.toml")],
        filter(p->endswith(p,".jl"),readdir(joinpath(OBS_ROOT,"numerics/src");join=true)))
    Dict(relpath(p,OBS_ROOT)=>bytes2hex(sha256(read(p))) for p in unique(files))
end

function observable_pair(context;interval=nothing,order=128,source="reconstructed validation fixture; no experimental admission")
    momentum=interval===nothing ? nothing : MomentumRule(interval...;order=order)
    if context.leptonic.process==:SIDIS
        [sidis_observable(context;definition=:structure_function_ratio,momentum=momentum,source=source),
         sidis_observable(context;definition=:sin_phi_sum_amplitude,y=OBS_Y,momentum=momentum,source=source)]
    else
        angular=SIAAngularFactor(OBS_ANGULAR;source="fixed validation ratio of angular means; not an experimental average")
        [sia_observable(context;ratio=kind,angular=angular,momentum=momentum,source=source) for kind in (:UL,:UC)]
    end
end
observable_values(pair,theta)=[predict_observable(o,theta) for o in pair]
all_observable_values(context,theta)=reduce(vcat,[observable_values(observable_pair(context;interval=interval),theta) for interval in (nothing,OBS_BINS...)])
observable_scaled(a,b,policy)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
function observable_check!(checks,name,error,tolerance)
    push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error)&&error<=tolerance))
end

function momentum_refinements(context,theta,policy,checks)
    result=[];process=context.leptonic.process
    for interval in OBS_BINS
        levels=[]
        for order in policy["momentum_orders"]
            values=observable_values(observable_pair(context;interval=interval,order=order),theta)
            error=isempty(levels) ? nothing : observable_scaled(values,last(levels)["values"],policy)
            accepted=order>=policy["minimum_momentum_order"] && error!==nothing && error<=policy["momentum_refinement_max_scaled"]
            push!(levels,Dict("order"=>order,"values"=>values,"previous_order_max_error"=>error,"refinement_accepted"=>accepted))
            accepted && break
        end
        observable_check!(checks,"$process P interval $interval momentum refinement admitted",last(levels)["refinement_accepted"] ? 0. : 1.,0.)
        observable_check!(checks,"$process P interval $interval adjacent momentum difference",last(levels)["previous_order_max_error"],policy["momentum_refinement_max_scaled"])
        push!(result,Dict("bounds_GeV"=>collect(interval),"levels"=>levels))
    end
    result
end

function observable_derivatives(context,theta,policy,checks)
    process=context.leptonic.process;interval=OBS_BINS[2]
    pair=observable_pair(context;interval=interval)
    predict(t)=observable_values(pair,t)
    generic=prepare_generic_reference(context;momentum_interval=interval,y=OBS_Y,angular=OBS_ANGULAR)
    reference=generic.predict
    observable_check!(checks,"$process scalar/adaptive reference at fixed rows",observable_scaled(predict(theta),reference(theta),policy),policy["algebra_max_scaled"])
    d=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,.001]
    J=ForwardDiff.jacobian(predict,theta)
    hv=reshape(ForwardDiff.jacobian(s->vec(ForwardDiff.jacobian(predict,theta+s[1]*d)),[0.])[:,1],size(J))
    observable_check!(checks,"$process full 13-parameter observable gradient/HVP finite",all(isfinite,J)&&all(isfinite,hv) ? 0. : 1.,0.)
    result=[]
    for step in policy["finite_difference_steps"]
        fd,fdh=setprecision(policy["finite_difference_reference_bits"]) do
            t=BigFloat.(theta);h=BigFloat(step);direction=BigFloat.(d)
            gradient(x)=ForwardDiff.jacobian(reference,x)
            (Float64.((-3reference(t)+4reference(t+h*direction)-reference(t+2h*direction))/(2h)),
             Float64.((-3gradient(t)+4gradient(t+h*direction)-gradient(t+2h*direction))/(2h)))
        end
        jvp=J*d;e1=observable_scaled(fd,jvp,policy);e2=observable_scaled(fdh,hv,policy)
        observable_check!(checks,"$process observable full-13 inward JVP step $step",e1,policy["derivative_max_scaled"])
        observable_check!(checks,"$process observable full-13 inward HVP step $step",e2,policy["derivative_max_scaled"])
        push!(result,Dict("process"=>string(process),"definitions"=>OBS_NAMES[process],"step"=>step,
            "JVP_error"=>e1,"HVP_error"=>e2,"actual_JVP"=>jvp,"reference_JVP"=>fd,
            "actual_HVP"=>hv,"reference_HVP"=>fdh,"direction"=>d,"parameters"=>theta,
            "momentum_bounds_GeV"=>collect(interval),"momentum_order"=>128,"radial_order"=>context.leptonic.radial.order,
            "radial_nodes"=>length(context.leptonic.radial.b),"intervals"=>context.first_operator.intervals,
            "reference_bits"=>policy["finite_difference_reference_bits"],"reference_scope"=>generic.scope))
        println(JSON.json(Dict("completed_process"=>string(process),"derivative_step"=>step,"JVP_error"=>e1,"HVP_error"=>e2)));flush(stdout)
    end
    result
end

function observable_result(checks,sources,contract;kwargs...)
    merge(Dict{String,Any}("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL",
        "scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "contract_id"=>contract["id"],"checks"=>checks,"sources"=>sources,"julia_version"=>string(VERSION),
        "completed_at"=>string(now(UTC)),"momentum_bounds_GeV"=>collect.(OBS_BINS),"physical_grid_momentum_order"=>128,
        "SIDIS_y"=>OBS_Y,"SIA_angular_factor"=>OBS_ANGULAR,"definitions"=>Dict(string(k)=>v for (k,v) in OBS_NAMES),
        "limitations"=>contract["limitations"]),Dict{String,Any}(string(k)=>v for (k,v) in kwargs))
end
