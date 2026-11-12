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
const OBS_MEAN_MEASURE="synthetic fixed event weights [0.6,0.4] at cos_squared(theta) [0,1]; fixed other kinematics"
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

function observable_pair(context;interval=nothing,order=128,source="reconstructed validation fixture; no experimental admission",angular_kind=:ratio_of_angular_means)
    momentum=interval===nothing ? nothing : MomentumRule(interval...;order=order)
    if context.leptonic.process==:SIDIS
        [sidis_observable(context;definition=:structure_function_ratio,momentum=momentum,source=source),
         sidis_observable(context;definition=:sin_phi_sum_amplitude,y=OBS_Y,momentum=momentum,source=source)]
    else
        angular_kind in (:ratio_of_angular_means,:mean_of_angular_ratio) || error("unknown validation angular source kind")
        angular=angular_kind==:ratio_of_angular_means ?
            SIAAngularFactor(OBS_ANGULAR;source="fixed validation ratio of angular means; not an experimental average") :
            SIAMeanOfRatio(OBS_ANGULAR;source="synthetic event angular average; no experimental admission",averaging_measure=OBS_MEAN_MEASURE)
        [sia_observable(context;ratio=kind,angular=angular,momentum=momentum,source=source) for kind in (:UL,:UC)]
    end
end
observable_values(pair,theta)=[predict_observable(o,theta) for o in pair]
all_observable_values(context,theta;angular_kind=:ratio_of_angular_means)=reduce(vcat,[observable_values(observable_pair(context;interval=interval,angular_kind=angular_kind),theta) for interval in (nothing,OBS_BINS...)])
observable_scaled(a,b,policy)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
function observable_check!(checks,name,error,tolerance)
    push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error)&&error<=tolerance))
end

function momentum_refinements(context,theta,policy,checks;angular_kind=:ratio_of_angular_means)
    result=[];process=context.leptonic.process
    label=angular_kind==:mean_of_angular_ratio ? "$process mean-of-ratio" : string(process)
    for interval in OBS_BINS
        levels=[]
        for order in policy["momentum_orders"]
            values=observable_values(observable_pair(context;interval=interval,order=order,angular_kind=angular_kind),theta)
            error=isempty(levels) ? nothing : observable_scaled(values,last(levels)["values"],policy)
            accepted=order>=policy["minimum_momentum_order"] && error!==nothing && error<=policy["momentum_refinement_max_scaled"]
            push!(levels,Dict("order"=>order,"values"=>values,"previous_order_max_error"=>error,"refinement_accepted"=>accepted))
            accepted && break
        end
        observable_check!(checks,"$label P interval $interval momentum refinement admitted",last(levels)["refinement_accepted"] ? 0. : 1.,0.)
        observable_check!(checks,"$label P interval $interval adjacent momentum difference",last(levels)["previous_order_max_error"],policy["momentum_refinement_max_scaled"])
        push!(result,Dict("bounds_GeV"=>collect(interval),"levels"=>levels))
    end
    result
end

function angular_source_level(context,theta,values,level,policy,checks)
    context.leptonic.process==:SIA || return nothing
    supplied=all_observable_values(context,theta;angular_kind=:mean_of_angular_ratio)
    observable_check!(checks,"SIA level $level equal effective angular factor across source kinds",
        observable_scaled(supplied,values,policy),policy["algebra_max_scaled"])
    supplied
end

"Independent weighted-angle arithmetic; this is an explicitly synthetic measure check."
function angular_fixture_reference(contract,policy,checks)
    fixture=contract["synthetic_angular_fixture"]
    fixture["scope"]=="explicitly synthetic mathematical averaging check" || error("synthetic angular fixture scope required")
    reference=setprecision(256) do
        c=BigFloat.(fixture["cos_squared"]);s=1 .-c
        uniform=BigFloat.(fixture["uniform_weights"]);event=BigFloat.(fixture["event_weights"])
        acceptance=BigFloat.(fixture["paired_acceptance_weights"])
        all(w->all(isfinite,w)&&all(>=(0),w)&&abs(sum(w)-1)<=BigFloat(policy["algebra_max_scaled"]),(uniform,event,acceptance)) || error("normalized nonnegative synthetic angular weights required")
        all(x->0<=x<=1,c) || error("synthetic cos-squared range")
        weighted=acceptance.*(1 .+c);weighted./=sum(weighted)
        (uniform_mean=sum(uniform.*s./(1 .+c)),uniform_ratio=sum(uniform.*s)/(1+sum(uniform.*c)),
         event_mean=sum(event.*s./(1 .+c)),acceptance_ratio=sum(acceptance.*s)/(1+sum(acceptance.*c)),
         derived_event_weights=weighted,event_weights=event)
    end
    identity_error=maximum(abs.(Float64.([reference.event_mean-reference.acceptance_ratio,
        reference.event_mean-OBS_ANGULAR,reference.event_mean-fixture["effective_factor"],
        (reference.derived_event_weights.-reference.event_weights)...])))
    observable_check!(checks,"SIA synthetic event/acceptance angular measure identity",identity_error,policy["algebra_max_scaled"])
    separate=maximum(abs.(Float64.([reference.uniform_mean-BigFloat(1)/2,reference.uniform_ratio-BigFloat(1)/3])))
    observable_check!(checks,"SIA same-weight angular means retain distinct definitions",separate,policy["algebra_max_scaled"])
    fixture["averaging_measure"]==OBS_MEAN_MEASURE || error("synthetic angular measure descriptor changed")
    merge(Dict{String,Any}(fixture),Dict("uniform_mean_of_ratio"=>Float64(reference.uniform_mean),
        "uniform_ratio_of_means"=>Float64(reference.uniform_ratio),"event_mean_of_ratio"=>Float64(reference.event_mean),
        "acceptance_ratio_of_means"=>Float64(reference.acceptance_ratio),
        "derived_event_weights"=>Float64.(reference.derived_event_weights),"reference_bits"=>256))
end

function angular_source_checks(context,coarse,theta,references,derivatives,contract,checks;reference_tolerance,reference_scope)
    context.leptonic.process==coarse.leptonic.process==:SIA || error("SIA source-kind context required")
    policy=contract["validation_policy"]
    fixture=angular_fixture_reference(contract,policy,checks)
    rows=filter(r->r["process"]=="SIA",references)
    length(rows)==4 || error("all point and finite-P source references required")
    comparisons=[]
    for row in rows
        interval=row["momentum_bounds_GeV"]
        pair=observable_pair(context;interval=interval,angular_kind=:mean_of_angular_ratio)
        actual=observable_values(pair,theta)
        observable_check!(checks,"SIA mean-of-ratio interval $interval independent composition",
            observable_scaled(actual,row["reference_values"],policy),reference_tolerance)
        metadata=[Dict(k=>o.provenance[k] for k in ("contract_id","scope","production_ready","physical_extraction",
            "process","definition","angular_definition","angular_source","angular_factor","angular_averaging_measure")) for o in pair]
        valid=all(m->m["contract_id"]==contract["id"] && m["scope"]=="reconstructed_theory" &&
            m["production_ready"]===false && m["physical_extraction"]===false &&
            m["angular_definition"]=="mean_of_angular_ratio" && m["angular_factor"]==OBS_ANGULAR &&
            m["angular_averaging_measure"]==fixture["averaging_measure"],metadata)
        observable_check!(checks,"SIA mean-of-ratio interval $interval source metadata",valid ? 0. : 1.,0.)
        push!(comparisons,Dict("process"=>"SIA","momentum_bounds_GeV"=>interval,"values"=>actual,
            "reference_values"=>row["reference_values"],"reference_scope"=>reference_scope,"metadata"=>metadata))
    end
    old=filter(r->r["process"]=="SIA",derivatives)
    length(old)==length(policy["finite_difference_steps"]) || error("complete source-kind derivative references required")
    pair=observable_pair(coarse;interval=OBS_BINS[2],angular_kind=:mean_of_angular_ratio)
    prediction(t)=observable_values(pair,t)
    J=ForwardDiff.jacobian(prediction,theta);d=old[1]["direction"]
    hv=reshape(ForwardDiff.jacobian(s->vec(ForwardDiff.jacobian(prediction,theta+s[1]*d)),[0.])[:,1],size(J))
    derived=[]
    for row in old
        row["parameters"]==theta && row["direction"]==d && row["momentum_order"]==128 &&
            row["momentum_bounds_GeV"]==collect(OBS_BINS[2]) &&
            row["reference_bits"]==policy["finite_difference_reference_bits"] || error("source-kind FD reference binding mismatch")
        step=row["step"];jvp=J*d
        e1=observable_scaled(row["reference_JVP"],jvp,policy);e2=observable_scaled(row["reference_HVP"],hv,policy)
        observable_check!(checks,"SIA mean-of-ratio full-13 inward JVP step $step",e1,policy["derivative_max_scaled"])
        observable_check!(checks,"SIA mean-of-ratio full-13 inward HVP step $step",e2,policy["derivative_max_scaled"])
        push!(derived,Dict("process"=>"SIA","angular_definition"=>"mean_of_angular_ratio","step"=>step,
            "actual_JVP"=>jvp,"reference_JVP"=>row["reference_JVP"],"actual_HVP"=>hv,"reference_HVP"=>row["reference_HVP"],
            "JVP_error"=>e1,"HVP_error"=>e2,"direction"=>d,"parameters"=>theta,"momentum_bounds_GeV"=>collect(OBS_BINS[2]),
            "reference_angular_factor"=>OBS_ANGULAR,"reference_scope"=>row["reference_scope"],
            "actual_angular_factor"=>pair[1].angular,"averaging_measure"=>OBS_MEAN_MEASURE,
            "reference_bits"=>row["reference_bits"],"momentum_order"=>128,
            "radial_order"=>coarse.leptonic.radial.order,"radial_nodes"=>length(coarse.leptonic.radial.b),
            "intervals"=>coarse.first_operator.intervals,
            "reference_reused_from_equal_effective_factor"=>true))
    end
    Dict("fixture_reference"=>fixture,"point_and_bin_comparisons"=>comparisons,"derivative_checks"=>derived,
        "momentum_refinements"=>momentum_refinements(context,theta,policy,checks;angular_kind=:mean_of_angular_ratio),
        "reference_reuse_policy"=>"Equal fixed effective angular factors share the independent scalar/adaptive FD arrays; both typed prediction paths are evaluated separately.")
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
