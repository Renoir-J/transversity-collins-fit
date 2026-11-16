"Prospectively frozen resource measurements, never a full132 diagnostic gate."
const HISTORICAL_SIDIS_PILOT_SCOPE="historical_reconstructed_sidis_resource_pilot"
const HISTORICAL_SIDIS_PILOT_POLICY=normpath(joinpath(@__DIR__,"../../config/historical-sidis-pilot-policy.json"))
const HISTORICAL_SIDIS_PILOT_SHA="395854b73f595fb361141f3993320ba5c127ff313dc0943e261346ae9ea19f4a"
const HISTORICAL_SIDIS_PILOT_FALSE_FLAGS=(HISTORICAL_SIDIS_FALSE_FLAGS...,
    "accepted_prediction_diagnostic","accepted_fitting_gate","all132_coverage_verified")

function historical_sidis_pilot_policy()
    bytes2hex(sha256(read(HISTORICAL_SIDIS_PILOT_POLICY)))==HISTORICAL_SIDIS_PILOT_SHA || error("changed frozen SIDIS pilot policy")
    policy=JSON.parsefile(HISTORICAL_SIDIS_PILOT_POLICY);contract=JSON.parsefile(HISTORICAL_SIDIS_CONTRACT)
    policy["scope"]==HISTORICAL_SIDIS_PILOT_SCOPE && all(policy[k]===false for k in HISTORICAL_SIDIS_PILOT_FALSE_FLAGS) ||
        error("resource pilot cannot approve a diagnostic or fitting")
    policy["parent_contract_sha256"]==bytes2hex(sha256(read(HISTORICAL_SIDIS_CONTRACT))) || error("changed SIDIS pilot parent")
    policy["numerics"]["native_level"]==4 && policy["numerics"]["setting"]==last(contract["numerical_policy"]["refinements"]) &&
        policy["numerics"]["physical_value_max_scaled"]==contract["numerical_policy"]["physical_value_max_scaled"] &&
        policy["numerics"]["relative_floor"]==contract["numerical_policy"]["relative_floor"] || error("changed inherited SIDIS precision")
    policy
end

"Canonical ASCII path/hash dictionary, matching Python provenance.digest."
function sidis_pilot_source_fingerprint(sources)
    all(isascii,keys(sources)) || error("unsupported non-ASCII source path in canonical binding")
    canonical="{"*join([JSON.json(k)*":"*JSON.json(sources[k]) for k in sort!(collect(keys(sources)))],",")*"}"
    bytes2hex(sha256(codeunits(canonical)))
end

"All computational files, including workflow/configuration; no raw cache or run output."
function sidis_pilot_current_sources(root)
    paths=String[]
    skipped_names=("github-state.json","project-tasks.json","actions.lock.json")
    skipped_suffixes=(".so",".dll",".dylib",".o",".mod")
    for directory in ("src","numerics/src","numerics/bin","fortran","bridges","workflow","config","fixtures")
        base=joinpath(root,directory);isdir(base) || continue
        for (folder,dirs,files) in walkdir(base)
            filter!(d->d!="__pycache__" && !endswith(d,".egg-info"),dirs)
            for file in files
                file in skipped_names && continue
                splitext(file)[2] in skipped_suffixes && continue
                push!(paths,joinpath(folder,file))
            end
        end
    end
    for name in ("pyproject.toml","requirements.lock","numerics/Project.toml","numerics/Manifest.toml",
                 "references/source-inventory.json","references/code-review.json")
        isfile(joinpath(root,name)) && push!(paths,joinpath(root,name))
    end
    Dict(replace(relpath(p,root),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in unique(paths))
end

function sidis_pilot_case(request,policy,case_id)
    selected=filter(c->c["id"]==case_id,policy["selection"]["cases"])
    length(selected)==1 || error("case is outside frozen SIDIS pilot selection")
    c=only(selected);points=filter(p->p["id"]==c["point_id"],request.data["points"])
    length(request.data["points"])==132 && length(points)==1 || error("full source request or identified pilot point missing")
    point=only(points)
    point["publication"]==c["publication"] && point["target"]==c["target"] || error("pilot source identity changed")
    point
end

"Separate explicit same-axis source factor; conditional mappings stay source-bound."
function sidis_pilot_source_factor(point)
    y=point["kinematics"]["y"]
    point["publication"]=="1006.4221v1" && return 2(1-y)/(1+(1-y)^2)
    point["publication"] in ("1408.4405","0802.2160") || error("unsupported pilot publication")
    -1.
end

"Saved arithmetic only; this neither reevolves rows nor validates AD accuracy."
function sidis_pilot_value_checks(actual,reference,point,policy)
    floor=policy["numerics"]["relative_floor"];tolerance=policy["numerics"]["physical_value_max_scaled"]
    checks=Dict{String,Any}[]
    for key in ("UU","FUT","amplitude")
        a,b=actual[key],reference[key]
        sidis_source_number(a) && sidis_source_number(b) || error("nonfinite SIDIS pilot operand")
        residual=abs(a-b)/max(abs(b),floor)
        push!(checks,Dict("name"=>"independent scalar "*key,"error"=>residual,"tolerance"=>tolerance,"passed"=>residual<=tolerance))
    end
    actual["UU"]>0 && reference["UU"]>0 || error("nonpositive SIDIS pilot denominator")
    algebra=sidis_pilot_source_factor(point)*actual["FUT"]/actual["UU"]
    residual=abs(actual["amplitude"]-algebra)/max(abs(algebra),floor)
    contract=JSON.parsefile(HISTORICAL_SIDIS_CONTRACT);tol=contract["numerical_policy"]["observable_algebra_max_scaled"]
    push!(checks,Dict("name"=>"source factor and ratio","error"=>residual,"tolerance"=>tol,"passed"=>residual<=tol))
    checks
end

function sidis_pilot_measure(call)
    result=nothing
    seconds=@elapsed result=call()
    result,seconds
end

"Julia-only evaluation samples on prepared values; no derivative-accuracy claim."
function sidis_pilot_evaluations(value,theta,direction,repetitions)
    length(theta)==length(direction)==13 && all(sidis_source_number,theta) && all(sidis_source_number,direction) ||
        error("finite all13 pilot coordinates/direction required")
    repetitions isa Integer && !(repetitions isa Bool) && repetitions>0 || error("positive fixed repetition count required")
    gradient(t)=ForwardDiff.gradient(value,t)
    hvp(t)=ForwardDiff.hessian(value,t)*direction
    first_value,tv=sidis_pilot_measure(()->value(theta))
    first_gradient,tg=sidis_pilot_measure(()->gradient(theta))
    first_hvp,th=sidis_pilot_measure(()->hvp(theta))
    warm=Dict{String,Any}[]
    for i in 1:repetitions
        v,sv=sidis_pilot_measure(()->value(theta));g,sg=sidis_pilot_measure(()->gradient(theta));h,sh=sidis_pilot_measure(()->hvp(theta))
        push!(warm,Dict("index"=>i,"value"=>v,"gradient"=>g,"HVP"=>h,
            "seconds"=>Dict("value"=>sv,"gradient"=>sg,"HVP"=>sh)))
    end
    finite=length(first_gradient)==length(first_hvp)==13 && all(isfinite,first_gradient) && all(isfinite,first_hvp)
    repeatable=all(r["value"]==first_value && r["gradient"]==first_gradient && r["HVP"]==first_hvp for r in warm)
    samples=Dict("gradient"=>first_gradient,"HVP"=>first_hvp,"direction"=>direction,
        "independent_derivative_reference_executed"=>false,"derivative_accuracy_accepted"=>false)
    times=Dict("first_call_seconds"=>Dict("value"=>tv,"gradient"=>tg,"HVP"=>th),
        "warm_repetitions"=>warm,"scope"=>"informational timing samples including first-call compilation; competing workloads not controlled")
    (;amplitude=first_value,samples,times,finite,repeatable)
end
