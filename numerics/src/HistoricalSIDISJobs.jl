const HISTORICAL_SIDIS_JOB_POLICY=normpath(joinpath(@__DIR__,"../../config/historical-sidis-job-policy.json"))
const HISTORICAL_SIDIS_JOB_POLICY_SHA="7ce5d19f131c573caceaf64d24b2cc09fe589aaa4b628137d73d2b56df995911"
const HISTORICAL_SIDIS_JOB_SCOPE="historical_reconstructed_sidis_complete_point_jobs"
const SIDIS_JOB_FALSE_FLAGS=("production_ready","physical_extraction","historical_fit_ready",
    "experimental_data_admitted","fit_performed","accepted_prediction_diagnostic",
    "accepted_fitting_gate","G2PERF_pass","interval_coverage_verified","all132_coverage_verified")

function historical_sidis_job_policy()
    bytes2hex(sha256(read(HISTORICAL_SIDIS_JOB_POLICY)))==HISTORICAL_SIDIS_JOB_POLICY_SHA || error("changed frozen full SIDIS job policy")
    policy=JSON.parsefile(HISTORICAL_SIDIS_JOB_POLICY)
    policy["scope"]==HISTORICAL_SIDIS_JOB_SCOPE && all(policy[k]===false for k in SIDIS_JOB_FALSE_FLAGS) || error("job policy cannot self-admit")
    bytes2hex(sha256(read(HISTORICAL_SIDIS_CONTRACT)))==policy["parent_contract_sha256"] || error("changed full SIDIS scientific contract")
    policy["execution"]["source_points"]==132 && policy["execution"]["maximum_parallel_jobs_per_checkout"]==1 || error("changed point-job coverage or concurrency")
    policy
end

function sidis_job_point_key(index)
    index isa Integer && !(index isa Bool) && 1<=index<=132 || error("full SIDIS point index must be1..132")
    "p"*lpad(string(index),4,'0')
end

function sidis_job_scaled(actual,reference,floor)
    length(actual)==length(reference)>0 && isfinite(floor) && floor>0 || error("incomplete comparison operands")
    all(v->v isa Real && !(v isa Bool) && isfinite(v),actual) &&
        all(v->v isa Real && !(v isa Bool) && isfinite(v),reference) || error("nonfinite comparison operand")
    maximum(abs(a-b)/max(abs(b),floor) for (a,b) in zip(actual,reference))
end

function sidis_job_decimal(value)
    value isa Real && !(value isa Bool) && isfinite(value) || error("finite decimal reference required")
    string(value)
end

function sidis_job_parse_decimal(value)
    value isa AbstractString && ncodeunits(value)<=512 || error("bounded decimal string required")
    result=parse(BigFloat,value)
    isfinite(result) || error("nonfinite decimal reference")
    result
end
sidis_job_decimal_vector(values)=sidis_job_decimal.(values)
sidis_job_parse_vector(values)=sidis_job_parse_decimal.(values)

"Fourth-order inward derivative, used only with the frozen five-sample policy."
function sidis_job_five_point(samples,h)
    length(samples)==5 && h isa Real && !(h isa Bool) && isfinite(h) && h>0 || error("five samples and a positive step required")
    sum(c*v for (c,v) in zip((-25,48,-36,16,-3),samples))/(12h)
end

"Complete Julia gradient and directional derivative of its full13 Jacobian."
function sidis_job_derivatives(value,theta,direction)
    length(theta)==length(direction)==13 && all(isfinite,theta) && all(isfinite,direction) || error("finite full13 coordinates and direction required")
    gradient=ForwardDiff.gradient(value,theta)
    hvp=vec(ForwardDiff.jacobian(s->ForwardDiff.gradient(value,theta+s[1]*direction),[zero(eltype(theta))]))
    all(isfinite,gradient) && all(isfinite,hvp) || error("nonfinite full13 gradient/HVP")
    (;gradient,HVP=hvp)
end

"Full decimal operands for one all13 fixed-center inward refinement; no admission."
function sidis_job_reference_step(value,jacobian,theta,direction,level,policy;
        base_value,base_jacobian,reference_HVP,candidate_gradient,candidate_HVP,previous=nothing)
    precision(BigFloat)==policy["reference_bits"]==256 || error("fixed256-bit reference arithmetic required")
    length(theta)==length(direction)==length(base_jacobian)==length(reference_HVP)==
        length(candidate_gradient)==length(candidate_HVP)==13 || error("all13 reference coordinates required")
    level isa Integer && !(level isa Bool) && 1<=level<=4 || error("invalid inward refinement level")
    h=BigFloat(policy["steps"][level]);t=BigFloat.(theta);d=BigFloat.(direction)
    coordinate_values=Vector{Vector{BigFloat}}()
    for coordinate in 1:13
        unit=zeros(BigFloat,13);unit[coordinate]=1
        push!(coordinate_values,vcat([BigFloat(base_value)],[BigFloat(value(t+m*h*unit)) for m in 1:4]))
    end
    jacobian_values=vcat([BigFloat.(base_jacobian)],[BigFloat.(jacobian(t+m*h*d)) for m in 1:4])
    all(length(row)==13 for row in jacobian_values) || error("incomplete directional Jacobian sample")
    gradient=[sidis_job_five_point(row,h) for row in coordinate_values]
    hvp=[sidis_job_five_point([row[j] for row in jacobian_values],h) for j in 1:13]
    floor=BigFloat(policy["relative_floor"])
    errors=Dict{String,Any}(
        "candidate_gradient_vs_inward"=>sidis_job_scaled(candidate_gradient,gradient,floor),
        "candidate_HVP_vs_inward"=>sidis_job_scaled(candidate_HVP,hvp,floor),
        "reference_gradient_vs_inward"=>sidis_job_scaled(base_jacobian,gradient,floor),
        "reference_HVP_vs_inward"=>sidis_job_scaled(reference_HVP,hvp,floor),
        "adjacent_gradient"=>nothing,"adjacent_HVP"=>nothing)
    if previous!==nothing
        previous["level"]==level-1 && previous["declared_step"]==policy["steps"][level-1] || error("missing immediate adjacent inward reference")
        errors["adjacent_gradient"]=sidis_job_scaled(gradient,sidis_job_parse_vector(previous["inward_gradient"]),floor)
        errors["adjacent_HVP"]=sidis_job_scaled(hvp,sidis_job_parse_vector(previous["inward_HVP"]),floor)
    elseif level!=1
        error("later inward refinements require their original adjacent operands")
    end
    Dict{String,Any}("level"=>level,"declared_step"=>policy["steps"][level],"step_decimal"=>sidis_job_decimal(h),
        "reference_bits"=>256,"coordinate_value_samples"=>sidis_job_decimal_vector.(coordinate_values),
        "directional_Jacobian_samples"=>sidis_job_decimal_vector.(jacobian_values),
        "inward_gradient"=>sidis_job_decimal_vector(gradient),"inward_HVP"=>sidis_job_decimal_vector(hvp),
        "errors_decimal"=>Dict(k=>(v===nothing ? nothing : sidis_job_decimal(v)) for (k,v) in errors),
        "adjacent_level"=>previous===nothing ? nothing : previous["level"],
        "scope"=>"fixed_center_reference_operands_not_a_fit_domain_or_production_certificate")
end

function sidis_job_ordered_packets(point_ids,packets,phase)
    length(point_ids)==length(packets)==132 && length(unique(point_ids))==132 &&
        all(id->id isa AbstractString && !isempty(id),point_ids) || error("exact full132 source coverage required")
    all(packet["point"]["id"]==id && packet["phase"]==phase for (packet,id) in zip(packets,point_ids)) ||
        error("missing, reordered, duplicate or wrong-phase SIDIS point packet")
    nothing
end

function sidis_job_native_values(row,level)
    row["level"] isa Integer && !(row["level"] isa Bool) && row["level"]==level && row["bound_inverse_GeV"]==24 || error("mixed native level or bound")
    values=[row["values"][key] for key in ("UU","FUT","amplitude")]
    all(v->v isa Real && !(v isa Bool) && isfinite(v),values) && values[1]>0 || error("invalid native point values")
    values
end

"Choose one common native level from all132 already source-verified point packets."
function sidis_job_native_decision(point_ids,baseline;fourth=nothing)
    sidis_job_ordered_packets(point_ids,baseline,"native123")
    contract=JSON.parsefile(HISTORICAL_SIDIS_CONTRACT);policy=contract["numerical_policy"]
    floor=policy["relative_floor"];tolerance=policy["physical_value_max_scaled"]
    comparisons=Dict{String,Any}[];maxima=Dict{String,Any}()
    all(length(p["levels"])==3 for p in baseline) || error("all first3 native levels required")
    for level in 2:3
        errors=Float64[]
        for (id,packet) in zip(point_ids,baseline)
            before=sidis_job_native_values(packet["levels"][level-1],level-1)
            after=sidis_job_native_values(packet["levels"][level],level)
            residual=sidis_job_scaled(after,before,floor)
            push!(errors,residual);push!(comparisons,Dict("point_id"=>id,"level"=>level,"max_scaled"=>residual))
        end
        maxima[string(level)]=maximum(errors)
    end
    selected=maxima["3"]<=tolerance ? 3 : 0
    status=selected==3 ? "COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" : "NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"
    if fourth!==nothing
        selected==0 || error("level4 must not replace an already selected common level3")
        sidis_job_ordered_packets(point_ids,fourth,"native4");errors=Float64[]
        for (id,packet,coarse) in zip(point_ids,fourth,baseline)
            length(packet["levels"])==1 || error("one fourth native level required")
            before=sidis_job_native_values(coarse["levels"][3],3)
            after=sidis_job_native_values(packet["levels"][1],4)
            residual=sidis_job_scaled(after,before,floor)
            push!(errors,residual);push!(comparisons,Dict("point_id"=>id,"level"=>4,"max_scaled"=>residual))
        end
        maxima["4"]=maximum(errors);selected=maxima["4"]<=tolerance ? 4 : 0
        status=selected==4 ? "COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" : "FAIL"
    end
    Dict{String,Any}("status"=>status,"selected_level"=>selected,"source_ids"=>copy(point_ids),
        "source_points"=>132,"maximum_scaled_by_level"=>maxima,"comparisons"=>comparisons,
        "unchanged_target"=>tolerance,"minimum_level"=>policy["minimum_refinement_level"],
        (k=>false for k in SIDIS_JOB_FALSE_FLAGS)...)
end

function sidis_job_reference_metrics(packet,level,policy)
    rows=packet["reference_steps"]
    all(r->r["level"] isa Integer && !(r["level"] isa Bool),rows) || error("integer derivative levels required")
    expected=packet["phase"]=="tailref123" ? [1,2,3] : [4]
    [r["level"] for r in rows]==expected || error("incomplete or reordered inward ladder")
    row=only(r for r in rows if r["level"]==level)
    length(row["inward_gradient"])==length(row["inward_HVP"])==13 &&
        length(row["coordinate_value_samples"])==13 && all(r->length(r)==5,row["coordinate_value_samples"]) &&
        length(row["directional_Jacobian_samples"])==5 && all(r->length(r)==13,row["directional_Jacobian_samples"]) || error("missing all13 retained derivative operands")
    row["declared_step"]==policy["steps"][level] && row["reference_bits"]==policy["reference_bits"] &&
        row["adjacent_level"]==level-1 || error("wrong common derivative step/precision/adjacent level")
    Set(keys(packet["centre_errors_decimal"]))==Set(("gradient","HVP")) || error("missing independent center comparisons")
    Set(keys(row["errors_decimal"]))==Set(("candidate_gradient_vs_inward","candidate_HVP_vs_inward",
        "reference_gradient_vs_inward","reference_HVP_vs_inward","adjacent_gradient","adjacent_HVP")) || error("incomplete all13 reference comparisons")
    values=[sidis_job_parse_decimal(v) for v in vcat(collect(Base.values(packet["centre_errors_decimal"])),collect(Base.values(row["errors_decimal"])))]
    all(>=(0),values) || error("negative derivative error")
    maximum(values)
end

"Global derivative selection requires the complete source order and a common native level."
function sidis_job_derivative_decision(point_ids,initial,native_level;fourth=nothing)
    native_level isa Integer && !(native_level isa Bool) && native_level in (3,4) || error("one selected common native level required")
    sidis_job_ordered_packets(point_ids,initial,"tailref123")
    all(p["native_level"]==native_level for p in initial) || error("mixed native levels in all132 derivative records")
    policy=JSON.parsefile(HISTORICAL_SIDIS_CONTRACT)["derivative_policy"]
    setprecision(policy["reference_bits"]) do
        target=BigFloat(policy["gradient_max_scaled"])
        target==BigFloat(policy["HVP_max_scaled"]) || error("unexpected derivative targets")
        errors=[sidis_job_reference_metrics(p,3,policy) for p in initial]
        maxima=Dict{String,Any}("3"=>sidis_job_decimal(maximum(errors)))
        selected=maximum(errors)<=target ? 3 : 0
        status=selected==3 ? "COMMON_DERIVATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" : "DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"
        if fourth!==nothing
            selected==0 || error("fourth derivative step must not replace accepted common step3")
            sidis_job_ordered_packets(point_ids,fourth,"reference4")
            all(p["native_level"]==native_level for p in fourth) || error("mixed native levels in fourth-step records")
            errors=[sidis_job_reference_metrics(p,4,policy) for p in fourth]
            maxima["4"]=sidis_job_decimal(maximum(errors));selected=maximum(errors)<=target ? 4 : 0
            status=selected==4 ? "COMMON_DERIVATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" : "FAIL"
        end
        Dict{String,Any}("status"=>status,"selected_reference_level"=>selected,"native_level"=>native_level,
            "source_ids"=>copy(point_ids),"source_points"=>132,"coordinates"=>13,
            "maximum_scaled_by_level_decimal"=>maxima,"unchanged_target"=>policy["gradient_max_scaled"],
            (k=>false for k in SIDIS_JOB_FALSE_FLAGS)...)
    end
end
