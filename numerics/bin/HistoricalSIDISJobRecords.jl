module HistoricalSIDISJobRecords
using CollinsFit, JSON, SHA, Dates
include("HistoricalSIDISSourceMath.jl")

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const RECORDED="POINT_JOB_RECORDED_NOT_DIAGNOSTIC_GATE"
const AUDITED="SAVED_JOB_ARITHMETIC_REVIEWED_NOT_PRODUCTION_GATE"
const POINT_PHASES=("native123","native4","tailref123","reference4")
const GLOBAL_PHASES=("native-select3","native-selection","derivative-select3","source-math","diagnostic")
const FLAGS=CollinsFit.SIDIS_JOB_FALSE_FLAGS
digest(path)=bytes2hex(sha256(read(path)))
false_flags()=Dict{String,Any}(k=>false for k in FLAGS)
point_folder(run,phase,index)=joinpath(run,phase,sidis_job_point_key(index))
point_report(run,phase,index)=joinpath(point_folder(run,phase,index),"report.json")
global_report(run,phase)=joinpath(run,phase,"report.json")
all_ids(ctx)=[p["id"] for p in ctx.request.data["points"]]

function bindings(requestfile,nativefile,inputfile)
    policy=historical_sidis_job_policy();request=load_historical_sidis_point_request(requestfile)
    sources=sidis_pilot_current_sources(ROOT);fp=sidis_pilot_source_fingerprint(sources)
    request.data["source_fingerprint"]==fp || error("stale complete SIDIS request")
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    any(basename(p)=="hoppet_bridge.f90" && h==digest(joinpath(ROOT,"fortran/hoppet_bridge.f90")) for (p,h) in native["source_hashes"]) || error("native bridge differs from current source")
    any(basename(p)=="lhapdf.cpp" && h==digest(joinpath(ROOT,"bridges/lhapdf.cpp")) for (p,h) in inputs["source_hashes"]) || error("parton bridge differs from current source")
    (;request,policy,sources,fp,native,inputs,requestfile=abspath(requestfile),nativefile=abspath(nativefile),inputfile=abspath(inputfile),
        contract=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT))
end

function stable(ctx)
    historical_sidis_request_valid(ctx.request)
    sidis_pilot_current_sources(ROOT)==ctx.sources &&
        CollinsFit.verified_hoppet_build(ctx.nativefile)==ctx.native &&
        CollinsFit.verify_parton_build(ctx.inputfile)==ctx.inputs || error("SIDIS code/build inputs changed during complete job")
    true
end

function header(ctx,phase;index=nothing)
    phase in (POINT_PHASES...,GLOBAL_PHASES...) || error("unknown complete SIDIS phase")
    record=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "theory_scope"=>"reconstructed_theory","observation_scope"=>"historical_benchmark_candidate",
        "phase"=>phase,"source_fingerprint"=>ctx.fp,"sources"=>ctx.sources,
        "scientific_scope"=>ctx.contract["scope"],
        "request_sha256"=>ctx.request.sha256,"source_points"=>132,
        "policy_sha256"=>digest(CollinsFit.HISTORICAL_SIDIS_JOB_POLICY),
        "scientific_contract_sha256"=>ctx.request.data["contract_sha256"],
        "parameters"=>deepcopy(ctx.request.data["parameters"]),
        "native_build_receipt"=>ctx.nativefile,"parton_build_receipt"=>ctx.inputfile,
        "native_build_scope"=>ctx.native["status"],"started_at"=>string(now(UTC)),
        "limitations"=>vcat(ctx.contract["limitations"],[ctx.policy["resource_basis"]["unmeasured_costs"]]),false_flags()...)
    if index!==nothing
        sidis_job_point_key(index);record["point_index"]=index
        record["point"]=deepcopy(ctx.request.data["points"][index])
    end
    record
end

function selection(ctx,run,phase)
    record=JSON.parsefile(global_report(run,phase))
    validate_header(record,ctx,phase)
    record["source_ids"]==all_ids(ctx) || error("selection lacks exact full132 source coverage")
    all(isfile(p) && digest(p)==h for (p,h) in record["dependency_sha256"]) || error("changed selection input bindings")
    audit=JSON.parsefile(joinpath(run,phase,"audit.json"))
    audit["status"]==AUDITED && audit["source_fingerprint"]==ctx.fp &&
        audit["report_sha256"]==digest(global_report(run,phase)) || error("missing or stale independent selection audit")
    record
end

function point_dependencies(ctx,run,phase,index)
    paths=[ctx.requestfile,ctx.nativefile,ctx.inputfile]
    if phase=="native4"
        push!(paths,global_report(run,"native-select3"),point_report(run,"native123",index))
    elseif phase in ("tailref123","reference4")
        chosen=selection(ctx,run,"native-selection")
        chosen["status"]=="COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" && chosen["selected_level"] in (3,4) || error("native level is not globally selected")
        push!(paths,global_report(run,"native-selection"),
            point_report(run,chosen["selected_level"]==3 ? "native123" : "native4",index))
        if phase=="reference4"
            push!(paths,global_report(run,"derivative-select3"),point_report(run,"tailref123",index))
        end
    elseif phase!="native123"
        error("unknown point phase")
    end
    unique(paths)
end

function validate_header(record,ctx,phase;index=nothing,final=false)
    record["scope"]==(final ? ctx.contract["scope"] : CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE) &&
        record["scientific_scope"]==ctx.contract["scope"] && record["phase"]==phase &&
        record["theory_scope"]=="reconstructed_theory" && record["observation_scope"]=="historical_benchmark_candidate" || error("wrong SIDIS report scope")
    for flag in FLAGS
        expected=final && flag in ("accepted_prediction_diagnostic","all132_coverage_verified")
        get(record,flag,nothing)===expected || error("unexpected SIDIS admission flag: "*flag)
    end
    record["source_fingerprint"]==ctx.fp && record["sources"]==ctx.sources &&
        record["request_sha256"]==ctx.request.sha256 && record["source_points"]==132 &&
        record["policy_sha256"]==digest(CollinsFit.HISTORICAL_SIDIS_JOB_POLICY) &&
        record["scientific_contract_sha256"]==ctx.request.data["contract_sha256"] &&
        record["parameters"]==ctx.request.data["parameters"] &&
        record["native_build_receipt"]==ctx.nativefile && record["parton_build_receipt"]==ctx.inputfile || error("changed SIDIS source/code/parameter bindings")
    if index!==nothing
        record["point_index"]==index && record["point"]==ctx.request.data["points"][index] || error("wrong exact SIDIS point identity")
    end
    true
end

function audit_resource(resource,point,bound,ctx)
    all(resource[k] isa Integer && !(resource[k] isa Bool) && resource[k]>0 for k in ("radial_nodes","row_coefficient_payload_bytes")) || error("invalid producer payload observation")
    seconds=resource["cold_preparation_seconds_informational"]
    seconds isa Real && !(seconds isa Bool) && isfinite(seconds) && seconds>=0 || error("invalid preparation timing observation")
    provenance=resource["context_provenance"];target=ctx.request.data["target_weights"][point["target"]]
    provenance["scope"]=="reconstructed_theory" && provenance["production_ready"]===provenance["physical_extraction"]===false &&
        provenance["process"]=="SIDIS" && provenance["range_policy"]=="reference" &&
        provenance["radial_bound_inverse_GeV"]==bound && provenance["target"]["name"]==point["target"] &&
        provenance["target"]["unpolarized"]==target["unpolarized"] && provenance["target"]["polarized"]==target["polarized"] || error("changed target/range/bound context provenance")
    true
end

function resource_text(path)
    lines=filter(l->startswith(l,"COLLINSFLOW_TIME_V1\t"),readlines(path))
    length(lines)==1 || error("missing or ambiguous GNU time measurement")
    fields=split(only(lines),'\t');length(fields)==6 || error("incomplete GNU time measurement")
    result=Dict{String,Any}("maximum_RSS_KiB"=>parse(Int,fields[2]),"wall_seconds"=>parse(Float64,fields[3]),
        "user_seconds"=>parse(Float64,fields[4]),"system_seconds"=>parse(Float64,fields[5]),
        "child_exit_code"=>parse(Int,fields[6]),"scope"=>"GNU_time_direct_child_maximum_not_simultaneous_process_tree_peak")
    result["maximum_RSS_KiB"]>0 && result["child_exit_code"]==0 &&
        all(isfinite(result[k]) && result[k]>=0 for k in ("wall_seconds","user_seconds","system_seconds")) || error("invalid or failed GNU time measurement")
    result
end

function invocation(ctx,folder,phase;index=nothing,final=false)
    inv=JSON.parsefile(joinpath(folder,"invocation.json"))
    inv["status"]=="COMPLETE_JOB_RECORDED_NOT_PRODUCTION_GATE" && inv["phase"]==phase &&
        inv["source_fingerprint"]==ctx.fp && inv["policy_sha256"]==digest(CollinsFit.HISTORICAL_SIDIS_JOB_POLICY) || error("incomplete or stale SIDIS job invocation")
    index===nothing || inv["point_index"]==index || error("wrong invoked source point")
    required=Set((ctx.requestfile,ctx.nativefile,ctx.inputfile,CollinsFit.HISTORICAL_SIDIS_JOB_POLICY,"/usr/bin/time",joinpath(ROOT,"scripts/run_historical_sidis_job.py")))
    issubset(required,Set(keys(inv["input_hashes"]))) && inv["time_binary"]=="/usr/bin/time" || error("missing complete invocation input bindings")
    for flag in FLAGS
        expected=final && flag in ("accepted_prediction_diagnostic","all132_coverage_verified")
        inv[flag]===expected || error("self-admitting SIDIS invocation")
    end
    [s["name"] for s in inv["steps"]]==["producer","audit"] || error("incomplete whole-job command ledger")
    resource_policy=ctx.policy["resources"][phase in POINT_PHASES ? phase : "global"]
    for step in inv["steps"]
        step["returncode"]==step["process_returncode"]==0 &&
            step["interruption"]===step["termination_signal"]===nothing || error("interrupted command cannot approve a point")
        name=step["name"]
        digest(joinpath(folder,name*".log"))==step["log_sha256"] &&
            digest(joinpath(folder,name*".resources.txt"))==step["resource_sha256"] || error("changed command evidence")
        measured=resource_text(joinpath(folder,name*".resources.txt"))
        measured==step["resource"] && measured["maximum_RSS_KiB"]<=resource_policy["maximum_command_RSS_KiB"] || error("failed or changed resource budget")
    end
    sum(s["wrapper_elapsed_seconds"] for s in inv["steps"])<=resource_policy["whole_Julia_jobs_time_budget_seconds"] || error("whole-job budget exhausted")
    all(isfile(p) && digest(p)==h for (p,h) in inv["input_hashes"]) &&
        all(isfile(p) && digest(p)==h for (p,h) in inv["outputs"]) || error("changed complete-job inputs or outputs")
    expected_outputs=["report.json","audit.json","producer.log","producer.resources.txt","audit.log","audit.resources.txt"]
    final && append!(expected_outputs,["point-comparison.csv","comparison-hermes.png","comparison-compass-proton.png","comparison-compass-deuteron.png"])
    Set(keys(inv["outputs"]))==Set(joinpath(folder,name) for name in expected_outputs) || error("incomplete job output inventory")
    initial=joinpath(folder,"initial-invocation.json")
    digest(initial)==inv["initial_record_sha256"] &&
        JSON.parsefile(initial)["record_role"]=="initial_sentinel_not_an_executed_failure" || error("changed initial sentinel")
    isfile(joinpath(folder,"job.completed")) || error("missing scheduling completion marker")
    inv
end

function scalar_checks(values,reference,point,ctx)
    floor=ctx.contract["numerical_policy"]["relative_floor"]
    tolerance=ctx.contract["numerical_policy"]["physical_value_max_scaled"]
    checks=Dict{String,Any}[]
    for key in ("UU","FUT","amplitude")
        a,b=values[key],reference[key]
        residual=abs(a-b)/max(abs(b),floor)
        isfinite(a) && isfinite(b) || error("nonfinite scalar reference")
        push!(checks,Dict("name"=>"independent scalar "*key,"error"=>residual,"tolerance"=>tolerance,"passed"=>residual<=tolerance))
    end
    values["UU"]>0 && reference["UU"]>0 || error("nonpositive SIDIS denominator")
    checks
end

function compare_decimal(actual,expected)
    a=CollinsFit.sidis_job_parse_decimal(actual)
    # This scaled256-bit saved-arithmetic fidelity bound is separate from the
    # physical1e-4/1e-7 targets. All full decimal samples remain available.
    abs(a-expected)<=big"1e-60"*max(one(a),abs(expected)) || error("changed cancellation-sensitive decimal arithmetic")
    true
end

"Join saved scalar, candidate and reference-center representations without new physics."
function audit_value_links(packet,ctx)
    value=packet["values"];factor=packet["source_factor"];policy=ctx.contract["numerical_policy"]
    checked=Dict{String,Any}[]
    if haskey(packet,"tails")
        packet["tails"][3]["values"]==value || error("final tail and comparison values differ")
        append!(checked,[row["values"] for row in packet["tails"]])
    else
        push!(checked,value)
    end
    for values in checked
        all(v->v isa Real && !(v isa Bool) && isfinite(v),Base.values(values)) && values["UU"]>0 || error("invalid saved source structure values")
        algebra=factor*values["FUT"]/values["UU"]
        abs(values["amplitude"]-algebra)/max(abs(algebra),policy["relative_floor"])<=policy["observable_algebra_max_scaled"] || error("saved source amplitude does not match signed structure ratio")
    end
    packet["candidate_derivatives"]["value"]==value["amplitude"] || error("candidate value and final source amplitude differ")
    checks=scalar_checks(value,packet["independent_scalar_reference"],nothing,ctx)
    checks==packet["scalar_checks"] && all(c["passed"] for c in checks) || error("changed independent scalar composition")
    setprecision(256) do
        centre=CollinsFit.sidis_job_parse_decimal(packet["reference_centre"]["value"])
        scalar=BigFloat(packet["independent_scalar_reference"]["amplitude"])
        abs(centre-scalar)/max(abs(scalar),BigFloat(policy["relative_floor"]))<=BigFloat(policy["physical_value_max_scaled"]) || error("scalar and high-precision reference centers differ")
    end
    true
end

"Independent reader formula; never calls the producer's stencil implementation."
function recompute_reference_step(row,packet,ctx;previous=nothing)
    policy=ctx.contract["derivative_policy"]
    setprecision(policy["reference_bits"]) do
        level=row["level"];level isa Integer && !(level isa Bool) && 1<=level<=4 || error("invalid stored reference level")
        row["declared_step"]==policy["steps"][level] && row["reference_bits"]==256 || error("changed inward step")
        h=CollinsFit.sidis_job_parse_decimal(row["step_decimal"])
        h==BigFloat(policy["steps"][level]) || error("changed decimal step")
        v=row["coordinate_value_samples"];j=row["directional_Jacobian_samples"]
        length(v)==13 && all(length(r)==5 for r in v) && length(j)==5 && all(length(r)==13 for r in j) || error("missing all13 decimal samples")
        scalar(r)=begin
            a,b,c,d,e=CollinsFit.sidis_job_parse_vector(r)
            (-25a+48b-36c+16d-3e)/(12h)
        end
        gradient=scalar.(v);hvp=[scalar([r[k] for r in j]) for k in 1:13]
        all(compare_decimal(a,b) for (a,b) in zip(row["inward_gradient"],gradient)) &&
            all(compare_decimal(a,b) for (a,b) in zip(row["inward_HVP"],hvp)) || error("changed all13 reconstructed derivatives")
        reference=packet["reference_centre"];candidate=packet["candidate_derivatives"]
        centre_j=CollinsFit.sidis_job_parse_vector(reference["gradient"])
        centre_h=CollinsFit.sidis_job_parse_vector(reference["HVP"])
        length(centre_j)==length(centre_h)==length(candidate["gradient"])==length(candidate["HVP"])==13 || error("missing full center derivatives")
        all(v->v isa Real&&isfinite(v)&&!(v isa Bool),candidate["gradient"]) &&
            all(v->v isa Real&&isfinite(v)&&!(v isa Bool),candidate["HVP"]) || error("nonfinite candidate derivative")
        for samples in v
            compare_decimal(samples[1],CollinsFit.sidis_job_parse_decimal(reference["value"]))
        end
        all(compare_decimal(a,b) for (a,b) in zip(j[1],centre_j)) || error("directional reference center changed")
        floor=BigFloat(policy["relative_floor"])
        residual(a,b)=maximum(abs(BigFloat(x)-BigFloat(y))/max(abs(BigFloat(y)),floor) for (x,y) in zip(a,b))
        errors=Dict{String,Any}("candidate_gradient_vs_inward"=>residual(candidate["gradient"],gradient),
            "candidate_HVP_vs_inward"=>residual(candidate["HVP"],hvp),
            "reference_gradient_vs_inward"=>residual(centre_j,gradient),"reference_HVP_vs_inward"=>residual(centre_h,hvp),
            "adjacent_gradient"=>nothing,"adjacent_HVP"=>nothing)
        if previous===nothing
            level==1 && row["adjacent_level"]===nothing || error("missing adjacent reference packet")
        else
            row["adjacent_level"]==previous["level"]==level-1 || error("wrong adjacent derivative level")
            errors["adjacent_gradient"]=residual(gradient,CollinsFit.sidis_job_parse_vector(previous["inward_gradient"]))
            errors["adjacent_HVP"]=residual(hvp,CollinsFit.sidis_job_parse_vector(previous["inward_HVP"]))
        end
        Set(keys(errors))==Set(keys(row["errors_decimal"])) || error("missing reference comparison metric")
        for (key,value) in errors
            value===nothing ? row["errors_decimal"][key]===nothing || error("invented adjacent comparison") : compare_decimal(row["errors_decimal"][key],value)
        end
        compare_decimal(packet["centre_errors_decimal"]["gradient"],residual(candidate["gradient"],centre_j))
        compare_decimal(packet["centre_errors_decimal"]["HVP"],residual(candidate["HVP"],centre_h))
    end
    true
end

function audit_point(ctx,run,phase,index;complete_invocation=false)
    path=point_report(run,phase,index);packet=JSON.parsefile(path)
    validate_header(packet,ctx,phase;index)
    packet["status"]==RECORDED || error("incomplete or failed physical point record")
    expected=Dict(p=>digest(p) for p in point_dependencies(ctx,run,phase,index))
    packet["dependency_sha256"]==expected || error("changed or incomplete point dependencies")
    point=packet["point"];factor=sidis_pilot_source_factor(point)
    if phase in ("native123","native4")
        expected_levels=phase=="native123" ? [1,2,3] : [4]
        [r["level"] for r in packet["levels"]]==expected_levels || error("missing declared native observations")
        for row in packet["levels"]
            values=CollinsFit.sidis_job_native_values(row,row["level"])
            row["setting"]==ctx.contract["numerical_policy"]["refinements"][row["level"]] || error("changed native setting")
            abs(values[3]-factor*values[2]/values[1])/max(abs(factor*values[2]/values[1]),1e-8)<=ctx.contract["numerical_policy"]["observable_algebra_max_scaled"] || error("wrong native source factor/ratio")
            audit_resource(row["resources"],point,24,ctx)
        end
    else
        chosen=selection(ctx,run,"native-selection")
        packet["native_level"]==chosen["selected_level"] || error("point does not use common native level")
        packet["source_factor"]==factor && packet["reference_scope"]["native_rows_shared"]===true &&
            packet["reference_scope"]["initial_functions_shared"]===true && packet["reference_scope"]["independent_evolution"]===false || error("wrong reference scope/source factor")
        packet["candidate_derivatives"]["direction"]==ctx.contract["derivative_policy"]["direction"] || error("changed physical HVP direction")
        setprecision(256) do
            centre=packet["reference_centre"]
            centre["reference_bits"]==256 && CollinsFit.sidis_job_parse_vector(centre["parameters_decimal"])==BigFloat.(ctx.request.data["parameters"]["values"]) &&
                CollinsFit.sidis_job_parse_vector(centre["direction_decimal"])==BigFloat.(ctx.contract["derivative_policy"]["direction"]) || error("changed decimal fixed-center coordinates")
        end
        expected_levels=phase=="tailref123" ? [1,2,3] : [4]
        [r["level"] for r in packet["reference_steps"]]==expected_levels || error("missing full inward ladder")
        previous=phase=="reference4" ? JSON.parsefile(point_report(run,"tailref123",index))["reference_steps"][3] : nothing
        if phase=="reference4"
            old=JSON.parsefile(point_report(run,"tailref123",index))
            packet["kernel_signature"]==old["kernel_signature"] && packet["reference_centre"]==old["reference_centre"] &&
                packet["candidate_derivatives"]==old["candidate_derivatives"] && packet["values"]==old["values"] &&
                packet["independent_scalar_reference"]==old["independent_scalar_reference"] || error("fourth reference does not share the identical fixed-center kernel")
        end
        for row in packet["reference_steps"]
            recompute_reference_step(row,packet,ctx;previous);previous=row
        end
        audit_value_links(packet,ctx)
        if phase=="tailref123"
            [r["bound_inverse_GeV"] for r in packet["tails"]]==[12,24,32] || error("missing finite tail observations")
            fine=packet["tails"][3]["values"];coarse=packet["tails"][2]["values"]
            native_phase=chosen["selected_level"]==3 ? "native123" : "native4"
            native_packet=JSON.parsefile(point_report(run,native_phase,index))
            native_row=only(row for row in native_packet["levels"] if row["level"]==chosen["selected_level"])
            coarse==native_row["values"] && packet["tails"][2]["resources"]==native_row["resources"] || error("tail24 operands differ from the selected native point")
            residual=maximum(abs(fine[k]-coarse[k])/max(abs(coarse[k]),1e-8) for k in ("UU","FUT","amplitude"))
            residual==packet["tail_error"] && residual<=ctx.contract["numerical_policy"]["physical_value_max_scaled"] || error("failed or changed finite tail comparison")
            checks=scalar_checks(packet["values"],packet["independent_scalar_reference"],point,ctx)
            checks==packet["scalar_checks"] && all(c["passed"] for c in checks) || error("failed or changed independent scalar arithmetic")
            all(audit_resource(row["resources"],point,row["bound_inverse_GeV"],ctx) for row in packet["tails"])
        end
        audit_resource(packet["final_preparation_resources"],point,32,ctx)
    end
    if complete_invocation
        inv=invocation(ctx,point_folder(run,phase,index),phase;index)
        audit=JSON.parsefile(joinpath(point_folder(run,phase,index),"audit.json"))
        audit["status"]==AUDITED && audit["source_fingerprint"]==ctx.fp && audit["report_sha256"]==digest(path) &&
            all(audit[k]===false for k in FLAGS) || error("missing or self-admitting point arithmetic audit")
    end
    packet
end

point_packets(ctx,run,phase)=[audit_point(ctx,run,phase,i;complete_invocation=true) for i in 1:132]

function global_packet(ctx,run,phase)
    path=global_report(run,phase);record=JSON.parsefile(path)
    validate_header(record,ctx,phase;final=phase=="diagnostic")
    inv=invocation(ctx,joinpath(run,phase),phase;final=phase=="diagnostic")
    audit=JSON.parsefile(joinpath(run,phase,"audit.json"))
    audit["status"]==AUDITED && audit["report_sha256"]==digest(path) && audit["source_fingerprint"]==ctx.fp || error("missing global arithmetic audit")
    record
end

function global_dependencies(ctx,run,phase)
    paths=[ctx.requestfile,ctx.nativefile,ctx.inputfile]
    add_points(name)=append!(paths,[joinpath(point_folder(run,name,i),"invocation.json") for i in 1:132])
    add_global(name)=push!(paths,joinpath(run,name,"invocation.json"))
    if phase=="source-math"
        return paths
    end
    add_points("native123")
    if phase!="native-select3"
        add_global("native-select3")
        first=selection(ctx,run,"native-select3")
        first["selected_level"]==0 && add_points("native4")
    end
    if phase in ("derivative-select3","diagnostic")
        add_global("native-selection");add_points("tailref123")
    end
    if phase=="diagnostic"
        add_global("derivative-select3");add_global("source-math")
        derivative=selection(ctx,run,"derivative-select3")
        derivative["selected_reference_level"]==0 && add_points("reference4")
    end
    unique(paths)
end

"Independent global reader: scalar native comparisons, separate from producer selection."
function native_reader(ctx,initial;fourth=nothing)
    ids=all_ids(ctx);policy=ctx.contract["numerical_policy"]
    length(initial)==132 && [p["point"]["id"] for p in initial]==ids || error("incomplete native reader coverage")
    maxima=Dict{String,Any}();comparisons=Dict{String,Any}[]
    for level in 2:(fourth===nothing ? 3 : 4)
        errors=Float64[]
        for i in 1:132
            coarse=initial[i]["levels"][level==4 ? 3 : level-1]["values"]
            fine=level==4 ? fourth[i]["levels"][1]["values"] : initial[i]["levels"][level]["values"]
            e=maximum(abs(fine[k]-coarse[k])/max(abs(coarse[k]),policy["relative_floor"]) for k in ("UU","FUT","amplitude"))
            push!(errors,e);push!(comparisons,Dict("point_id"=>ids[i],"level"=>level,"max_scaled"=>e))
        end
        maxima[string(level)]=maximum(errors)
    end
    level=maxima["3"]<=policy["physical_value_max_scaled"] ? 3 : fourth!==nothing && maxima["4"]<=policy["physical_value_max_scaled"] ? 4 : 0
    status=level>0 ? "COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" : fourth===nothing ? "NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" : "FAIL"
    (;level,status,maxima,comparisons)
end

"Independent full132 common-step decision after decimal operands were recomputed."
function derivative_reader(ctx,initial,native_level;fourth=nothing)
    length(initial)==132 && [p["point"]["id"] for p in initial]==all_ids(ctx) || error("incomplete derivative reader coverage")
    all(p["native_level"]==native_level for p in initial) || error("mixed native level in independent derivative reader")
    setprecision(256) do
        maxerror(packet,level)=begin
            row=only(r for r in packet["reference_steps"] if r["level"]==level)
            all(r["level"] isa Integer && !(r["level"] isa Bool) for r in packet["reference_steps"]) || error("invalid reference level")
            numbers=vcat(collect(values(packet["centre_errors_decimal"])),collect(values(row["errors_decimal"])))
            maximum(CollinsFit.sidis_job_parse_decimal(v) for v in numbers)
        end
        maxima=Dict{String,Any}("3"=>string(maximum(maxerror(p,3) for p in initial)))
        threshold=BigFloat(ctx.contract["derivative_policy"]["gradient_max_scaled"])
        level=parse(BigFloat,maxima["3"])<=threshold ? 3 : 0
        if fourth!==nothing
            length(fourth)==132 && [p["point"]["id"] for p in fourth]==all_ids(ctx) &&
                all(p["native_level"]==native_level for p in fourth) || error("incomplete or mixed fourth-step coverage")
            maxima["4"]=string(maximum(maxerror(p,4) for p in fourth))
            level==0 && parse(BigFloat,maxima["4"])<=threshold && (level=4)
        end
        status=level>0 ? "COMMON_DERIVATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE" : fourth===nothing ? "DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" : "FAIL"
        (;level,status,maxima)
    end
end

function ledger(ctx,run,phase)
    rows=Dict{String,Any}[]
    for path in global_dependencies(ctx,run,phase)
        endswith(path,"/invocation.json") || continue
        inv=JSON.parsefile(path)
        push!(rows,Dict("path"=>path,"sha256"=>digest(path),"phase"=>inv["phase"],
            "point_index"=>get(inv,"point_index",nothing),"source_fingerprint"=>inv["source_fingerprint"],
            "complete_command_resources"=>Dict(step["name"]=>step["resource"] for step in inv["steps"])))
    end
    rows
end

function expected_rows(ctx,initial,native_level,reference_level)
    [merge(deepcopy(point),Dict("UU"=>packet["values"]["UU"],"FUT"=>packet["values"]["FUT"],
        "amplitude"=>packet["values"]["amplitude"],
        "unweighted_observed_minus_predicted"=>point["measured"]["central_value"]-packet["values"]["amplitude"],
        "native_level"=>native_level,"reference_level"=>reference_level)) for (point,packet) in zip(ctx.request.data["points"],initial)]
end

function audit_global(ctx,run,phase,record)
    validate_header(record,ctx,phase;final=phase=="diagnostic")
    record["dependency_sha256"]==Dict(p=>digest(p) for p in global_dependencies(ctx,run,phase)) || error("incomplete global dependency ledger")
    if phase=="source-math"
        return HistoricalSIDISSourceMath.audit(record,ctx)
    end
    initial=point_packets(ctx,run,"native123")
    first=native_reader(ctx,initial)
    fourth=phase!="native-select3" && first.level==0 ? point_packets(ctx,run,"native4") : nothing
    native=native_reader(ctx,initial;fourth)
    if phase in ("native-select3","native-selection")
        expected=phase=="native-select3" ? first : native
        record["status"]==expected.status && record["selected_level"]==expected.level &&
            record["maximum_scaled_by_level"]==expected.maxima && record["comparisons"]==expected.comparisons &&
            record["source_ids"]==all_ids(ctx) || error("changed common native selection")
        return Dict("native_comparisons_recomputed"=>length(expected.comparisons),"native_level"=>expected.level)
    end
    native.level in (3,4) || error("unaccepted common native level")
    selected=global_packet(ctx,run,"native-selection")
    selected["selected_level"]==native.level && selected["maximum_scaled_by_level"]==native.maxima || error("stale common native report")
    refs=point_packets(ctx,run,"tailref123")
    first_derivative=derivative_reader(ctx,refs,native.level)
    if phase=="derivative-select3"
        record["status"]==first_derivative.status && record["selected_reference_level"]==first_derivative.level &&
            record["maximum_scaled_by_level_decimal"]==first_derivative.maxima && record["source_ids"]==all_ids(ctx) &&
            record["native_level"]==native.level && record["coordinates"]==13 || error("changed common derivative selection")
        return Dict("all132_all13_saved_reference_sets_reviewed"=>132,"selected_reference_level"=>first_derivative.level)
    end
    phase=="diagnostic" || error("unknown global reader phase")
    fourth_refs=first_derivative.level==0 ? point_packets(ctx,run,"reference4") : nothing
    derivative=derivative_reader(ctx,refs,native.level;fourth=fourth_refs)
    derivative.level in (3,4) || error("unaccepted full132/all13 derivative references")
    choice=global_packet(ctx,run,"derivative-select3")
    choice["selected_reference_level"]==first_derivative.level && choice["maximum_scaled_by_level_decimal"]==first_derivative.maxima || error("stale initial derivative decision")
    math=global_packet(ctx,run,"source-math");mathreview=HistoricalSIDISSourceMath.audit(math,ctx)
    record["status"]=="CONDITIONAL_FULL132_POINT_DIAGNOSTIC_RECORDED_NOT_FITTING" &&
        record["native_level"]==native.level && record["reference_level"]==derivative.level &&
        record["source_ids"]==all_ids(ctx) && record["coordinates"]==13 &&
        record["full132_native_tail_reference_derivative_policy_completed"]===true &&
        record["points"]==expected_rows(ctx,refs,native.level,derivative.level) &&
        record["job_ledger"]==ledger(ctx,run,phase) || error("incomplete, stale or changed full132 diagnostic ledger")
    Dict("native_point_records_reviewed"=>fourth===nothing ? 132 : 264,
        "all132_all13_saved_reference_sets_reviewed"=>fourth_refs===nothing ? 132 : 264,
        "selected_reference_level"=>derivative.level,"source_math"=>mathreview,
        "native_evolution_physical_references_AD_and_RSS_reexecuted"=>false)
end

end
