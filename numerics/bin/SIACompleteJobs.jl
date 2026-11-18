"Complete-job provenance and explicit frozen-library dispatch for SIA."
module SIACompleteJobs
using CollinsFit, JSON, SHA, Dates
include("SIAJobBindings.jl")
include("SIAGeometryGlobalAudit.jl")
include("SIASourceControls.jl")
const B=SIAJobBindings
const G=SIAGeometryGlobalAudit
const R=G.R
const N=G.N
const S=G.S
const FLAGS=(N.FLAGS...,)
const PHYSICAL=("native123","native4","tailref","derivative")
const GLOBALS=("native-select3","native-selection","tail-selection","diagnostic")
const PRODUCTS=("report.json","audit.json","producer.log","audit.log","producer.resources.txt","audit.resources.txt")
const COMPLETE="COMPLETE_SIA_JULIA_JOB_PAIR_NOT_FITTING_GATE"
const REVIEWED="COMPLETE_SIA_SAVED_JOB_ARITHMETIC_REVIEWED_NOT_GATE"
const COMPUTED="COMPLETE_SIA_JOB_VALUES_COMPUTED_NOT_GATE"
require(value,message)=value || error(message)
digest(path)=B.digest(path)
readjson(path)=JSON.parsefile(path)

function hashes(inventory)
    require(inventory isa AbstractDict && !isempty(inventory) && all(isfile(p) && digest(p)==h for (p,h) in inventory),
        "missing or changed complete file inventory")
end

function canonical(path)
    path=String(path)
    require(isabspath(path) && normpath(path)==path && realpath(path)==path,"canonical non-aliased run file required")
    path
end

function load_manifest(path,external_sha)
    canonical(path);require(digest(path)==external_sha,"changed externally bound execution manifest")
    m=readjson(path)
    require(m["schema_version"]===1 && m["status"]=="PREPARED_COMPLETE_SIA_JOB_RUN_NOT_PHYSICS_OR_FITTING_ADMISSION" &&
        m["scope"]==S.SCOPE && m["theory_scope"]=="reconstructed_theory" &&
        all(get(m,k,nothing)===false for k in (FLAGS...,"physical_dispatch_started")),"wrong manifest scope or attempted admission")
    require(m["manifest"]==path && basename(path)=="execution.json" && dirname(path)==m["run"] &&
        dirname(m["source_plan"])==m["run"],"wrong execution/source-plan namespace")
    plan=readjson(m["source_plan"])
    context=B.load(plan["frozen"],m["source_plan"],m["source_plan_sha256"])
    require(S.declaration(context.request.data)=="reconstructed_theory","synthetic operands cannot enter the physical dispatch loader")
    policy=readjson(joinpath(B.ADAPTER,"config/historical-sia-job-contract.json"))
    require(m["adapter_head"]==context.plan["adapter_head"] && S.same_json(m["resources"],policy["resources"]) &&
        S.same_json(m["limitations"],context.plan["limitations"]),"changed adapter, limits or scientific caveats")
    native=CollinsFit.verified_hoppet_build(m["native"]);parton=CollinsFit.verify_parton_build(m["parton"])
    require(any(basename(p)=="hoppet_bridge.f90" && h==digest(joinpath(context.frozen,"fortran/hoppet_bridge.f90"))
        for (p,h) in native["source_hashes"]) &&
        any(basename(p)=="lhapdf.cpp" && h==digest(joinpath(context.frozen,"bridges/lhapdf.cpp"))
        for (p,h) in parton["source_hashes"]),"native or parton bridge differs from frozen numerical source")
    expected=Dict{String,Any}()
    for (file,receipt) in ((m["native"],native),(m["parton"],parton))
        for key in ("source_hashes","output_hashes"),(name,h) in receipt[key]
            require(!haskey(expected,name) || expected[name]==h,"conflicting build identities");expected[name]=h
        end
        expected[file]=digest(file)
    end
    require(S.same_json(m["build_inputs"],expected) && length(m["retained_files"])==9,"incomplete build or original9 inventory")
    hashes(expected);hashes(m["retained_files"]);hashes(m["source_control_files"])
    require(digest(m["sidis_owner_file"])==m["sidis_owner_file_sha256"] &&
        S.same_json(readjson(m["sidis_owner_file"]),m["sidis_owner"]),"changed bound current SIDIS owner")
    binding=merge(deepcopy(context.binding),Dict("execution_manifest_sha256"=>external_sha,
        "native_receipt_sha256"=>digest(m["native"]),"parton_receipt_sha256"=>digest(m["parton"])))
    (source=context,m=m,path=path,sha=external_sha,binding=binding)
end

function recheck(c)
    B.recheck(c.source)
    require(digest(c.path)==c.sha,"execution manifest changed during job")
    hashes(c.m["build_inputs"]);hashes(c.m["retained_files"]);hashes(c.m["source_control_files"])
    require(digest(c.m["sidis_owner_file"])==c.m["sidis_owner_file_sha256"],"SIDIS owner changed during job")
end

function folder(c,phase,index=0)
    require(index isa Int && ((phase in GLOBALS && index==0) ||
        (phase in PHYSICAL && index in (phase=="derivative" ? (1,15,36) : 1:36))),"undeclared phase or geometry")
    root=phase in GLOBALS ? joinpath(c.m["run"],"global",phase) : joinpath(c.m["run"],phase,"g"*lpad(index,4,'0'))
    if ispath(root);canonical(root);end
    root
end

function needs4(c)
    readjson(joinpath(folder(c,"native-select3"),"report.json"))["selection"]["status"]=="NATIVE4_REQUIRED_NOT_DIAGNOSTIC_GATE"
end

function dependencies(c,phase)
    first=[("native123",i) for i in 1:36]
    phase=="native123" && return Tuple{String,Int}[]
    phase=="native-select3" && return first
    phase=="native4" && return [("native-select3",0)]
    native=vcat(first,[("native-select3",0)])
    needs4(c) && append!(native,[("native4",i) for i in 1:36])
    phase=="native-selection" && return native
    push!(native,("native-selection",0))
    phase=="tailref" && return native
    append!(native,[("tailref",i) for i in 1:36])
    phase=="tail-selection" && return native
    push!(native,("tail-selection",0))
    phase=="derivative" && return native
    phase=="diagnostic" && return vcat(native,[("derivative",i) for i in (1,15,36)])
    error("undeclared dependency phase")
end

function plan_inputs(c)
    plan=c.source.plan;result=Dict{String,Any}()
    for path in (c.source.planfile,plan["request"],plan["source_request"],"/usr/bin/time",c.m["julia_launcher"])
        result[path]=digest(path)
    end
    for (root,inventory) in ((plan["adapter"],plan["adapter_files"]),(plan["frozen"],plan["frozen_files"]))
        for (name,h) in inventory;result[joinpath(root,name)]=h;end
    end
    merge!(result,plan["source_inputs"])
    result
end

function base_inputs(c)
    result=plan_inputs(c)
    result[c.path]=c.sha;result[c.m["sidis_owner_file"]]=c.m["sidis_owner_file_sha256"]
    for inventory in (c.m["build_inputs"],c.m["retained_files"],c.m["source_control_files"])
        merge!(result,inventory)
    end
    require(S.same_json(c.m["execution_inputs"],filter(pair->first(pair)!=c.path,result)),"changed immutable execution input set")
    result
end

function job_inputs(c,phase,index,base)
    result=copy(base)
    for (p,i) in dependencies(c,phase),name in (PRODUCTS...,"initial-invocation.json","invocation.json","job.completed","lease.json")
        file=joinpath(folder(c,p,i),name);canonical(file);result[file]=digest(file)
    end
    lease=joinpath(folder(c,phase,index),"lease.json");result[lease]=digest(lease)
    result
end

function command(c,phase,index,output,name;report_sha=nothing)
    root=folder(c,phase,index)
    args=[c.path,c.sha,phase,string(index)]
    program=name=="producer" ? "run_sia_geometry_job.jl" : "audit_sia_geometry_job.jl"
    name=="audit" && append!(args,[joinpath(root,"report.json"),report_sha])
    push!(args,joinpath(output,name=="producer" ? "report.json" : "audit.json"))
    ["/usr/bin/time","--output",joinpath(output,name*".resources.txt"),"--format",
        "COLLINSFLOW_TIME_V1\t%M\t%e\t%U\t%S\t%x","--",c.m["julia_launcher"],"--startup-file=no",
        "--project="*joinpath(c.source.frozen,"numerics"),joinpath(B.ADAPTER,"numerics/bin",program),args...]
end

function resource(path)
    rows=filter(line->startswith(line,"COLLINSFLOW_TIME_V1\t"),readlines(path))
    require(length(rows)==1,"one complete GNU-time record required")
    fields=split(only(rows),'\t');require(length(fields)==6,"incomplete GNU-time record")
    rss=parse(Int,fields[2]);wall,user,system=parse.(Float64,fields[3:5]);code=parse(Int,fields[6])
    require(rss>0 && 0<=code<=255 && all(x->isfinite(x)&&x>=0,(wall,user,system)),"invalid GNU-time operands")
    Dict("maximum_RSS_KiB"=>rss,"wall_seconds"=>wall,"user_seconds"=>user,"system_seconds"=>system,
        "child_exit_code"=>code,"scope"=>"GNU_time_direct_child_maximum_not_simultaneous_process_tree_peak")
end

function steps(record,budget,expected_commands,root)
    require([step["name"] for step in record["steps"]]==["producer","audit"] && S.same_json(record["resources"],budget),
        "exact measured producer/auditor pair required")
    for (step,expected) in zip(record["steps"],expected_commands)
        name=step["name"];measured=resource(joinpath(root,name*".resources.txt"))
        require(step["command"]==expected && step["returncode"]===0 && step["process_returncode"]===0 &&
            step["interruption"]===nothing && step["termination_signal"]===nothing && measured["child_exit_code"]===0 &&
            S.same_json(step["resource"],measured) && measured["maximum_RSS_KiB"]<=budget["maximum_command_RSS_KiB"] &&
            S.finite_number(step["wrapper_elapsed_seconds"]) && step["wrapper_elapsed_seconds"]>=0 &&
            step["minimum_available_memory_observed_before_command"] isa Int &&
            step["minimum_available_memory_observed_before_command"]>=budget["minimum_available_memory_bytes"] &&
            step["log_sha256"]==digest(joinpath(root,name*".log")) &&
            step["resource_sha256"]==digest(joinpath(root,name*".resources.txt")),"changed command, raw/effective exit or resource proof")
    end
    require(sum(s["resource"]["wall_seconds"] for s in record["steps"])<=budget["whole_Julia_jobs_time_budget_seconds"] &&
        sum(s["wrapper_elapsed_seconds"] for s in record["steps"])<=budget["whole_Julia_jobs_time_budget_seconds"],"whole Julia pair budget exceeded")
    nothing
end

function verify_pair(c,phase,index,base)
    root=folder(c,phase,index);path=joinpath(root,"invocation.json");canonical(path)
    record=readjson(path);sha=digest(path)
    require(record["status"]==COMPLETE && record["scope"]=="source_bound_complete_sia_Julia_job_pair" &&
        record["theory_scope"]==S.declaration(c.source.request.data) && record["record_role"]=="complete_job_execution_invocation" &&
        record["manifest_sha256"]==c.sha && record["phase"]==phase && record["index"]===index &&
        record["completed_at"] isa String && !isempty(record["completed_at"]) && all(get(record,k,nothing)===false for k in FLAGS) &&
        strip(read(joinpath(root,"job.completed"),String))==sha,"incomplete or substituted whole job")
    expected=job_inputs(c,phase,index,base)
    require(S.same_json(record["inputs"],expected),"incomplete or stale job input/dependency ledger")
    hashes(expected)
    initial=readjson(joinpath(root,"initial-invocation.json"))
    budget=c.m["resources"][phase in PHYSICAL ? phase : "global"]
    require(record["initial_record_sha256"]==digest(joinpath(root,"initial-invocation.json")) && initial["status"]=="FAIL" &&
        initial["record_role"]=="initial_sentinel_not_an_executed_failure" && initial["completed_at"]===nothing &&
        isempty(initial["steps"]) && initial["manifest_sha256"]==c.sha && initial["scope"]=="source_bound_complete_sia_Julia_job_pair" &&
        initial["theory_scope"]==S.declaration(c.source.request.data) && initial["phase"]==phase && initial["index"]===index &&
        S.same_json(initial["resources"],budget) && initial["started_at"] isa String && !isempty(initial["started_at"]) &&
        all(get(initial,k,nothing)===false for k in FLAGS),"changed initial nonterminal sentinel")
    require(Set(keys(record["outputs"]))==Set(joinpath(root,name) for name in PRODUCTS),"incomplete whole-job output inventory")
    hashes(record["outputs"])
    reportsha=digest(joinpath(root,"report.json"))
    steps(record,c.m["resources"][phase in PHYSICAL ? phase : "global"],
        [command(c,phase,index,root,"producer"),command(c,phase,index,root,"audit";report_sha=reportsha)],root)
    report=readjson(joinpath(root,"report.json"));audit=readjson(joinpath(root,"audit.json"))
    validate_header(c,phase,index,report,COMPUTED);validate_header(c,phase,index,audit,REVIEWED)
    require(all(get(row,"completed_at",nothing) isa String && !isempty(row["completed_at"]) for row in (report,audit)) &&
        audit["report_sha256"]==reportsha && audit["complete_prerequisite_pair_ledgers_verified"]===length(dependencies(c,phase)) &&
        all(get(audit,key,nothing)===false for key in ("native_evolution_reexecuted_by_reader",
            "cold_AD_or_integrals_reexecuted_by_reader","raw_primary_tables_replayed_in_Julia")),
        "missing complete independent saved-record identity or enlarged reader claim")
    phase in GLOBALS && require(get(audit,"complete72_saved_arithmetic_verified",nothing)===(phase=="diagnostic"),
        "changed complete72 saved arithmetic scope")
    report
end

"Verify an actual Linux descriptor; dispatch supplies its fixed expected paths."
function resource_lease_identity(proof,expected_path)
    require(proof["pid"] isa Int && proof["pid"]>0 && proof["start_ticks"] isa Int && proof["start_ticks"]>0 &&
        proof["file_descriptor"] isa Int && proof["file_descriptor"]>=0,"typed positive resource process identity required")
    boot=strip(read("/proc/sys/kernel/random/boot_id",String))
    pid=proof["pid"]
    text=read("/proc/$pid/stat",String);fields=split(text[findlast(')',text)+1:end])
    require(startswith(text,string(pid)*" (") && fields[1]!="Z" && proof["boot_id"]==boot &&
        parse(Int,fields[20])==proof["start_ticks"] && proof["lock_path"]==expected_path &&
        isabspath(expected_path) && normpath(expected_path)==expected_path && realpath(expected_path)==expected_path,
        "expired, aliased or substituted live resource lease")
    fd="/proc/$pid/fd/$(proof["file_descriptor"])";info=stat(fd)
    require(realpath(fd)==expected_path && isfile(fd) && info.inode==proof["inode"] && info.device==proof["device"] &&
        info.uid==ccall(:getuid,Cuint,()) && info.nlink==1 && info.mode & 0o077==0 &&
        occursin(Regex("FLOCK\\s+ADVISORY\\s+WRITE\\s+"*string(pid)*"\\s"),read("/proc/$pid/fdinfo/$(proof["file_descriptor"])",String)),
        "resource lease is not a private live exclusive flock")
    (pid=pid,start_ticks=proof["start_ticks"],boot_id=boot,device=info.device,inode=info.inode)
end

function live_lease(c,phase,index)
    lease=readjson(joinpath(folder(c,phase,index),"lease.json"))
    boot=strip(read("/proc/sys/kernel/random/boot_id",String))
    for (key,path) in (("family","/tmp/collinsflow-sia-family-$(ccall(:getuid,Cuint,())).lock"),
        ("geometry","/tmp/collinsflow-sia-geometry-$(ccall(:getuid,Cuint,())).lock"))
        resource_lease_identity(lease[key],path)
    end
    require(lease["family"]["pid"]==lease["geometry"]["pid"] &&
        lease["host_observation"]["status"]=="HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION" &&
        lease["host_observation"]["live_host_observed"]===true,"owner or host observation cannot authorize this job")
    # The actual worker must descend from the process that holds both leases.
    ancestor=getpid();found=false
    for _ in 1:8
        ancestor==lease["family"]["pid"] && (found=true;break)
        ancestor>1 || break
        text=read("/proc/$ancestor/stat",String);fields=split(text[findlast(')',text)+1:end]);ancestor=parse(Int,fields[2])
    end
    require(found,"worker is not a child of the live lease holder")
    owner=c.m["sidis_owner"]
    if owner["boot_id"]==boot && isfile("/proc/$(owner["pid"])/stat")
        text=read("/proc/$(owner["pid"])/stat",String);fields=split(text[findlast(')',text)+1:end])
        require(parse(Int,fields[20])!=owner["start_ticks"] || fields[1]=="Z","bound SIDIS owner still owns host")
    end
end

module FrozenLibrary end
const LIBRARY=Ref{Union{Nothing,String}}(nothing)
function load_library(frozen)
    if LIBRARY[]===nothing
        Base.include(FrozenLibrary,joinpath(frozen,"numerics/bin/predict_historical_sia_integrals.jl"))
        Base.include(FrozenLibrary,joinpath(frozen,"numerics/bin/HistoricalSIAIntegralAudit.jl"))
        LIBRARY[]=frozen
    end
    require(LIBRARY[]==frozen,"do not mix frozen numerical libraries")
    FrozenLibrary
end

function source_controls(c)
    root=c.m["source_controls"];record=readjson(joinpath(root,"invocation.json"));sha=digest(joinpath(root,"invocation.json"))
    require(record["status"]=="COMPLETE_SYNTHETIC_SIA_SOURCE_CONTROL_PAIR_NOT_PHYSICAL_GATE" &&
        record["scope"]==SIASourceControls.SCOPE && record["record_role"]=="complete_source_control_pair_invocation" &&
        record["completed_at"] isa String && record["theory_scope"]=="synthetic_closure" &&
        record["source_plan"]==c.source.planfile && record["source_plan_sha256"]==c.source.plan_sha &&
        record["julia_launcher"]==c.m["julia_launcher"] && record["adapter_head"]==c.source.plan["adapter_head"] &&
        strip(read(joinpath(root,"job.completed"),String))==sha &&
        all(get(record,k,nothing)===false for k in FLAGS),"missing complete synthetic source-control pair")
    require(S.same_json(record["input_hashes"],plan_inputs(c)),"source-control complete input inventory differs from current plan")
    hashes(record["input_hashes"]);hashes(record["outputs"])
    require(Set(keys(record["outputs"]))==Set(joinpath(root,name) for name in PRODUCTS) &&
        record["initial_record_sha256"]==digest(joinpath(root,"initial-invocation.json")),"incomplete source-control retained files")
    initial=readjson(joinpath(root,"initial-invocation.json"))
    require(initial["status"]=="FAIL" && initial["record_role"]=="initial_sentinel_not_an_executed_failure" &&
        isempty(initial["steps"]) && !haskey(initial,"completed_at") && initial["source_plan_sha256"]==c.source.plan_sha,
        "changed source-control initial nonterminal sentinel")
    function source_command(name,output)
        args=[c.source.frozen,c.source.planfile,c.source.plan_sha]
        if name=="audit";append!(args,[joinpath(root,"report.json"),digest(joinpath(root,"report.json"))]);end
        push!(args,joinpath(output,name=="producer" ? "report.json" : "audit.json"))
        ["/usr/bin/time","--output",joinpath(output,name*".resources.txt"),"--format",
            "COLLINSFLOW_TIME_V1\t%M\t%e\t%U\t%S\t%x","--",c.m["julia_launcher"],"--startup-file=no",
            "--project="*joinpath(c.source.frozen,"numerics"),joinpath(B.ADAPTER,"numerics/bin",
                name=="producer" ? "run_sia_source_controls.jl" : "audit_sia_source_controls.jl"),args...]
    end
    steps(record,c.m["resources"]["global"],[source_command("producer",root),source_command("audit",root)],root)
    report=readjson(joinpath(root,"report.json"));audit=readjson(joinpath(root,"audit.json"))
    expected=SIASourceControls.audit(c.source.request.data,c.source.parent,c.source.binding,report)
    require(all(haskey(audit,k) && S.same_json(audit[k],v) for (k,v) in expected) &&
        audit["source_control_report_sha256"]==digest(joinpath(root,"report.json")),"changed source-control arithmetic review")
    fresh=c.m["source_control_fresh_review"]
    require(digest(fresh["path"])==fresh["sha256"] && haskey(c.m["source_control_files"],fresh["path"]),"unbound fresh source-control cache review")
    cached=readjson(fresh["path"]);cachedroot=dirname(fresh["path"])
    require(cached["status"]=="FRESH_CACHED_SYNTHETIC_SOURCE_CONTROL_REVIEW_COMPLETE_NOT_PHYSICAL_GATE" &&
        cached["scope"]==SIASourceControls.SCOPE && cached["theory_scope"]=="synthetic_closure" &&
        cached["completed_at"] isa String && all(get(cached,k,nothing)===false for k in FLAGS) &&
        cached["command"]==source_command("audit",cachedroot) && cached["audit_sha256"]==digest(joinpath(cachedroot,"audit.json")) &&
        cached["source_control_report_sha256"]==digest(joinpath(root,"report.json")),"stale or absent source-control cached review")
    hashes(cached["input_hashes"])
    measured=resource(joinpath(cachedroot,"audit.resources.txt"));budget=c.m["resources"]["global"]
    require(S.same_json(cached["resource"],measured) && cached["returncode"]===0 && cached["process_returncode"]===0 &&
        cached["interruption"]===nothing && cached["termination_signal"]===nothing && measured["child_exit_code"]===0 &&
        measured["maximum_RSS_KiB"]<=budget["maximum_command_RSS_KiB"] && measured["wall_seconds"]<=budget["whole_Julia_jobs_time_budget_seconds"] &&
        S.finite_number(cached["wrapper_elapsed_seconds"]) && 0<=cached["wrapper_elapsed_seconds"]<=budget["whole_Julia_jobs_time_budget_seconds"] &&
        cached["minimum_available_memory_observed_before_command"]>=budget["minimum_available_memory_bytes"] &&
        cached["log_sha256"]==digest(joinpath(cachedroot,"audit.log")) &&
        cached["resource_sha256"]==digest(joinpath(cachedroot,"audit.resources.txt")),"failed or changed cached-control resource proof")
    cached_audit=readjson(joinpath(cachedroot,"audit.json"))
    require(all(haskey(cached_audit,k) && S.same_json(cached_audit[k],v) for (k,v) in expected),"cached source-control arithmetic changed")
    expected
end

function common_header(c,phase,index,status)
    theory=S.declaration(c.source.request.data)
    result=merge(N.flags(),Dict{String,Any}("schema_version"=>1,"status"=>status,"scope"=>S.SCOPE,
        "theory_scope"=>theory,"source_context"=>deepcopy(c.binding),"manifest_sha256"=>c.sha,
        "phase"=>phase,"index"=>index))
    theory=="synthetic_closure" && (result["fixture_scope"]="explicitly_synthetic_complete_job_arithmetic_not_physical_dispatch")
    result
end

function validate_header(c,phase,index,report,status)
    expected=common_header(c,phase,index,status)
    require(all(haskey(report,k) && S.same_json(report[k],v) for (k,v) in expected),"changed complete report identity or premature admission")
end

function collect_dependencies(c,phase,legacy)
    base=base_inputs(c);hashes(base);result=Dict{Tuple{String,Int},Any}()
    inventory=copy(base)
    for (p,i) in dependencies(c,phase),name in (PRODUCTS...,"initial-invocation.json","invocation.json","job.completed","lease.json")
        path=joinpath(folder(c,p,i),name);canonical(path);inventory[path]=digest(path)
    end
    for (p,i) in dependencies(c,phase)
        report=verify_pair(c,p,i,base)
        validate_header(c,p,i,report,COMPUTED)
        result[(p,i)]=report
    end
    # Independently check every leaf in this complete prerequisite ledger.
    for ((p,i),report) in result
        if p in PHYSICAL
            leaf=report["leaf_report"];spec=leaf["definition"]
            checked=if p in ("native123","native4")
                N.audit(c.source.request.data,c.source.parent,spec,c.binding,leaf)
            elseif p=="tailref"
                R.audit_tail(c.source.request.data,c.source.parent,spec,c.binding,leaf)
            else
                R.audit_derivative(c.source.request.data,c.source.parent,spec,c.binding,leaf,legacy.HistoricalSIAIntegralAudit.inward_checks)
            end
            require(spec["phase"]==p && spec["geometry_index"]===i,"leaf identity differs from complete pair")
            stored=readjson(joinpath(folder(c,p,i),"audit.json"))
            require(S.same_json(stored["leaf_arithmetic_review"],checked),"saved leaf audit differs from the fresh arithmetic reader")
        end
    end
    hashes(inventory)
    (records=result,input_hashes=inventory)
end

function native_reports(c,records;first_only=false)
    result=Any[[records[("native123",i)]["leaf_report"] for i in 1:36]]
    !first_only && needs4(c) && push!(result,[records[("native4",i)]["leaf_report"] for i in 1:36])
    result
end

function selected_native(c,records)
    first=native_reports(c,records;first_only=true)
    coarse=records[("native-select3",0)]["selection"]
    G.audit_native(c.source.request.data,c.source.parent,first,coarse)
    all=native_reports(c,records)
    selection=records[("native-selection",0)]["selection"]
    chosen=G.audit_native(c.source.request.data,c.source.parent,all,selection)
    require(chosen.selected!==nothing,"native refinement exhausted; do not start later phases")
    chosen
end

function selected_tails(c,records)
    selected_native(c,records)
    tails=[records[("tailref",i)]["leaf_report"] for i in 1:36]
    selected=G.audit_tails(c.source.request.data,c.source.parent,native_reports(c,records),tails,
        records[("tail-selection",0)]["selection"])
    errors=G.reference_errors(c.source.request.data,c.source.parent,tails,selected.hi)
    require(selected.error<=1e-4 && all(<=(1e-4),errors),"all72 finite-tail or scalar adaptive reference target failed")
    for index in (1,15,36)
        if haskey(records,("derivative",index))
            leaf=records[("derivative",index)]["leaf_report"]
            expected=R.derivative_definition(c.source.request.data,c.source.parent,index,selected.native.selected,selected.hi.order)
            require(S.same_json(leaf["definition"],expected),"derivative job differs from global common native/order selection")
            groups,_=S.geometry_groups(c.source.request.data)
            require(S.same_json(leaf["value_points"],selected.hi.points[groups[index]]),"derivative value operands differ from selected B32 outputs")
        end
    end
    selected,errors
end

function definition(c,phase,index,records)
    data=c.source.request.data;parent=c.source.parent
    if phase in ("native123","native4")
        if phase=="native4"
            # This dependency includes all36 native123 operands through the
            # completed global pair's input inventory. Recompute its selection.
            first=[readjson(joinpath(folder(c,"native123",i),"report.json"))["leaf_report"] for i in 1:36]
            G.audit_native(data,parent,[first],records[("native-select3",0)]["selection"])
            require(needs4(c),"native4 was not requested by the complete native3 reader")
        end
        return N.definition(data,parent,phase,index)
    end
    native=selected_native(c,records)
    phase=="tailref" && return R.tail_definition(data,parent,index,native.selected)
    tail,_=selected_tails(c,records)
    R.derivative_definition(data,parent,index,native.selected,tail.hi.order)
end

function produce_global(c,phase,records)
    data=c.source.request.data;parent=c.source.parent
    if phase=="native-select3"
        return Dict("selection"=>G.produce_native(data,parent,native_reports(c,records;first_only=true)))
    elseif phase=="native-selection"
        first=native_reports(c,records;first_only=true)
        G.audit_native(data,parent,first,records[("native-select3",0)]["selection"])
        selection=G.produce_native(data,parent,native_reports(c,records))
        require(selection["selection_complete"],"all declared native refinements exhausted")
        return Dict("selection"=>selection)
    elseif phase=="tail-selection"
        selected_native(c,records)
        tails=[records[("tailref",i)]["leaf_report"] for i in 1:36]
        selection=G.produce_tails(data,parent,native_reports(c,records),tails)
        # Producer records the selected-order reference operands; the fresh
        # reader independently repeats selection and comparison arithmetic.
        checked=G.audit_tails(data,parent,native_reports(c,records),tails,selection)
        errors=G.reference_errors(data,parent,tails,checked.hi)
        require(checked.error<=1e-4 && all(<=(1e-4),errors),"finite-tail or independent scalar target failed")
        return Dict("selection"=>selection,"selected_order_reference_errors"=>errors)
    end
    tail,errors=selected_tails(c,records)
    Dict("points"=>tail.hi.points,"selected_native_level"=>tail.native.selected,"selected_momentum_order"=>tail.hi.order,
        "selected_order_reference_errors"=>errors,"finite24to32_max_scaled"=>tail.error,
        "derivative_geometry_indices"=>[1,15,36],"representative_saved_derivative_checks"=>27,
        "limitations"=>deepcopy(c.m["limitations"]),"infinite_tail_proved"=>false,
        "independent_evolution_verified"=>false,"all72_production_derivative_certificate"=>false)
end

function audit_global(c,phase,records,report)
    data=c.source.request.data;parent=c.source.parent
    if phase=="native-select3"
        G.audit_native(data,parent,native_reports(c,records;first_only=true),report["selection"])
    elseif phase=="native-selection"
        G.audit_native(data,parent,native_reports(c,records;first_only=true),records[("native-select3",0)]["selection"])
        selection=G.audit_native(data,parent,native_reports(c,records),report["selection"])
        require(selection.selected!==nothing,"common native refinement exhausted")
    elseif phase=="tail-selection"
        selected_native(c,records)
        tails=[records[("tailref",i)]["leaf_report"] for i in 1:36]
        selection=G.audit_tails(data,parent,native_reports(c,records),tails,report["selection"])
        errors=G.reference_errors(data,parent,tails,selection.hi)
        require(selection.error<=1e-4 && all(<=(1e-4),errors) &&
            S.same_json(report["selected_order_reference_errors"],errors),"all72 reference/tail target failed or changed")
    else
        tail,errors=selected_tails(c,records)
        require(S.same_json(report["points"],tail.hi.points) && report["selected_native_level"]===tail.native.selected &&
            report["selected_momentum_order"]===tail.hi.order && S.same_json(report["selected_order_reference_errors"],errors) &&
            report["finite24to32_max_scaled"]===tail.error && report["derivative_geometry_indices"]==[1,15,36] &&
            report["representative_saved_derivative_checks"]===27 && S.same_json(report["limitations"],c.m["limitations"]) &&
            all(get(report,k,nothing)===false for k in ("infinite_tail_proved","independent_evolution_verified","all72_production_derivative_certificate")),
            "changed complete72 diagnostic or enlarged scientific scope")
    end
end

function compute_leaf(view,parent,spec,binding,legacy,prepare;progress=stage->nothing,checkpoint=report->nothing)
    phase=spec["phase"]
    require(phase in PHYSICAL,"undeclared leaf computation")
    phase in ("native123","native4") ?
        N.compute(view,parent,spec,binding,legacy,(v,p,l)->prepare(v,p,l,24);progress=progress,checkpoint=checkpoint) :
        phase=="tailref" ? R.compute_tail(view,parent,spec,binding,legacy,prepare;progress=progress,checkpoint=checkpoint) :
        R.compute_derivative(view,parent,spec,binding,legacy,prepare;progress=progress,checkpoint=checkpoint)
end

function audit_leaf(data,parent,spec,binding,leaf,legacy)
    phase=spec["phase"]
    require(phase in PHYSICAL,"undeclared leaf arithmetic audit")
    phase in ("native123","native4") ? N.audit(data,parent,spec,binding,leaf) :
        phase=="tailref" ? R.audit_tail(data,parent,spec,binding,leaf) :
        R.audit_derivative(data,parent,spec,binding,leaf,legacy.HistoricalSIAIntegralAudit.inward_checks)
end

function execute_loaded(c,phase,index,output,legacy;input=nothing,input_sha=nothing)
    proof=collect_dependencies(c,phase,legacy);records=proof.records
    source_controls(c)
    if input===nothing
        live_lease(c,phase,index)
        report=common_header(c,phase,index,COMPUTED)
        if phase in PHYSICAL
            spec=definition(c,phase,index,records);view=N.M.SelectedSIAGeometry(c.source.request,index)
            inputs=N.native_inputs(c.m["native"],c.m["parton"])
            checkpoint=leaf->B.atomic_report(output,merge(common_header(c,phase,index,"IN_PROGRESS_NOT_A_GATE"),Dict("leaf_report"=>leaf)))
            progress=stage->B.atomic_report(joinpath(dirname(output),"progress.json"),stage)
            prepare=(view,parent,level,B)->N.prepare_geometry(view,parent,inputs,level,B)
            leaf=compute_leaf(view,c.source.parent,spec,c.binding,legacy,prepare;progress=progress,checkpoint=checkpoint)
            N.inputs_stable(inputs);report["leaf_report"]=leaf
        else
            merge!(report,produce_global(c,phase,records))
        end
        live_lease(c,phase,index)
    else
        require(digest(input)==input_sha && canonical(input)==joinpath(folder(c,phase,index),"report.json"),"changed or substituted complete producer report")
        if phase=="diagnostic" && dirname(output)!=folder(c,phase,index)
            verify_pair(c,phase,index,base_inputs(c))
        end
        saved=readjson(input);validate_header(c,phase,index,saved,COMPUTED)
        report=common_header(c,phase,index,REVIEWED)
        if phase in PHYSICAL
            spec=definition(c,phase,index,records);leaf=saved["leaf_report"]
            audit=audit_leaf(c.source.request.data,c.source.parent,spec,c.binding,leaf,legacy)
            report["leaf_arithmetic_review"]=audit
        else
            audit_global(c,phase,records,saved)
            report["complete72_saved_arithmetic_verified"]=(phase=="diagnostic")
        end
        require(digest(input)==input_sha,"producer report changed during independent audit")
        report["report_sha256"]=input_sha
        report["native_evolution_reexecuted_by_reader"]=false
        report["cold_AD_or_integrals_reexecuted_by_reader"]=false
        report["raw_primary_tables_replayed_in_Julia"]=false
        report["complete_prerequisite_pair_ledgers_verified"]=length(records)
    end
    recheck(c)
    hashes(proof.input_hashes)
    report
end

function main(args=ARGS;audit=false)
    length(args)==(audit ? 7 : 5) || error("usage: WORKER MANIFEST MANIFEST_SHA PHASE INDEX [REPORT REPORT_SHA] FRESH_OUTPUT")
    manifest,sha,phase,index=args[1:4];index=parse(Int,index);output=last(args)
    c=load_manifest(manifest,sha);folder(c,phase,index)
    output=B.fresh_output(output,c.source.planfile)
    # Failures are persisted even when the library load or the complete saved
    # ledger fails. An initial sentinel is not an executed failure.
    result=common_header(c,phase,index,"FAIL")
    B.atomic_report(output,result)
    try
        legacy=load_library(c.source.frozen)
        result=audit ? Base.invokelatest(execute_loaded,c,phase,index,output,legacy;input=args[5],input_sha=args[6]) :
            Base.invokelatest(execute_loaded,c,phase,index,output,legacy)
    catch exception
        result=readjson(output)
        result["status"]="FAIL";result["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        result["completed_at"]=string(now(UTC));B.atomic_report(output,result)
    end
    println(JSON.json(Dict("status"=>result["status"],"phase"=>phase,"index"=>index,"fit_performed"=>false)))
    result["status"]==(audit ? REVIEWED : COMPUTED) ? 0 : 2
end
end
