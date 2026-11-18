module SIDISRecoveryBindings
using JSON, SHA, Dates
include("SIDISRecoveryPaths.jl")
include("SIDISRecoveryLocations.jl")
include("SIDISContinuationPaths.jl")
const RP=SIDISRecoveryPaths
const RL=SIDISRecoveryLocations
const CP=SIDISContinuationPaths
const ADAPTER=dirname(dirname(@__DIR__))
const SCOPE="same_attempt_reconstructed_full132_recovery_not_fitting"
const PREDECESSOR="b524304c005091a5717de33141fc842fa07fbee1"
const DECLARATION="a88e9ba2d0059796a3dadc6c7b86707f36ded616"
const FIRST_MANIFEST="8743c5052a7b950812f21db35c83ae2008c43f2d32dbd6cfd12c602808fea83a"
const FIRST_LOSS="2afc9d062c0124897b45ed7507d09603fc9cfb0c6b7c781bb9064caa8bf7d1b7"
const ANCHOR_FIELDS=("original","original_inventory","review","integrity","loss","request","request_sha256",
    "native","parton","dependency_receipts","numerical_source_fingerprint","frozen","frozen_head","frozen_files")
const COMPLETED="COMPLETE_JOB_RECORDED_NOT_PRODUCTION_GATE"
const PREPARED="PREPARED_IMMUTABLE_RECOVERY_NOT_A_GATE"
digest(path)=RP.digest(path)
false_flags()=Dict{String,Any}(key=>false for key in RP.FLAGS)
require(ok,message)=ok || error(message)

struct LegacyRun
    root::String
    original::String
    transport::Dict{String,Any}
end

"Owned current locations plus immutable earlier producing contexts."
struct RecoveryRun
    paths::RL.ScopedLocations
    origins::Tuple
    transport::Dict{String,Any}
    function RecoveryRun(paths::RP.RecoveryRunPaths,origins::Tuple,transport::Dict{String,Any})
        new(RL.snapshot(paths,RP),origins,deepcopy(transport))
    end
end
RP.retained(run::RL.ScopedLocations,phase,index)=RL.retained(run,phase,index,RP)
RP.folder(run::RL.ScopedLocations,phase,index=0)=RL.folder(run,phase,index,RP)
RP.fresh_folder(run::RL.ScopedLocations,phase,index=0)=RL.fresh_folder(run,phase,index,RP)
RP.metadata(run::RL.ScopedLocations)=merge(false_flags(),Dict{String,Any}("scope"=>"immutable_whole_job_location_snapshot_not_admission",
    "saved_numerical_arithmetic_verified"=>false,"complete_job_provenance_verified"=>false,"numerical_dispatch_ready"=>false))
root(run::RecoveryRun)=run.paths.root
root(run::LegacyRun)=run.root
Base.joinpath(run::RecoveryRun,parts::AbstractString...)=RL.current_path(run.paths,parts...)
Base.joinpath(run::LegacyRun,parts::AbstractString...)=joinpath(run.root,parts...)

function point_folder(run::LegacyRun,phase,index)
    RP.point_index(phase,index)
    joinpath(phase=="native123" && index<=13 ? run.original : run.root,phase,"p"*lpad(string(index),4,'0'))
end
point_folder(run::RecoveryRun,phase,index)=RP.folder(run.paths,phase,index)

function origin(run::RecoveryRun,phase,index)
    row=RP.retained(run.paths,phase,index)
    row===nothing && return run
    row["protocol"]=="frozen_original" && return row["producer_root"]
    matches=filter(old->root(old)==row["producer_root"],run.origins)
    require(length(matches)==1,"missing or ambiguous actual producing context")
    only(matches)
end
origin(run::LegacyRun,phase,index)=phase=="native123" && index<=13 ? run.original : run

function transport(run,phase,index)
    owner=origin(run,phase,index)
    owner isa AbstractString ? nothing : owner.transport
end

function sha_file(path,sha)
    RP.canonical(path;directory=false)
    require(sha isa String && occursin(r"^[0-9a-f]{64}$",sha) && digest(path)==sha,"changed externally bound immutable file")
    path
end

function inventory(root)
    RP.canonical(root;directory=true)
    rows=Dict{String,Any}[]
    for (current,directories,files) in walkdir(root;follow_symlinks=false)
        for name in sort(vcat(directories,files))
            path=joinpath(current,name)
            if islink(path)
                push!(rows,Dict("path"=>path,"kind"=>"symlink_retained_as_link_not_followed","target"=>readlink(path)))
            elseif isfile(path)
                require(realpath(path)==path,"ordinary ancestor files cannot traverse aliases")
                push!(rows,Dict("path"=>path,"kind"=>"ordinary_file","bytes"=>filesize(path),"sha256"=>digest(path)))
            else
                require(isdir(path),"nonordinary ancestor output")
            end
        end
    end
    sort(rows;by=row->row["path"])
end

function verify_inventory(root,rows,workspace)
    expected=Dict{String,Any}[]
    for row in rows
        path=isabspath(row["path"]) ? row["path"] : joinpath(workspace,row["path"])
        require(normpath(path)==path && RP.inside(path,root) && path!=root,"ancestor inventory path outside its root")
        kind=get(row,"kind","ordinary_file");item=Dict{String,Any}("path"=>path,"kind"=>kind)
        if kind=="ordinary_file"
            require(row["bytes"] isa Int && row["bytes"]>=0,"exact nonnegative inventory size required")
            item["bytes"]=row["bytes"];item["sha256"]=row["sha256"]
        else
            require(kind=="symlink_retained_as_link_not_followed","unknown immutable inventory entry")
            item["target"]=row["target"]
        end
        push!(expected,item)
    end
    sort!(expected;by=row->row["path"])
    require(length(unique(row["path"] for row in expected))==length(expected) && inventory(root)==expected,"ancestor namespace inventory changed")
    expected
end

entry(phase,index,folder,producer,manifest,sha,protocol)=Dict{String,Any}("phase"=>phase,"index"=>index,
    "folder"=>folder,"producer_root"=>producer,"producer_manifest"=>manifest,"producer_manifest_sha256"=>sha,"protocol"=>protocol)

function original_locations(first,policy)
    producer=joinpath(first["original"],"full");manifest=joinpath(first["original"],"invocation.json")
    [entry("native123",index,joinpath(producer,"native123","p"*lpad(string(index),4,'0')),producer,
        manifest,policy["logical_attempt"]["original_outer_invocation_sha256"],"frozen_original") for index in 1:13]
end

function complete_locations(producer,protocol,manifest,sha)
    rows=Dict{String,Any}[]
    for phase in RP.POINT_PHASES
        phase_root=joinpath(producer,phase)
        require(!islink(phase_root) && (!ispath(phase_root) || realpath(phase_root)==phase_root),"aliased ancestor phase")
        isdir(phase_root) || continue
        for name in sort(readdir(phase_root))
            require(occursin(r"^p\d{4}$",name),"undeclared ancestor point folder")
            index=parse(Int,name[2:end]);RP.point_index(phase,index)
            folder=joinpath(phase_root,name);RP.canonical(folder;directory=true)
            path=joinpath(folder,"invocation.json");isfile(path) || continue
            inv=JSON.parsefile(path)
            if get(inv,"status",nothing)!=COMPLETED
                require(get(inv,"completed_at",nothing)===nothing && !ispath(joinpath(folder,"job.completed")),"terminal failed point cannot become a recovery candidate")
                continue
            end
            require(get(inv,"completed_at",nothing)!==nothing && inv["phase"]==phase && inv["point_index"] isa Int &&
                inv["point_index"]==index && isfile(joinpath(folder,"job.completed")),"false or incomplete point completion")
            require([step["name"] for step in inv["steps"]]==["producer","audit"],"complete producer/auditor pair required")
            for step in inv["steps"]
                require(all(step[key] isa Int && step[key]==0 for key in ("returncode","process_returncode")) &&
                    step["interruption"]===nothing && step["termination_signal"]===nothing,"failed/interrupted command cannot approve recovery")
            end
            push!(rows,entry(phase,index,folder,producer,manifest,sha,protocol))
        end
    end
    for phase in RP.GLOBAL_PHASES
        path=joinpath(producer,phase,"invocation.json");isfile(path) || continue
        inv=JSON.parsefile(path)
        require(get(inv,"status",nothing)==COMPLETED || get(inv,"completed_at",nothing)===nothing,"terminal failed global gate cannot be silently retried")
    end
    rows
end

function ordered(rows)
    result=sort(rows;by=row->(findfirst(==(row["phase"]),RP.POINT_PHASES),row["index"]))
    require(length(unique((row["phase"],row["index"]) for row in result))==length(result),"conflicting retained phase/index origins")
    for phase in RP.POINT_PHASES
        indices=[row["index"] for row in result if row["phase"]==phase]
        require(indices==collect(1:length(indices)),"retained phase must be a complete source-ordered prefix")
    end
    result
end

provenance(data,sha)=Dict{String,Any}("scope"=>SCOPE,"manifest_sha256"=>sha,
    "adapter_commit"=>data["adapter_commit"],"contract_sha256"=>RP.POLICY_SHA)
provenance(m)=provenance(m.data,m.sha)

function ancestors(nodes,frozen,predecessor,current_root,policy)
    require(nodes isa AbstractVector && 1<=length(nodes)<=64,"one to64 ordered ancestor nodes required")
    workspace=dirname(dirname(dirname(frozen)));seen=String[];contexts=Any[];rows=Dict{String,Any}[]
    first=nothing
    for (position,node) in enumerate(nodes)
        require(Set(keys(node))==Set(("protocol","manifest","manifest_sha256","inventory","inventory_sha256")),"complete externally bound ancestry schema required")
        path=sha_file(node["manifest"],node["manifest_sha256"]);evidence=sha_file(node["inventory"],node["inventory_sha256"])
        producer=dirname(path);RP.canonical(producer;directory=true)
        require(RP.inside(producer,joinpath(frozen,"runs")) && !RP.within_directory(producer,current_root) &&
            !RP.within_directory(current_root,producer) && !any(RP.within_directory(producer,parent) || RP.within_directory(parent,producer) for parent in seen),
            "repeated, nested, aliased or nonfamily ancestor root")
        parent=JSON.parsefile(path);loss=JSON.parsefile(evidence)
        if position==1
            require(node["protocol"]=="continuation_v1" && node["manifest_sha256"]==FIRST_MANIFEST && node["inventory_sha256"]==FIRST_LOSS &&
                basename(path)=="continuation.json" && parent["adapter"]==predecessor && parent["adapter_commit"]==PREDECESSOR &&
                parent["contract_sha256"]==CP.CONTRACT_SHA && parent["scope"]==CP.SCOPE && parent["status"]=="PREPARED_IMMUTABLE_CONTINUATION_NOT_A_GATE",
                "the exact supported stopped continuation is the first ancestor")
            first=parent;rows=original_locations(first,policy)
            original=joinpath(frozen,"runs/historical-sidis-jobs-real-data/attempt-2")
            require(first["original"]==original && !RP.within_directory(producer,original) &&
                !RP.within_directory(current_root,original) && !RP.within_directory(original,current_root),"wrong/overlapping original attempt2 anchor")
            require(loss["old_owner_absent_verified"]===true && loss["source_run"]==relpath(producer,workspace) &&
                loss["recorded_complete_prefix"]==37 && length(loss["retained_files"])==404,"wrong supported interruption inventory")
            context=LegacyRun(producer,joinpath(original,"full"),Dict{String,Any}("scope"=>CP.SCOPE,"manifest_sha256"=>node["manifest_sha256"],
                "adapter_commit"=>PREDECESSOR,"contract_sha256"=>CP.CONTRACT_SHA))
        else
            require(node["protocol"]=="recovery_v2" && basename(path)=="recovery.json" && parent["status"]==PREPARED &&
                parent["scope"]==SCOPE && parent["contract_sha256"]==RP.POLICY_SHA && parent["ancestors"]==nodes[1:position-1] &&
                parent["retained_jobs"]==ordered(rows) && loss["status"]=="STOPPED_SIDIS_NAMESPACE_INVENTORIED_NOT_A_GATE" &&
                loss["manifest_sha256"]==node["manifest_sha256"] && loss["run"]==producer && loss["old_owner_absent_verified"]===true,
                "inconsistent, cyclic or unreviewed recovery ancestry")
            files=CP.verify_git(parent["adapter"],parent["adapter_commit"])
            require(files==parent["adapter_files"] && all(get(files,k,nothing)==v for (k,v) in first["adapter_files"]),"changed ancestor adapter source")
            context=RecoveryRun(RP.RecoveryRunPaths(producer,parent["retained_jobs"]),Tuple(contexts),provenance(parent,node["manifest_sha256"]))
        end
        require(parent["run"]==producer && all(get(parent,key,nothing)==first[key] for key in ANCHOR_FIELDS) &&
            all(get(parent,key,nothing)===false for key in RP.FLAGS),"different/self-admitting numerical family")
        verify_inventory(producer,loss["retained_files"],workspace)
        append!(rows,complete_locations(producer,node["protocol"],path,node["manifest_sha256"]));rows=ordered(rows)
        push!(seen,producer);push!(contexts,context)
    end
    require(count(row->row["protocol"]=="continuation_v1",rows)==24,"predecessor must retain exactly its24 complete point pairs")
    (;first,rows,contexts=Tuple(contexts))
end

function load_manifest(frozen_arg,predecessor_arg,path_arg,sha_arg)
    frozen=RP.canonical(frozen_arg;directory=true);predecessor=RP.canonical(predecessor_arg;directory=true)
    path=sha_file(path_arg,sha_arg);m=JSON.parsefile(path);policy=RP.policy()
    require(basename(path)=="recovery.json" && m["status"]==PREPARED && m["scope"]==SCOPE && m["schema_version"]===2 &&
        m["contract_sha256"]==RP.POLICY_SHA && m["run"]==dirname(path) && m["frozen"]==frozen && m["predecessor"]==predecessor &&
        m["adapter"]==ADAPTER && all(get(m,key,nothing)===false for key in RP.FLAGS),"wrong or self-admitting recovery manifest")
    ffiles=CP.verify_git(frozen,CP.FROZEN_HEAD);pfiles=CP.verify_git(predecessor,PREDECESSOR);afiles=CP.verify_git(ADAPTER,m["adapter_commit"])
    require(length(ffiles)==335 && length(pfiles)==348 && m["frozen_head"]==CP.FROZEN_HEAD && m["predecessor_head"]==PREDECESSOR &&
        m["frozen_files"]==ffiles && m["predecessor_files"]==pfiles && m["adapter_files"]==afiles &&
        all(get(pfiles,k,nothing)==v for (k,v) in ffiles) && all(get(afiles,k,nothing)==v for (k,v) in pfiles),"changed335/348/additive source bindings")
    require(success(`git -C $ADAPTER merge-base --is-ancestor $DECLARATION $(m["adapter_commit"])`),"implementation must descend from its declaration")
    require(RP.inside(m["run"],joinpath(frozen,"runs")),"current namespace must be inside the explicit frozen checkout's runs")
    state=ancestors(m["ancestors"],frozen,predecessor,m["run"],policy)
    require(m["retained_jobs"]==state.rows && all(m[key]==state.first[key] for key in ANCHOR_FIELDS),"changed original family or retained-job map")
    anchor=policy["logical_attempt"]
    for (key,field) in (("review","independent_original13_review_sha256"),("integrity","post_review_original141_integrity_sha256"),("loss","original_loss_inventory_sha256"))
        sha_file(m[key],anchor[field])
    end
    workspace=dirname(dirname(dirname(frozen)));loss=JSON.parsefile(m["loss"])
    rows=[merge(row,Dict("path"=>joinpath(workspace,row["path"]))) for row in loss["interrupted_attempts"]["sidis132"]["retained_files"]]
    require(length(rows)==141 && rows==m["original_inventory"],"changed original141-file anchor")
    CP.verify_inventory(m["original"],rows)
    sha_file(joinpath(m["original"],"invocation.json"),anchor["original_outer_invocation_sha256"])
    deps=JSON.parsefile(joinpath(m["original"],"invocation.json"))["dependency_receipts"]
    require(deps==m["dependency_receipts"] && Set(keys(deps))==Set((m["native"],m["parton"])) && all(digest(p)==h for (p,h) in deps),"changed native/provider receipts")
    sha_file(m["request"],m["request_sha256"]);request=JSON.parsefile(m["request"])
    require(m["request"]==joinpath(m["original"],"full/request.json") && m["numerical_source_fingerprint"]==policy["frozen_numerics"]["source_fingerprint"] &&
        request["source_fingerprint"]==m["numerical_source_fingerprint"] && length(request["points"])==132 &&
        length(unique(row["id"] for row in request["points"]))==132,"wrong full132 source request")
    for (p,h) in request["source_inputs"];sha_file(p,h);end
    replay=m["source_replay"];replay_root=joinpath(m["run"],"source-replay")
    require(replay["request"]==joinpath(replay_root,"request.json") && replay["receipt"]==joinpath(replay_root,"replay.json"),"source replay belongs to the current namespace")
    sha_file(replay["request"],replay["request_sha256"]);sha_file(replay["receipt"],replay["receipt_sha256"])
    receipt=JSON.parsefile(replay["receipt"])
    require(JSON.parsefile(replay["request"])==request && receipt["status"]=="CURRENT_RAW_SOURCE_REPLAY_VERIFIED_NEW_OUTPUT_ONLY_NOT_NUMERICAL_GATE" &&
        receipt["original_request_sha256"]==m["request_sha256"] && receipt["fresh_request_sha256"]==replay["request_sha256"] &&
        all(get(receipt,key,nothing)===false for key in RP.FLAGS),"changed or failed raw source replay")
    runtime=m["runtime"];actual=realpath("/proc/self/exe")
    require(runtime["version"]==string(VERSION)=="1.12.6" && actual==runtime["executable"] && digest(actual)==runtime["executable_sha256"] &&
        isabspath(runtime["launcher"]) && digest(runtime["launcher"])==runtime["launcher_sha256"],"wrong Julia launcher or actual runtime")
    run=RecoveryRun(RP.RecoveryRunPaths(m["run"],m["retained_jobs"]),state.contexts,provenance(m,sha_arg))
    (;data=m,path,sha=sha_arg,frozen,predecessor,run)
end

function require_new_output(m,output;canonical=nothing)
    require(output isa String && abspath(output)==output && normpath(output)==output && RP.inside(output,root(m.run)) &&
        output!=root(m.run) && realpath(dirname(output))==dirname(output) && !ispath(output) && !islink(output) &&
        !ispath(output*".progress.json") && !islink(output*".progress.json"),"fresh output inside the current run required")
    require(canonical===nothing || output==canonical,"wrong canonical new job output")
    output
end

function stable(m)
    require(JSON.parsefile(m.path)==m.data,"prepared recovery manifest changed during its owned lifetime")
    load_manifest(m.frozen,m.predecessor,m.path,m.sha)
    true
end
write_record(path,record)=CP.write_record(path,record)
end
