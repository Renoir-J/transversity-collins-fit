"Owned locations of retained jobs; no dispatch or numerical admission."
module SIDISRecoveryPaths
using JSON, SHA
const ROOT=dirname(dirname(@__DIR__))
const POLICY_SHA="b9b704ab85b473021de062ef555a1979d2522a9760de1122ffb3fa350b4bbcf2"
const POINT_PHASES=("native123","native4","tailref123","reference4")
const GLOBAL_PHASES=("native-select3","native-selection","source-math","derivative-select3","diagnostic")
const FIELDS=("phase","index","folder","producer_root","producer_manifest","producer_manifest_sha256","protocol")
const PROTOCOLS=("frozen_original","continuation_v1","recovery_v2")
const FLAGS=("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted",
    "fit_performed","accepted_prediction_diagnostic","accepted_fitting_gate","G2PERF_pass","interval_coverage_verified","all132_coverage_verified")
require(value,message)=value || error(message)
digest(path)=bytes2hex(sha256(read(path)))

function policy()
    path=joinpath(ROOT,"config/sidis-recovery-contract.json")
    require(digest(path)==POLICY_SHA,"prospective recovery policy changed")
    data=JSON.parsefile(path)
    require(data["retained_jobs"]["point_phases"]==collect(POINT_PHASES) && all(get(data,key,nothing)===false for key in FLAGS),
        "wrong or self-admitting recovery policy")
    data
end

function canonical(path;directory)
    require(path isa AbstractString && isabspath(path) && normpath(path)==path &&
        !any(occursin(character,path) for character in ('\0','\n','\r')) &&
        (directory ? isdir(path) : isfile(path)) && !islink(path) && realpath(path)==path,"ordinary canonical file or directory required")
    String(path)
end

function point_index(phase,index)
    require(phase in POINT_PHASES && index isa Int && 1<=index<=132,"declared point phase and exact index1..132 required")
end

inside(path,root)=path==root || startswith(path,root*"/")

function within_directory(path,root)
    target=stat(root);identity=(target.device,target.inode);current=path
    while true
        value=stat(current)
        (value.device,value.inode)==identity && return true
        parent=dirname(current);parent==current && return false;current=parent
    end
end

function entry(root,row)
    require(row isa AbstractDict && Set(keys(row))==Set(FIELDS),"complete location schema required; no admission fields allowed")
    point_index(row["phase"],row["index"])
    require(row["protocol"] in PROTOCOLS,"unsupported producing protocol")
    producer=canonical(row["producer_root"];directory=true);folder=canonical(row["folder"];directory=true)
    manifest=canonical(row["producer_manifest"];directory=false)
    require(folder==joinpath(producer,row["phase"],"p"*lpad(string(row["index"]),4,'0')),"folder does not match its exact producing phase/index")
    protected=row["protocol"]=="frozen_original" ? dirname(producer) : producer
    require(!(within_directory(root,protected) || within_directory(protected,root)),"current and protected producing roots must be disjoint")
    if row["protocol"]=="frozen_original"
        require(row["phase"]=="native123" && row["index"]<=13 && manifest==joinpath(dirname(producer),"invocation.json"),
            "only the declared original13 native123 folders belong to the original anchor")
    else
        require(manifest==joinpath(producer,row["protocol"]=="continuation_v1" ? "continuation.json" : "recovery.json"),
            "manifest does not belong to the declared producing root/protocol")
    end
    sha=row["producer_manifest_sha256"]
    require(sha isa AbstractString && occursin(r"^[0-9a-f]{64}$",sha) && digest(manifest)==sha,"changed externally bound producing manifest")
    deepcopy(Dict{String,Any}(row))
end

signature(root,entries)=bytes2hex(sha256(JSON.json(Any[root,[[row[key] for key in FIELDS] for row in entries]])))

struct RecoveryRunPaths
    root::String
    entries::Vector{Dict{String,Any}}
    signature::String
    function RecoveryRunPaths(root,entries)
        policy();root=canonical(root;directory=true)
        require(entries isa AbstractVector,"ordered retained-job entries required")
        owned=Dict{String,Any}[entry(root,row) for row in entries]
        keys=[(findfirst(==(row["phase"]),POINT_PHASES),row["index"]) for row in owned]
        require(keys==sort(unique(keys)),"unique source-ordered retained phase/index entries required")
        producers=Dict{String,Any}()
        for row in owned
            producer=row["producer_root"];identity=(row["protocol"],row["producer_manifest"],row["producer_manifest_sha256"])
            protected=row["protocol"]=="frozen_original" ? dirname(producer) : producer
            if haskey(producers,producer)
                require(producers[producer][1]==identity,"one producing root cannot claim different protocols/manifests")
            else
                require(!any(within_directory(protected,old[2]) || within_directory(old[2],protected) for old in values(producers)),
                    "protected producing roots must be disjoint")
                producers[producer]=(identity,protected)
            end
        end
        new(root,owned,signature(root,owned))
    end
end

function valid(run::RecoveryRunPaths)
    policy()
    require(signature(run.root,run.entries)==run.signature && canonical(run.root;directory=true)==run.root,"mutated owned recovery path view")
    for row in run.entries;entry(run.root,row);end
    run
end

function retained(run::RecoveryRunPaths,phase,index)
    point_index(phase,index);valid(run)
    rows=filter(row->row["phase"]==phase && row["index"]==index,run.entries)
    isempty(rows) ? nothing : deepcopy(only(rows))
end

function current_path(run::RecoveryRunPaths,parts::AbstractString...)
    valid(run);path=joinpath(run.root,parts...)
    require(normpath(path)==path && inside(path,run.root),"current-run paths cannot escape")
    ancestor=path
    while !ispath(ancestor) && !islink(ancestor);ancestor=dirname(ancestor);end
    require(!islink(ancestor) && realpath(ancestor)==ancestor,"current-run paths cannot traverse aliases")
    path
end
Base.joinpath(run::RecoveryRunPaths,parts::AbstractString...)=current_path(run,parts...)

function folder(run::RecoveryRunPaths,phase,index=0)
    if phase in POINT_PHASES
        row=retained(run,phase,index)
        return row===nothing ? current_path(run,phase,"p"*lpad(string(index),4,'0')) : row["folder"]
    end
    require(phase in GLOBAL_PHASES && index isa Int && index==0,"declared fresh global phase with index0 required")
    current_path(run,phase)
end

function fresh_folder(run::RecoveryRunPaths,phase,index=0)
    path=folder(run,phase,index)
    require(inside(path,run.root) && !ispath(path) && !islink(path),"never dispatch through retained, aliased or occupied folders")
    path
end

function metadata(run::RecoveryRunPaths)
    valid(run)
    merge(Dict{String,Any}(key=>false for key in FLAGS),Dict{String,Any}("scope"=>"owned_whole_job_path_view_not_admission",
        "root"=>run.root,"retained_jobs"=>deepcopy(run.entries),"policy_sha256"=>POLICY_SHA,
        "numerical_dispatch_ready"=>false,"complete_job_provenance_verified"=>false,"saved_numerical_arithmetic_verified"=>false))
end
end
