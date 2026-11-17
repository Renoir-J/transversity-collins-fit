module SIDISContinuationPaths
using JSON, SHA, Dates

const ADAPTER=dirname(dirname(@__DIR__))
const FROZEN_HEAD="bec12db0d65ed67a74de80bca949e21b80a99757"
const CONTRACT_SHA="47059f3236115d0c171bb7b487929b13ca579d876e61839cb33632f0a62a2417"
const POLICY_COMMIT="e9f86df76ab5f1bedbea52290a0aad5f70be7a2e"
const SCOPE="same_attempt_reconstructed_full132_continuation_not_fitting"
const FLAGS=("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted",
    "fit_performed","accepted_prediction_diagnostic","accepted_fitting_gate","G2PERF_pass",
    "interval_coverage_verified","all132_coverage_verified")
const POINT_PHASES=("native123","native4","tailref123","reference4")
const GLOBAL_PHASES=("native-select3","native-selection","source-math","derivative-select3","diagnostic")
digest(path)=bytes2hex(sha256(read(path)))
false_flags()=Dict{String,Any}(k=>false for k in FLAGS)

"A new argument type; existing String path methods and frozen source stay intact."
struct ContinuationRun
    root::String
    original::String
    function ContinuationRun(root::AbstractString,original::AbstractString)
        r=abspath(root);o=abspath(original)
        r==root && o==original && realpath(r)==r && realpath(o)==o || error("canonical real run roots required")
        (r==o || startswith(r,o*"/") || startswith(o,r*"/")) && error("run roots must be disjoint")
        new(r,o)
    end
end
Base.joinpath(run::ContinuationRun,paths::AbstractString...)=joinpath(run.root,paths...)

function point_folder(run::ContinuationRun,phase,index)
    phase in POINT_PHASES && index isa Integer && !(index isa Bool) && 1<=index<=132 || error("exact phase/index required")
    base=phase=="native123" && index<=13 ? run.original : run.root
    joinpath(base,phase,"p"*lpad(string(index),4,'0'))
end

function verify_git(root,head)
    readchomp(`git -C $root rev-parse HEAD`)==head || error("wrong committed source head")
    isempty(read(`git -C $root status --porcelain`,String)) || error("dirty source checkout")
    files=Dict{String,String}()
    for entry in split(read(`git -C $root ls-tree -rz HEAD`,String),'\0';keepempty=false)
        meta,name=split(entry,'\t';limit=2);_,kind,oid=split(meta)
        kind=="blob" && !occursin('\n',name) && !occursin('\r',name) || error("non-file Git source")
        path=joinpath(root,name)
        isfile(path) && !islink(path) || error("missing/aliased source")
        bytes=read(path)
        bytes2hex(sha1(vcat(Vector{UInt8}(codeunits("blob $(length(bytes))\0")),bytes)))==oid || error("working bytes differ from Git: "*name)
        files[name]=bytes2hex(sha256(bytes))
    end
    files
end

function verify_inventory(original,inventory)
    expected=Dict{String,Any}()
    for row in inventory
        p=row["path"]
        startswith(p,original*"/") && abspath(p)==p && realpath(p)==p && !haskey(expected,p) || error("invalid original inventory")
        expected[p]=row
    end
    current=Set(joinpath(folder,name) for (folder,_,files) in walkdir(original) for name in files)
    current==Set(keys(expected)) || error("original file inventory changed")
    for (path,row) in expected
        filesize(path)==row["bytes"] && digest(path)==row["sha256"] || error("original file changed")
    end
    true
end

function load_manifest(frozen_arg,path_arg,sha_arg)
    frozen=abspath(frozen_arg);path=abspath(path_arg)
    frozen==frozen_arg && realpath(frozen)==frozen && digest(path)==sha_arg || error("changed explicit continuation inputs")
    m=JSON.parsefile(path);cpath=joinpath(ADAPTER,"config/sidis-continuation-contract.json")
    digest(cpath)==CONTRACT_SHA || error("changed prospective continuation contract")
    c=JSON.parsefile(cpath)
    m["status"]=="PREPARED_IMMUTABLE_CONTINUATION_NOT_A_GATE" && m["scope"]==SCOPE &&
        m["contract_sha256"]==CONTRACT_SHA && m["frozen"]==frozen && m["adapter"]==ADAPTER &&
        m["frozen_head"]==FROZEN_HEAD && m["numerical_source_fingerprint"]==c["frozen_numerical_fingerprint"] &&
        all(get(m,k,nothing)===false && get(c,k,nothing)===false for k in FLAGS) || error("wrong/self-admitting continuation manifest")
    frozen_files=verify_git(frozen,FROZEN_HEAD)
    adapter_files=verify_git(ADAPTER,m["adapter_commit"])
    length(frozen_files)==335 && frozen_files==m["frozen_files"] && adapter_files==m["adapter_files"] &&
        all(get(adapter_files,k,nothing)==v for (k,v) in frozen_files) || error("changed frozen or adapter source bindings")
    success(`git -C $ADAPTER merge-base --is-ancestor $POLICY_COMMIT $(m["adapter_commit"])`) &&
        m["adapter_commit"]!=POLICY_COMMIT || error("implementation must follow frozen addendum")
    original=joinpath(frozen,"runs/historical-sidis-jobs-real-data/attempt-2")
    m["original"]==original && m["run"]==dirname(path) && startswith(m["run"],joinpath(frozen,"runs")*"/") || error("wrong original or continuation root")
    origin=c["supported_origin"]
    for (key,field) in (("review","independent_prefix_review_sha256"),("integrity","post_review_integrity_sha256"),("loss","original_loss_inventory_sha256"))
        digest(m[key])==origin[field] || error("changed pinned original evidence")
    end
    workspace=dirname(dirname(dirname(frozen)))
    loss=JSON.parsefile(m["loss"])
    inventory=[merge(row,Dict("path"=>joinpath(workspace,row["path"]))) for row in loss["interrupted_attempts"]["sidis132"]["retained_files"]]
    inventory==m["original_inventory"] && length(inventory)==141 || error("exact original141-file inventory required")
    verify_inventory(original,inventory)
    digest(joinpath(original,"invocation.json"))==origin["outer_invocation_sha256"] || error("changed original outer invocation")
    deps=JSON.parsefile(joinpath(original,"invocation.json"))["dependency_receipts"]
    deps==m["dependency_receipts"] && Set(keys(deps))==Set((m["native"],m["parton"])) &&
        all(digest(p)==h for (p,h) in deps) || error("changed original build receipts")
    m["request"]==joinpath(original,"full/request.json") && digest(m["request"])==m["request_sha256"] || error("changed original request")
    run=ContinuationRun(m["run"],joinpath(original,"full"))
    (;data=m,path,sha=sha_arg,frozen,run)
end

provenance(m)=Dict("scope"=>SCOPE,"manifest_sha256"=>m.sha,"adapter_commit"=>m.data["adapter_commit"],"contract_sha256"=>CONTRACT_SHA)

function require_new_output(m,output;canonical=nothing)
    p=abspath(output)
    p==output && startswith(p,m.run.root*"/") && realpath(dirname(p))==dirname(p) &&
        !ispath(p) && !ispath(p*".progress.json") || error("fresh output inside the new run required")
    canonical===nothing || p==canonical || error("wrong canonical new job output")
    p
end

function stable(m)
    load_manifest(m.frozen,m.path,m.sha)
    true
end

function write_record(path,record)
    temporary=path*".tmp"
    open(temporary,"w") do io
        JSON.print(io,record,2);println(io)
    end
    mv(temporary,path)
end
end
