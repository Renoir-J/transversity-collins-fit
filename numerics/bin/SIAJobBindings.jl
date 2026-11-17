"Verify source-plan bindings before including frozen numerical helper libraries."
module SIAJobBindings
using CollinsFit, JSON, SHA
const ADAPTER=dirname(dirname(@__DIR__))
const FROZEN_HEAD="6c88d12516074889f49a0dfeae456d6cc63491b4"
const NUMERICAL_FP="01be55cc1449679f6dc26bd43ea9788604938c1fecbb25d276acaa35bc4bae47"
const POLICY_HEAD="655eee588b7bf8c69e120f3bb62fb707b794dfa2"
const POLICY_SHA="da14c5b667c2866aea2aebe546baedb45d97ad7d513c6bed1b1aa884a23bf337"
const PARENT_SHA="29830e78c4e48d2eb8e5a70711885fc961b2b94e0ec3303315befd59a65fdf97"
const FLAGS=("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted",
    "fit_performed","G2PERF_pass","interval_coverage_verified","accepted_fitting_gate","numerical_dispatch_ready","physical_jobs_started")
require(value,message)=value || error(message)
digest(path)=bytes2hex(sha256(read(path)))

function checkout(root,head,inventory)
    require(realpath(root)==root && readchomp(`git -C $root rev-parse HEAD`)==head,"wrong canonical source checkout")
    actual=Dict{String,String}()
    for entry in split(read(`git -C $root ls-tree -r -z HEAD`,String),'\0';keepempty=false)
        metadata,path=split(entry,'\t';limit=2);mode,kind,blob=split(metadata)
        require(mode in ("100644","100755") && kind=="blob","ordinary committed source files required")
        filename=joinpath(root,path)
        require(realpath(filename)==filename && isfile(filename),"source aliases or missing files are not allowed")
        bytes=read(filename)
        require(bytes2hex(sha1(vcat(Vector{UInt8}(codeunits("blob $(length(bytes))\0")),bytes)))==blob,"raw source blob changed")
        actual[path]=bytes2hex(sha256(bytes))
    end
    require(actual==inventory,"missing or changed complete source inventory")
    require(isempty(read(`git -C $root ls-files --others --exclude-standard -z`,String)),"commit complete source before running controls")
    actual
end

function load(frozen,planfile,expected_sha)
    frozen,planfile=abspath(frozen),abspath(planfile)
    require(realpath(frozen)==frozen && realpath(planfile)==planfile && digest(planfile)==expected_sha,"changed externally bound source plan")
    plan=JSON.parsefile(planfile)
    require(plan["schema_version"]===1 && plan["status"]=="SOURCE_PLAN_PREPARED_NOT_NUMERICAL_DISPATCH_OR_GATE" &&
        plan["scope"]=="historical_reconstructed_sia_complete_geometry_jobs" && plan["theory_scope"]=="reconstructed_theory" &&
        all(get(plan,key,nothing)===false for key in FLAGS),"source plan cannot authorize physical work")
    require(plan["frozen"]==frozen && plan["frozen_head"]==FROZEN_HEAD && plan["adapter"]==ADAPTER &&
        plan["numerical_source_fingerprint"]==NUMERICAL_FP && plan["policy_commit"]==POLICY_HEAD &&
        plan["policy_sha256"]==POLICY_SHA && plan["adapter_head"]!=POLICY_HEAD,"wrong source-plan origin")
    require(startswith(planfile,ADAPTER*"/runs/") && plan["request"]==joinpath(dirname(planfile),"request.json"),"wrong source-plan namespace")
    require(dirname(dirname(dirname(pathof(CollinsFit))))==frozen,"use the explicit frozen Julia project")
    checkout(frozen,FROZEN_HEAD,plan["frozen_files"]);checkout(ADAPTER,plan["adapter_head"],plan["adapter_files"])
    require(length(plan["frozen_files"])==302 && all(get(plan["adapter_files"],name,nothing)==sha for (name,sha) in plan["frozen_files"]),"all302 inherited numerical files must stay unchanged")
    require(success(`git -C $ADAPTER merge-base --is-ancestor $POLICY_HEAD $(plan["adapter_head"])`),"prospective policy must precede implementation")
    require(digest(joinpath(ADAPTER,"config/historical-sia-job-contract.json"))==POLICY_SHA &&
        digest(joinpath(frozen,"config/historical-sia-integral-contract.json"))==PARENT_SHA,"source or scientific policy changed")
    request=load_historical_sia_integral_request(plan["request"])
    require(request.sha256==plan["request_sha256"] && request.data["source_fingerprint"]==NUMERICAL_FP &&
        request.data["source_inputs"]==plan["source_inputs"] && [point["id"] for point in request.data["points"]]==plan["point_order"],"changed full72 source request")
    require(digest(plan["source_request"])==plan["source_request_sha256"],"changed original source request bytes")
    parent=JSON.parsefile(joinpath(frozen,"config/historical-sia-integral-contract.json"))
    require(VERSION==v"1.12.6","use the validated Julia1.12.6 runtime")
    executable=realpath("/proc/self/exe")
    binding=Dict{String,Any}("source_plan_sha256"=>expected_sha,"source_request_sha256"=>request.sha256,
        "frozen_head"=>FROZEN_HEAD,"numerical_source_fingerprint"=>NUMERICAL_FP,
        "adapter_head"=>plan["adapter_head"],"adapter_files"=>plan["adapter_files"],
        "parent_contract_sha256"=>PARENT_SHA,"job_policy_sha256"=>POLICY_SHA,
        "Julia_version"=>string(VERSION),"Julia_executable"=>executable,"Julia_executable_sha256"=>digest(executable),
        "source_inputs"=>request.data["source_inputs"],"raw_primary_replay"=>"Python source replay attestation; Julia verifies current full source bindings")
    (request=request,parent=parent,binding=binding,plan=plan,frozen=frozen,planfile=planfile,plan_sha=expected_sha)
end

function recheck(context)
    require(digest(context.planfile)==context.plan_sha && digest(context.plan["request"])==context.request.sha256 &&
        digest(context.plan["source_request"])==context.plan["source_request_sha256"],"source request or plan bytes changed during controls")
    CollinsFit.historical_sia_integral_request_valid(context.request)
    checkout(context.frozen,FROZEN_HEAD,context.plan["frozen_files"])
    checkout(ADAPTER,context.plan["adapter_head"],context.plan["adapter_files"])
    require(digest(context.planfile)==context.plan_sha && digest(context.plan["request"])==context.request.sha256 &&
        digest(context.plan["source_request"])==context.plan["source_request_sha256"],"source bytes changed during final scan")
    context
end

function fresh_output(path,planfile)
    path=abspath(path);run=dirname(abspath(planfile))
    require(normpath(path)==path && startswith(path,run*"/") && !(ispath(path) || islink(path)),"fresh output below this source-plan run required")
    ancestor=dirname(path)
    while !ispath(ancestor)
        ancestor=dirname(ancestor)
    end
    require(realpath(ancestor)==ancestor,"output ancestors cannot be aliases")
    mkpath(dirname(path));require(realpath(dirname(path))==dirname(path),"output parent cannot be aliased")
    require(!(ispath(path*".writing") || islink(path*".writing")),"retain occupied partial outputs")
    path
end

function atomic_report(path,data)
    temporary=path*".writing";writejson(temporary,data);mv(temporary,path;force=true)
end
end
