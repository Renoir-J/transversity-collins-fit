"Explicit synthetic numerical rows and saved-protocol fixtures; no dispatch loader."
module SIACompleteJobFixture
using CollinsFit, JSON, SHA, Dates
include("../bin/SIACompleteJobs.jl")
const C=SIACompleteJobs
const M=C.N.M
const ROOT=dirname(dirname(@__DIR__))
const FIXTURE="explicitly_synthetic_complete72_saved_job_fixture_not_dispatch"
include("sia_geometry_fixtures.jl")

function create(directory;force4=false)
    full=geometry_fixture_request(directory)
    parent=JSON.parsefile(joinpath(ROOT,"config/historical-sia-integral-contract.json"))
    markers=Dict{String,String}()
    for name in ("adapter","frozen","build","controls",("retained$i" for i in 1:9)...)
        path=joinpath(directory,"fixture-sources",name,"identity.json");mkpath(dirname(path))
        writejson(path,Dict("fixture_scope"=>FIXTURE,"kind"=>name,"native_or_experimental_input"=>false))
        markers[name]=path
    end
    owner=joinpath(directory,"synthetic-owner.json")
    writejson(owner,Dict("fixture_scope"=>FIXTURE,"live_process_claimed"=>false))
    planfile=joinpath(directory,"plan.json")
    plan=Dict{String,Any}("fixture_scope"=>FIXTURE,"adapter"=>dirname(markers["adapter"]),"frozen"=>dirname(markers["frozen"]),
        "adapter_files"=>Dict("identity.json"=>C.digest(markers["adapter"])),"frozen_files"=>Dict("identity.json"=>C.digest(markers["frozen"])),
        "request"=>full.path,"source_request"=>full.path,"source_inputs"=>Dict{String,String}())
    writejson(planfile,plan)
    m=Dict{String,Any}("schema_version"=>1,"status"=>"EXPLICIT_SYNTHETIC_CONTEXT_NOT_A_DISPATCH_MANIFEST","fixture_scope"=>FIXTURE,
        "theory_scope"=>"synthetic_closure","run"=>directory,"force_native4_fixture"=>force4,
        "julia_launcher"=>first(Base.julia_cmd().exec),"resources"=>JSON.parsefile(joinpath(ROOT,"config/historical-sia-job-contract.json"))["resources"],
        "sidis_owner_file"=>owner,"sidis_owner_file_sha256"=>C.digest(owner),
        "build_inputs"=>Dict(markers["build"]=>C.digest(markers["build"])),
        "source_control_files"=>Dict(markers["controls"]=>C.digest(markers["controls"])),
        "retained_files"=>Dict(markers["retained$i"]=>C.digest(markers["retained$i"]) for i in 1:9),
        "limitations"=>["All rows and command/resource records are explicit synthetic fixtures.",
            "No native evolution, raw experimental intake, production resource lock or physical dispatch occurs."])
    source=(request=full,parent=parent,plan=plan,planfile=planfile,plan_sha=C.digest(planfile),frozen=plan["frozen"])
    partial=(source=source,m=m)
    inputs=C.plan_inputs(partial);inputs[owner]=C.digest(owner)
    for key in ("build_inputs","source_control_files","retained_files");merge!(inputs,m[key]);end
    m["execution_inputs"]=inputs
    path=joinpath(directory,"execution.json");writejson(path,m)
    load(directory)
end

function load(directory)
    path=joinpath(directory,"execution.json");m=JSON.parsefile(path)
    C.require(m["fixture_scope"]==FIXTURE && m["theory_scope"]=="synthetic_closure","explicit fixture scope required")
    planfile=joinpath(directory,"plan.json");plan=JSON.parsefile(planfile)
    data=JSON.parsefile(plan["request"])
    C.require(M.S.declaration(data)=="synthetic_closure","fixture cannot relabel primary data")
    request=HistoricalSIAIntegralRequest(data,plan["request"],C.digest(plan["request"]),CollinsFit.numerical_snapshot(data))
    sha=C.digest(path)
    source=(request=request,parent=JSON.parsefile(joinpath(ROOT,"config/historical-sia-integral-contract.json")),
        plan=plan,planfile=planfile,plan_sha=C.digest(planfile),frozen=plan["frozen"])
    binding=Dict{String,Any}("fixture_scope"=>FIXTURE,"execution_manifest_sha256"=>sha,"request_sha256"=>request.sha256,
        "native_evolution_executed"=>false,"raw_experimental_inputs_admitted"=>false)
    (source=source,m=m,path=path,sha=sha,binding=binding)
end

function seal(c,phase,index;verify=true)
    root=C.folder(c,phase,index);budget=c.m["resources"][phase in C.PHYSICAL ? phase : "global"]
    reportpath=joinpath(root,"report.json");auditpath=joinpath(root,"audit.json")
    audit=JSON.parsefile(auditpath);audit["report_sha256"]=C.digest(reportpath);writejson(auditpath,audit)
    steps=Any[]
    for name in ("producer","audit")
        log=joinpath(root,name*".log");resource=joinpath(root,name*".resources.txt")
        isfile(log) || write(log,"Explicit synthetic protocol fixture; per-job command not executed.\n")
        isfile(resource) || write(resource,"COLLINSFLOW_TIME_V1\t100\t1.0\t0.5\t0.1\t0\n")
        push!(steps,Dict{String,Any}("name"=>name,"command"=>C.command(c,phase,index,root,name;report_sha=C.digest(reportpath)),
            "returncode"=>0,"process_returncode"=>0,"interruption"=>nothing,"termination_signal"=>nothing,
            "wrapper_elapsed_seconds"=>1.1,"minimum_available_memory_observed_before_command"=>2^40,
            "resource"=>C.resource(resource),"resource_sha256"=>C.digest(resource),"log_sha256"=>C.digest(log),
            "fixture_scope"=>"synthetic_measurement_shape_not_actual_job_measurement"))
    end
    initial=JSON.parsefile(joinpath(root,"initial-invocation.json"))
    record=merge(deepcopy(initial),Dict{String,Any}("status"=>C.COMPLETE,"record_role"=>"complete_job_execution_invocation",
        "fixture_scope"=>FIXTURE,"per_job_commands_executed"=>false,"completed_at"=>string(now(UTC)),"steps"=>steps,
        "initial_record_sha256"=>C.digest(joinpath(root,"initial-invocation.json")),"resources"=>budget,
        "inputs"=>C.job_inputs(c,phase,index,C.base_inputs(c)),
        "outputs"=>Dict(joinpath(root,name)=>C.digest(joinpath(root,name)) for name in C.PRODUCTS)))
    writejson(joinpath(root,"invocation.json"),record)
    write(joinpath(root,"job.completed"),C.digest(joinpath(root,"invocation.json"))*"\n")
    verify && C.verify_pair(c,phase,index,C.base_inputs(c))
end

function save_job(c,phase,index,report,audit)
    root=C.folder(c,phase,index);C.require(!ispath(root),"retain occupied synthetic job directory")
    mkpath(root)
    writejson(joinpath(root,"lease.json"),Dict("fixture_scope"=>FIXTURE,"live_production_lease_claimed"=>false))
    initial=merge(C.N.flags(),Dict{String,Any}("status"=>"FAIL","scope"=>"source_bound_complete_sia_Julia_job_pair",
        "theory_scope"=>"synthetic_closure","fixture_scope"=>FIXTURE,"record_role"=>"initial_sentinel_not_an_executed_failure",
        "phase"=>phase,"index"=>index,"manifest_sha256"=>c.sha,"started_at"=>string(now(UTC)),"completed_at"=>nothing,
        "steps"=>Any[],"resources"=>c.m["resources"][phase in C.PHYSICAL ? phase : "global"]))
    writejson(joinpath(root,"initial-invocation.json"),initial)
    report["completed_at"]=string(now(UTC));audit["completed_at"]=string(now(UTC))
    writejson(joinpath(root,"report.json"),report);writejson(joinpath(root,"audit.json"),audit)
    seal(c,phase,index)
end

function produce_all(c,legacy)
    data=c.source.request.data;parent=c.source.parent
    groups,keys=M.S.geometry_groups(data)
    contexts=Dict(key=>geometry_fixture_context(key[1],key[2]) for key in keys)
    fine=Dict(key=>polarized_leptonic_plan(base.leptonic,1.01.*base.unpolarized_first,base.unpolarized_second,
        base.first_operator,base.first_samples;charges=(1,-1)) for (key,base) in contexts)
    function prepare(view,parent,level,B)
        key=historical_sia_integral_geometry_key(first(view.data["points"]))
        context=c.m["force_native4_fixture"] && level>=3 ? fine[key] : contexts[key]
        context,Dict("geometry"=>collect(key),"bound_inverse_GeV"=>B,"radial_nodes"=>length(context.leptonic.radial.b),
            "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes,
            "cold_seconds_informational"=>0.,"context_provenance"=>deepcopy(context.provenance))
    end
    jobs=Tuple{String,Int}[]
    function complete(phase,index)
        proof=C.collect_dependencies(c,phase,legacy);records=proof.records
        report=C.common_header(c,phase,index,C.COMPUTED);audit=C.common_header(c,phase,index,C.REVIEWED)
        if phase in C.PHYSICAL
            spec=C.definition(c,phase,index,records);view=M.SelectedSIAGeometry(c.source.request,index)
            leaf=C.compute_leaf(view,parent,spec,c.binding,legacy,prepare)
            report["leaf_report"]=leaf
            audit["leaf_arithmetic_review"]=C.audit_leaf(data,parent,spec,c.binding,leaf,legacy)
        else
            merge!(report,C.produce_global(c,phase,records));C.audit_global(c,phase,records,report)
            audit["complete72_saved_arithmetic_verified"]=(phase=="diagnostic")
        end
        audit["complete_prerequisite_pair_ledgers_verified"]=length(records)
        for key in ("native_evolution_reexecuted_by_reader","cold_AD_or_integrals_reexecuted_by_reader","raw_primary_tables_replayed_in_Julia")
            audit[key]=false
        end
        C.recheck_dependencies(proof)
        save_job(c,phase,index,report,audit);push!(jobs,(phase,index))
        println("Explicit synthetic saved-job fixture: ",phase," ",index);flush(stdout)
    end
    for i in 1:36;complete("native123",i);end
    complete("native-select3",0)
    if C.needs4(c);for i in 1:36;complete("native4",i);end;end
    complete("native-selection",0)
    for i in 1:36;complete("tailref",i);end
    complete("tail-selection",0)
    for i in (1,15,36);complete("derivative",i);end
    complete("diagnostic",0)
    writejson(joinpath(c.m["run"],"fixture-jobs.json"),Dict("fixture_scope"=>FIXTURE,"jobs"=>[collect(job) for job in jobs]))
    jobs
end

function build(directory;force4=false)
    c=create(directory;force4=force4)
    legacy=C.load_library(ROOT)
    jobs=Base.invokelatest(produce_all,c,legacy)
    c,jobs
end

function inspect_loaded(c,legacy)
    report=C.verify_pair(c,"diagnostic",0,C.base_inputs(c))
    proof=C.collect_dependencies(c,"diagnostic",legacy)
    C.audit_global(c,"diagnostic",proof.records,report);C.recheck_dependencies(proof)
    Dict("fixture_scope"=>FIXTURE,"theory_scope"=>"synthetic_closure","full72_saved_fixture_arithmetic_verified"=>true,
        "prerequisite_pairs"=>length(proof.records),"selected_native_level"=>report["selected_native_level"],
        "source_point_ids"=>[p["id"] for p in report["points"]],"per_job_commands_measured"=>false,
        "production_locks_acquired"=>false,"native_evolution_executed"=>false,"physical_SIA_jobs_started"=>false,
        "fit_performed"=>false,"production_ready"=>false,"accepted_fitting_gate"=>false)
end

function inspect(directory)
    c=load(directory);legacy=C.load_library(ROOT)
    Base.invokelatest(inspect_loaded,c,legacy)
end
end
