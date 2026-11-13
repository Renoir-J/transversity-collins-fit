# Prospective comparator only. It does not select a production coupling or reevaluate physical moments.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra
include("NativeCouplingReference.jl")
using .NativeCouplingReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/native-coupling-contract.json")
const SOURCES=vcat([@__FILE__,joinpath(@__DIR__,"NativeCouplingReference.jl"),CONTRACT,
    joinpath(ROOT,"fortran/native_coupling_reference.f90"),joinpath(ROOT,"fortran/hoppet_bridge.f90"),
    joinpath(ROOT,"src/collinsflow/hoppet.py"),joinpath(ROOT,"src/collinsflow/provenance.py"),
    joinpath(ROOT,"scripts/prepare_native_coupling.py"),joinpath(ROOT,"scripts/prepare_hoppet.py"),
    joinpath(ROOT,"config/hoppet-contract.json"),joinpath(ROOT,"config/evolution-contract.json"),
    joinpath(ROOT,"config/css-contract.json"),joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in SOURCES)

function coupling_checks!(report,buildfile,output)
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"];settings=contract["settings"]
    hashes=source_hashes();build=CollinsFit.verified_hoppet_build(buildfile)
    buildhash=bytes2hex(sha256(read(buildfile)));contracthash=bytes2hex(sha256(read(CONTRACT)))
    reference=build["coupling_reference"];exe=reference["executable"]
    reference["protocol"]==contract["native_reference"]["protocol"] && reference["contract_sha256"]==contracthash ||
        error("reference protocol/contract binding mismatch")
    reference["paper_sha256"]==contract["primary_source"]["sha256"] || error("reference paper binding mismatch")
    paperfile=joinpath(dirname(buildfile),"paper.html")
    bytes2hex(sha256(read(paperfile)))==contract["primary_source"]["sha256"] || error("locked primary paper changed")
    haskey(build["output_hashes"],exe) && isfile(exe) || error("reference executable lacks output hash")
    for path in ("src/collinsflow/hoppet.py","src/collinsflow/provenance.py","fortran/hoppet_bridge.f90",
        "fortran/native_coupling_reference.f90","scripts/prepare_hoppet.py","scripts/prepare_native_coupling.py",
        "config/native-coupling-contract.json","config/hoppet-contract.json","config/evolution-contract.json",
        "config/matched-rows-contract.json","config/hoppet-batch-contract.json")
        matches=filter(p->endswith(replace(p,'\\'=>'/'),"/"*path),collect(keys(build["source_hashes"])))
        length(matches)==1 && build["source_hashes"][only(matches)]==bytes2hex(sha256(read(joinpath(ROOT,path)))) ||
            error("reference build does not bind current project source: $path")
    end
    report["sources"]=hashes;report["contract_id"]=contract["id"];report["contract_sha256"]=contracthash
    report["native_build"]=build;report["build_file"]=buildfile;report["build_sha256"]=buildhash
    report["primary_source"]=contract["primary_source"];report["limitations"]=contract["limitations"]
    report["settings"]=settings;report["policy"]=policy
    checks=report["checks"]
    check(name,a,b,tol;absolute=false)=push!(checks,comparison(name,a,b,tol;floor=policy["relative_floor"],absolute=absolute))
    q0=settings["Q02_GeV2"];nf=settings["nf"]
    probes=vcat(Float64.(settings["Q2_probes_GeV2"]),Float64(settings["MZ_probe_GeV"])^2)
    probes[1]==q0 && length(unique(probes))==4 || error("four distinct forward probes including Q0 required")
    report["Q2_probes_GeV2"]=probes
    coupling=AsymptoticTwoLoopCoupling(settings["Lambda_GeV"],nf)
    baseline=Dict("id"=>contract["models"][1]["id"],"alpha0"=>alpha_s(coupling,q0),
        "alpha"=>alpha_s.(Ref(coupling),probes),"tau"=>evolution_time.(Ref(coupling),probes,q0))
    report["asymptotic_baseline"]=baseline
    directory=joinpath(dirname(output),"batches");index=Ref(0)
    cache=Dict{Tuple{Float64,Float64,Float64,Int,Int},Any}()
    function independent(a,q,qstart,flavors,loops)
        key=(Float64(a),Float64(q),Float64(qstart),Int(flavors),Int(loops))
        get!(cache,key) do
            rg_reference(a,q,qstart,flavors,loops;bits=contract["independent_reference"]["precision_bits"])
        end
    end
    function trial(a,qstart,qs,flavors,loops,dt,n)
        rules=[log_rule(qstart,q,n) for q in qs]
        queries=sort(unique(vcat(qs,[v for r in rules for v in r.nodes])))
        index[]+=1
        native=native_batch(exe,a,qstart,queries,flavors,loops,dt,directory,index[])
        push!(report["native_batches"],native.evidence)
        refs=[independent(a,q,qstart,flavors,loops) for q in queries]
        alpha_refs=[v.alpha for v in refs]
        probe_indices=searchsortedfirst.(Ref(queries),qs)
        alpha=native.values[probe_indices];reference_alpha=alpha_refs[probe_indices]
        times=[isempty(r.nodes) ? 0. : dot(r.weights,native.values[searchsortedfirst.(Ref(queries),r.nodes)]) for r in rules]
        reference_times=[independent(a,q,qstart,flavors,loops).tau for q in qs]
        Dict{String,Any}("dt"=>dt,"quadrature_order"=>n,"Q02_GeV2"=>qstart,"alpha0"=>a,"nf"=>flavors,"nloop"=>loops,
            "query_Q2_GeV2"=>queries,"native_alpha"=>native.values,"reference_alpha"=>alpha_refs,
            "probe_alpha"=>alpha,"probe_reference_alpha"=>reference_alpha,"tau"=>times,"reference_tau"=>reference_times,
            "rules"=>[Dict("nodes"=>r.nodes,"weights"=>r.weights) for r in rules],"batch_index"=>index[],
            "sampled_alpha_error"=>scaled_error(native.values,alpha_refs,policy["relative_floor"]),
            "anchor_error"=>scaled_error(alpha[1:1],[a],policy["relative_floor"]),
            "time_reference_error"=>scaled_error(times,reference_times,policy["relative_floor"]))
    end
    models=report["native_models"]
    for spec in contract["models"][2:end]
        id=spec["id"]
        a=id=="native_literal_alpha0_0p327" ? settings["alpha0_literal"] :
            id=="native_matched_asymptotic_alpha0" ? baseline["alpha0"] : error("unsupported native diagnostic model")
        model=Dict{String,Any}("id"=>id,"role"=>spec["role"],"alpha0"=>a,"dt_trials"=>Dict{String,Any}[])
        push!(models,model)
        for (dti,dt) in enumerate(policy["native_dt_refinements"])
            level=Dict{String,Any}("dt"=>dt,"quadrature_trials"=>Dict{String,Any}[])
            push!(model["dt_trials"],level)
            for (qi,n) in enumerate(policy["quadrature_orders"])
                row=trial(a,q0,probes,nf,2,dt,n)
                push!(level["quadrature_trials"],row)
                if qi>1
                    row["quadrature_refinement_error"]=scaled_error(row["tau"],level["quadrature_trials"][end-1]["tau"],policy["relative_floor"])
                end
                if qi>=3
                    # Fix alpha's own native discretization first; higher-order quadrature cannot repair it.
                    alpha_ok=row["sampled_alpha_error"]<=policy["alpha_max_scaled"] && row["anchor_error"]<=policy["anchor_max_scaled"]
                    quadrature_ok=row["quadrature_refinement_error"]<=policy["quadrature_refinement_max_scaled"] && row["time_reference_error"]<=policy["time_max_scaled"]
                    (!alpha_ok || quadrature_ok) && break
                end
            end
            final=level["quadrature_trials"][end]
            level["selected_order"]=final["quadrature_order"]
            if dti>1
                previous=model["dt_trials"][end-1]["quadrature_trials"][end]
                level["dt_refinement_error"]=scaled_error(final["probe_alpha"],previous["probe_alpha"],policy["relative_floor"])
            end
            accepted=dti>=3 && final["sampled_alpha_error"]<=policy["alpha_max_scaled"] &&
                final["anchor_error"]<=policy["anchor_max_scaled"] && level["dt_refinement_error"]<=policy["dt_refinement_max_scaled"] &&
                final["quadrature_refinement_error"]<=policy["quadrature_refinement_max_scaled"] && final["time_reference_error"]<=policy["time_max_scaled"]
            level["accepted"]=accepted
            level["status"]=accepted ? "SELECTED_NUMERICAL_RESOLUTION" : "UNACCEPTED_REFINEMENT_TRIAL"
            println(JSON.json(Dict("model"=>id,"dt"=>dt,"order"=>final["quadrature_order"],
                "alpha_error"=>final["sampled_alpha_error"],"anchor_error"=>final["anchor_error"],
                "time_error"=>final["time_reference_error"],"accepted"=>accepted)));flush(stdout)
            accepted && break
        end
        selected=model["dt_trials"][end];final=selected["quadrature_trials"][end]
        previous=model["dt_trials"][end-1]["quadrature_trials"][end]
        check("$id native alpha/reference",final["native_alpha"],final["reference_alpha"],policy["alpha_max_scaled"])
        check("$id initial anchor",final["probe_alpha"][1:1],[a],policy["anchor_max_scaled"])
        check("$id final dt refinement",final["probe_alpha"],previous["probe_alpha"],policy["dt_refinement_max_scaled"])
        check("$id native time/reference",final["tau"],final["reference_tau"],policy["time_max_scaled"])
        check("$id final quadrature refinement",final["tau"],selected["quadrature_trials"][end-1]["tau"],policy["quadrature_refinement_max_scaled"])
        moments=Float64.(policy["Mellin_N_controls"])
        actual=[exp(transversity_moment(Int(N))*t) for t in final["tau"],N in moments]
        PN=[-2/3,-2.,-26/9]
        expected=[exp(P*t) for t in final["reference_tau"],P in PN]
        check("$id analytic Mellin factors",actual,expected,policy["mellin_factor_max_scaled"])
        # Reinitialize at Q2=10 with its native value and check the same fixed-nf flow.
        chain=trial(final["probe_alpha"][2],probes[2],probes[2:3],nf,2,final["dt"],final["quadrature_order"])
        model["composition_control"]=chain
        check("$id native alpha composition",chain["probe_alpha"][2:2],final["probe_alpha"][3:3],policy["composition_max_scaled"])
        check("$id native time composition",[final["tau"][2]+chain["tau"][2]],final["tau"][3:3],policy["composition_max_scaled"])
        model["selected_dt"]=selected["dt"];model["selected_order"]=selected["selected_order"]
        model["alpha"]=final["probe_alpha"];model["tau"]=final["tau"]
        model["difference_from_asymptotic_alpha"]=final["probe_alpha"]-baseline["alpha"]
        model["difference_from_asymptotic_tau"]=final["tau"]-baseline["tau"]
        model["difference_from_printed_alpha_MZ"]=last(final["probe_alpha"])-settings["reported_alpha_MZ"]
    end
    for flavors in policy["one_loop_nf_controls"]
        row=trial(settings["alpha0_literal"],q0,probes,flavors,1,policy["one_loop_control_dt"],64)
        push!(report["one_loop_controls"],row)
        check("one-loop nf$flavors alpha",row["native_alpha"],row["reference_alpha"],policy["one_loop_max_scaled"])
        check("one-loop nf$flavors time",row["tau"],row["reference_tau"],policy["one_loop_max_scaled"])
    end
    check("independent implicit-root residuals",[maximum(v.implicit_residual for v in values(cache))],[0.],policy["implicit_residual_max_abs"];absolute=true)
    good=native_request(.327,2.4,[2.4,10.],5,2,.05)
    bad=[("wrong protocol",replace(good,"CF-COUPLING-1"=>"OTHER")),
         ("zero query count",replace(good,"2 5 2"=>"0 5 2")),
         ("unsupported nf",replace(good,"2 5 2"=>"2 7 2")),
         ("unsupported loops",replace(good,"2 5 2"=>"2 5 3")),
         ("zero dt",replace(good,"0.327 2.4 0.05"=>"0.327 2.4 0.0")),
         ("nonfinite alpha",replace(good,"0.327 2.4 0.05"=>"NaN 2.4 0.05")),
         ("backward query",replace(good,"10.0\n"=>"1.0\n")),
         ("truncated queries",first(split(good,"10.0"))),
         ("unexpected trailer",good*"TRAILER\n")]
    for (label,request) in bad
        index[]+=1;raw=run_native_raw(exe,request,directory,index[]);push!(report["native_batches"],raw.evidence)
        rejected=raw.evidence["exit_code"]!=0 && isempty(raw.text) && occursin("REFERENCE_INPUT_OR_RESULT_REJECTED",raw.stderr)
        push!(report["protocol_controls"],Dict("name"=>label,"batch_index"=>index[],"rejected"=>rejected))
        check("native invalid request: $label",[rejected ? 0. : 1.],[0.],0.;absolute=true)
    end
    check("zero native precision diagnostics",[sum(v["precision_errors"] for v in report["native_batches"])],[0.],0.;absolute=true)
    check("source/build stability",[hashes==source_hashes() && bytes2hex(sha256(read(buildfile)))==buildhash &&
        CollinsFit.verified_hoppet_build(buildfile)==build && bytes2hex(sha256(read(paperfile)))==contract["primary_source"]["sha256"] ? 0. : 1.],[0.],0.;absolute=true)
    report["independent_reference_cache_entries"]=length(cache)
    report["status"]=!isempty(checks) && all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    report
end

function main()
    length(ARGS)==2 || error("usage: check_native_coupling.jl BUILD_JSON OUTPUT_JSON")
    report=Dict{String,Any}("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
        "physical_extraction"=>false,"G2PERF_pass"=>false,"historical_coupling_recovered"=>false,
        "physical_partial_moments_evaluated"=>false,"failed_paper_benchmark_replaced"=>false,
        "checks"=>Dict{String,Any}[],"native_models"=>Dict{String,Any}[],"native_batches"=>Dict{String,Any}[],
        "one_loop_controls"=>Dict{String,Any}[],"protocol_controls"=>Dict{String,Any}[],
        "component_scope"=>"Native alpha/time comparator and analytic Mellin controls only; not historical prescription admission",
        "julia_version"=>string(VERSION))
    try
        coupling_checks!(report,ARGS[1],ARGS[2])
    catch exception
        report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(ARGS[2],report)
    end
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "error"=>get(report,"error",nothing),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
end

if abspath(PROGRAM_FILE)==(@__FILE__)
    main()
end
