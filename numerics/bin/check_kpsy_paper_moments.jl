# A scoped publication benchmark. No optimizer, production gate or synthetic substitute.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra
include("PaperMoments.jl")
using .PaperMoments

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/kpsy-paper-moment-contract.json")
const SOURCE_FILES=vcat([@__FILE__,joinpath(@__DIR__,"PaperMoments.jl"),CONTRACT,
    joinpath(ROOT,"scripts/prepare_kpsy_paper.py"),joinpath(ROOT,"config/initial-contract.json"),
    joinpath(ROOT,"config/evolution-contract.json"),joinpath(ROOT,"config/parton-contract.json"),
    joinpath(ROOT,"config/parton-sources.json"),joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"numerics/Project.toml"),
    joinpath(ROOT,"numerics/Manifest.toml")],filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in SOURCE_FILES)

function check_moments!(report,nativefile,partonfile,paperfile)
    contract=JSON.parsefile(CONTRACT);policy=contract["numerical_policy"];definition=contract["definition"]
    sources=source_hashes();report["sources"]=sources
    report["contract_id"]=contract["id"];report["contract_sha256"]=bytes2hex(sha256(read(CONTRACT)))
    report["source_paper"]=contract["primary_source"]
    report["matrix_encoding"]="Julia matrices are serialized as outer arrays of columns"
    bytes2hex(sha256(read(paperfile)))==contract["primary_source"]["sha256"] &&
        filesize(paperfile)==contract["primary_source"]["bytes"] || error("primary paper bytes differ from the frozen contract")
    native=CollinsFit.verified_hoppet_build(nativefile);parton=CollinsFit.verify_parton_build(partonfile)
    build_hashes=Dict(p=>bytes2hex(sha256(read(p))) for p in (nativefile,partonfile))
    # Reusing a build from another checkout is safe only with identical current adapter/input contracts.
    for (build,paths) in ((native,["config/hoppet-contract.json","config/evolution-contract.json",
        "config/matched-rows-contract.json","config/hoppet-batch-contract.json","fortran/hoppet_bridge.f90",
        "src/collinsflow/hoppet.py","src/collinsflow/provenance.py","scripts/prepare_hoppet.py"]),
        (parton,["config/parton-sources.json","src/collinsflow/parton.py",
        "scripts/prepare_parton_inputs.py","bridges/lhapdf.cpp"]))
        for path in paths
            matches=filter(p->endswith(replace(p,'\\'=>'/'),"/"*path),collect(keys(build["source_hashes"])))
            length(matches)==1 && build["source_hashes"][only(matches)]==bytes2hex(sha256(read(joinpath(ROOT,path)))) ||
                error("build does not bind the current project source: $path")
        end
    end
    report["build_sha256"]=build_hashes;report["native_build"]=native;report["parton_build"]=parton
    theta=collect(kpsy_published_parameters().values)
    initial_contract=JSON.parsefile(joinpath(ROOT,"config/initial-contract.json"))
    theta==Float64.(initial_contract["published_parameters"]) && initial_contract["Q02_GeV2"]==definition["Q02_GeV2"] ||
        error("initial parameter/scale contract mismatch")
    report["printed_parameters"]=theta;report["definition"]=definition
    report["printed_value_policy"]=contract["printed_value_policy"]
    report["limitations"]=contract["limitations"]
    checks=report["checks"];native_ladder=report["native_refinements"];dense_ladder=report["dense_refinements"]
    check(name,a,b,tol;metric=:scaled)=push!(checks,compare_check(name,a,b,tol;metric=metric,floor=policy["relative_floor"]))
    lower,upper=Float64.(definition["x_range"]);q2,q0=definition["Q2_GeV2"],definition["Q02_GeV2"]
    coupling=AsymptoticTwoLoopCoupling(lambda_GeV=contract["coupling"]["Lambda_GeV"],nf=contract["coupling"]["nf"])
    pdf=CT10Input(partonfile);hel=load_dssv08(parton["input_root"]);ff=load_dss14(parton["input_root"])
    last_native=nothing;last_dense=nothing
    function native_run(resolution)
        plan=nothing
        try
            elapsed=@elapsed begin
                plan=HoppetPlan(nativefile,[q2],q0,coupling;xmin=lower,dy=-log(lower)/resolution,
                    layout=:single,order=policy["native_order"],eps=policy["native_eps"],du=policy["native_du"])
                xs,soffer=CollinsFit.prepare_initial_soffer(pdf,hel,plan.grid)
                p=KPSYParameters(theta)
                initial=[initial_transversity(p,pid,x,soffer[f,j]) for (j,x) in enumerate(xs),(f,pid) in enumerate((2,1))]
                independent=direct_initial(xs,soffer,theta;bits=policy["rounding_reference_bits"])
                evolved=apply_operator(plan,initial)[:,:,1]
                weights=moment_weights(xs,lower,upper)
                moments=vec(transpose(evolved)*weights)
                functional=vec(operator_adjoint(plan,reshape(weights,:,1,1)))
                check("native $resolution initial formula",initial,independent,policy["initial_formula_max_scaled"])
                check("native $resolution moment adjoint",vec(transpose(initial)*functional),moments,policy["adjoint_max_scaled"])
                check("native $resolution exact endpoint",evolved[end,:],[0.,0.],policy["endpoint_max_abs"];metric=:absolute)
            end
            result=Dict{String,Any}("intervals"=>resolution,"grid"=>xs,"soffer"=>soffer,
                "initial"=>initial,"independent_initial"=>independent,"evolved"=>evolved,"weights"=>weights,
                "moment_functional"=>functional,"moments"=>moments,"elapsed_seconds_informational"=>elapsed,
                "production_performance_admission"=>false)
            if !isempty(native_ladder)
                previous=native_ladder[end]
                result["refinement_error"]=scaled_error(moments,previous["moments"],policy["relative_floor"])
                result["refinement_within_policy"]=result["refinement_error"]<=policy["refinement_max_scaled"]
            end
            push!(native_ladder,result)
            println(JSON.json(Dict("stage"=>"native","intervals"=>resolution,"moments"=>moments,
                "refinement_error"=>get(result,"refinement_error",nothing))));flush(stdout)
            result
        finally
            plan===nothing || close(plan)
        end
    end
    function dense_run(resolution)
        elapsed=@elapsed begin
            grids=initial_source_grids(pdf,hel,ff;background_nodes=resolution,x_probes=[lower,upper])
            grid=sort(unique(vcat([lower,upper,1.],filter(x->x>=lower,grids.x))))
            xs,soffer=CollinsFit.prepare_initial_soffer(pdf,hel,grid)
            initial=direct_initial(xs,soffer,theta;bits=policy["rounding_reference_bits"])
            op=collinear_lo_operator(xs,q2,q0,coupling)
            evolved=apply_operator(op,initial);weights=moment_weights(xs,lower,upper)
            moments=vec(transpose(evolved)*weights)
            reference=[reference_integral(xs,evolved[:,f],lower,upper;rtol=policy["reference_quad_rtol"],
                atol=policy["reference_quad_atol_per_cell"]) for f in 1:2]
            check("dense $resolution independent log-x integration",moments,[v.value for v in reference],policy["integration_max_scaled"])
            check("dense $resolution quadrature error estimate",[v.estimated_error/max(abs(v.value),policy["relative_floor"]) for v in reference],
                [0.,0.],policy["integration_max_scaled"];metric=:absolute)
            check("dense $resolution exact endpoint",evolved[end,:],[0.,0.],policy["endpoint_max_abs"];metric=:absolute)
        end
        result=Dict{String,Any}("background_nodes"=>resolution,"grid"=>xs,"soffer"=>soffer,"initial"=>initial,
            "evolved"=>evolved,"weights"=>weights,"moments"=>moments,
            "independent_log_quadrature"=>[Dict("value"=>v.value,"estimated_error"=>v.estimated_error) for v in reference],
            "elapsed_seconds_informational"=>elapsed,"production_performance_admission"=>false)
        if !isempty(dense_ladder)
            previous=dense_ladder[end]
            result["refinement_error"]=scaled_error(moments,previous["moments"],policy["relative_floor"])
            result["refinement_within_policy"]=result["refinement_error"]<=policy["refinement_max_scaled"]
        end
        push!(dense_ladder,result)
        println(JSON.json(Dict("stage"=>"dense","background_nodes"=>resolution,"moments"=>moments,
            "refinement_error"=>get(result,"refinement_error",nothing))));flush(stdout)
        GC.gc();result
    end
    try
        for n in policy["native_intervals"][1:3];last_native=native_run(n);end
        for n in policy["dense_background_nodes"][1:3];last_dense=dense_run(n);end
        ni=di=3
        while true
            cross=scaled_error(last_native["moments"],last_dense["moments"],policy["relative_floor"])
            extend_native=(!last_native["refinement_within_policy"] || cross>policy["cross_max_scaled"]) && ni<length(policy["native_intervals"])
            extend_dense=(!last_dense["refinement_within_policy"] || cross>policy["cross_max_scaled"]) && di<length(policy["dense_background_nodes"])
            extend_native || extend_dense || break
            if extend_native;ni+=1;last_native=native_run(policy["native_intervals"][ni]);end
            if extend_dense;di+=1;last_dense=dense_run(policy["dense_background_nodes"][di]);end
        end
        check("native final refinement",last_native["moments"],native_ladder[end-1]["moments"],policy["refinement_max_scaled"])
        check("dense final refinement",last_dense["moments"],dense_ladder[end-1]["moments"],policy["refinement_max_scaled"])
        check("independent evolution cross comparison",last_native["moments"],last_dense["moments"],policy["cross_max_scaled"])
        numerical_end=length(checks)
        centers=[contract["reported_centers"][f] for f in ("u","d")]
        for f in 1:2
            check("paper printed raw center $(f==1 ? "u" : "d")",[last_native["moments"][f]],[centers[f]],
                contract["printed_value_policy"]["raw_center_max_abs"];metric=:absolute)
        end
        rounding=contract["printed_value_policy"]["parameter_rounding"]
        vertex_parameters=rounding_vertices(theta,rounding["bounds"])
        vertex_moments=[vec(transpose(direct_initial(last_native["grid"],last_native["soffer"],t;
            bits=policy["rounding_reference_bits"]))*last_native["moment_functional"]) for t in vertex_parameters]
        report["parameter_rounding_diagnostic"]=Dict("status"=>"DIAGNOSTIC_ONLY_NOT_AN_ACCEPTANCE_GATE",
            "certified_continuous_extrema"=>false,"bounds"=>rounding["bounds"],"center_moments"=>last_native["moments"],
            "vertex_parameters"=>vertex_parameters,"vertex_moments"=>vertex_moments,
            "sampled_min"=>[minimum(v[f] for v in vertex_moments) for f in 1:2],
            "sampled_max"=>[maximum(v[f] for v in vertex_moments) for f in 1:2])
        check("source/build hashes stable",[source_hashes()==sources && all(bytes2hex(sha256(read(p)))==d for (p,d) in build_hashes) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(partonfile)==parton &&
            bytes2hex(sha256(read(paperfile)))==contract["primary_source"]["sha256"] ? 0. : 1.],[0.],0.;metric=:absolute)
        report["numerical_checks_passed"]=all(c["passed"] for c in checks[1:numerical_end])
        report["raw_paper_centers_passed"]=all(c["passed"] for c in checks[numerical_end+1:numerical_end+2])
        report["selected_native_intervals"]=last_native["intervals"]
        report["selected_dense_background_nodes"]=last_dense["background_nodes"]
    finally
        close(pdf)
    end
    report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    report
end

function main()
    length(ARGS)==4 || error("usage: check_kpsy_paper_moments.jl NATIVE_BUILD PARTON_BUILD LOCKED_PAPER_HTML OUTPUT_JSON")
    report=Dict{String,Any}("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "G2PERF_pass"=>false,"status"=>"FAIL","checks"=>Dict{String,Any}[],
        "native_refinements"=>Dict{String,Any}[],"dense_refinements"=>Dict{String,Any}[],"julia_version"=>string(VERSION))
    try
        check_moments!(report,ARGS[1],ARGS[2],ARGS[3])
    catch exception
        report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(ARGS[4],report)
    end
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "error"=>get(report,"error",nothing),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
end

if abspath(PROGRAM_FILE)==(@__FILE__)
    main()
end
