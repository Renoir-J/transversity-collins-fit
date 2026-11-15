# Frozen coupling-only physical sensitivity, not a refit or production receipt.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff
include("CouplingMomentReference.jl")
include("PaperMoments.jl")
using .CouplingMomentReference
using .PaperMoments: moment_weights, reference_integral, direct_initial, compare_check, scaled_error

function sensitivity!(report,nativefile,partonfile,paperfile,couplingfile,analyticfile,output)
    contract=JSON.parsefile(CouplingMomentReference.CONTRACT);p=contract["numerical_policy"];d=contract["definition"]
    report["contract_id"]=contract["id"];report["contract_sha256"]=bytes2hex(sha256(read(CouplingMomentReference.CONTRACT)))
    sources=study_sources();report["sources"]=sources;report["limitations"]=contract["limitations"]
    coupling=verify_coupling_report(couplingfile,nativefile)
    analytic=verify_analytic_report(analyticfile,couplingfile,nativefile)
    profiles=moment_profiles(coupling,couplingfile)
    [v.model_id for v in profiles]==contract["models"] || error("predefined model identities changed")
    report["profiles"]=prepared_time_metadata.(profiles)
    report["dependency_sha256"]=Dict(path=>bytes2hex(sha256(read(path))) for path in
        (nativefile,partonfile,paperfile,couplingfile,analyticfile))
    papercontract=JSON.parsefile(joinpath(CouplingMomentReference.ROOT,"config/kpsy-paper-moment-contract.json"))
    bytes2hex(sha256(read(paperfile)))==papercontract["primary_source"]["sha256"] &&
        filesize(paperfile)==papercontract["primary_source"]["bytes"] || error("primary paper mismatch")
    d["x_range"]==papercontract["definition"]["x_range"] && d["Q2_GeV2"]==papercontract["definition"]["Q2_GeV2"] &&
        d["Q02_GeV2"]==papercontract["definition"]["Q02_GeV2"] &&
        d["raw_center_max_abs"]==papercontract["printed_value_policy"]["raw_center_max_abs"] || error("historical benchmark was altered")
    parton=CollinsFit.verify_parton_build(partonfile)
    for path in ("config/parton-sources.json","src/collinsflow/parton.py","scripts/prepare_parton_inputs.py","bridges/lhapdf.cpp")
        matches=filter(k->endswith(replace(k,'\\'=>'/'),"/"*path),collect(keys(parton["source_hashes"])))
        length(matches)==1 && parton["source_hashes"][only(matches)]==
            bytes2hex(sha256(read(joinpath(CouplingMomentReference.ROOT,path)))) || error("current physical input source mismatch: $path")
    end
    report["parton_build"]=parton
    theta=collect(kpsy_published_parameters().values)
    theta==JSON.parsefile(joinpath(CouplingMomentReference.ROOT,"config/initial-contract.json"))["published_parameters"] || error("literal parameter mismatch")
    report["printed_parameters"]=theta;report["definition"]=d;report["policy"]=p
    pdf=CT10Input(partonfile);hel=load_dssv08(parton["input_root"]);ff=load_dss14(parton["input_root"])
    lower,upper=Float64.(d["x_range"]);q2,q0=d["Q2_GeV2"],d["Q02_GeV2"]
    checks=report["checks"];native=report["native_refinements"];dense=report["dense_refinements"]
    check(name,a,b,tol;metric=:scaled)=push!(checks,compare_check(name,a,b,tol;metric=metric,floor=p["relative_floor"]))
    function progress(stage,resolution,row)
        record=Dict("status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>"reconstructed_theory","production_ready"=>false,
            "stage"=>stage,"resolution"=>resolution,"model_moments"=>[r["moments"] for r in row["models"]])
        writejson(output*".progress.json",record);println(JSON.json(record));flush(stdout)
    end
    function native_run(resolution)
        row=Dict{String,Any}("intervals"=>resolution,"models"=>Dict{String,Any}[])
        xs=soffer=initial=independent=weights=nothing
        for (mi,profile) in enumerate(profiles)
            plan=HoppetPlan(nativefile,[q2],q0,profile;nf=d["nf"],xmin=lower,dy=-log(lower)/resolution,
                layout=:single,order=p["native_order"],eps=p["native_eps"],du=p["native_du"])
            try
                if mi==1
                    xs,soffer=CollinsFit.prepare_initial_soffer(pdf,hel,plan.grid)
                    initial=[initial_transversity(KPSYParameters(theta),pid,x,soffer[f,j]) for (j,x) in enumerate(xs),(f,pid) in enumerate((2,1))]
                    independent=direct_initial(xs,soffer,theta;bits=p["reference_bits"])
                    weights=moment_weights(xs,lower,upper)
                    row["grid"]=xs;row["soffer"]=soffer;row["initial"]=initial;row["independent_initial"]=independent;row["weights"]=weights
                    row["shared_input_snapshot"]=CollinsFit.numerical_snapshot((xs,soffer,initial,theta))
                    check("native $resolution independent initial formula",initial,independent,p["initial_formula_max_scaled"])
                end
                xs==plan.grid && row["shared_input_snapshot"]==CollinsFit.numerical_snapshot((plan.grid,soffer,initial,theta)) || error("cross-model native inputs differ")
                evolved=apply_operator(plan,initial)[:,:,1]
                moments=vec(transpose(evolved)*weights)
                functional=vec(operator_adjoint(plan,reshape(weights,:,1,1)))
                check("$(profile.model_id) native $resolution moment adjoint",vec(transpose(initial)*functional),moments,p["adjoint_max_scaled"])
                check("$(profile.model_id) native $resolution endpoint",evolved[end,:],[0.,0.],p["endpoint_max_abs"];metric=:absolute)
                entry=Dict{String,Any}("model_id"=>profile.model_id,"evolved"=>evolved,"moments"=>moments,
                    "moment_functional"=>functional,"tau"=>plan.tau,"signature"=>plan.signature,
                    "shared_input_snapshot"=>row["shared_input_snapshot"])
                if !isempty(native)
                    entry["refinement_error"]=scaled_error(moments,native[end]["models"][mi]["moments"],p["relative_floor"])
                end
                push!(row["models"],entry)
            finally;close(plan);end
        end
        push!(native,row);progress("native",resolution,row);row
    end
    function dense_run(resolution)
        grids=initial_source_grids(pdf,hel,ff;background_nodes=resolution,x_probes=[lower,upper])
        grid=sort(unique(vcat([lower,upper,1.],filter(x->x>=lower,grids.x))))
        xs,soffer=CollinsFit.prepare_initial_soffer(pdf,hel,grid)
        initial=direct_initial(xs,soffer,theta;bits=p["reference_bits"])
        weights=moment_weights(xs,lower,upper)
        snapshot=CollinsFit.numerical_snapshot((xs,soffer,initial,theta))
        row=Dict{String,Any}("background_nodes"=>resolution,"grid"=>xs,"soffer"=>soffer,"initial"=>initial,
            "weights"=>weights,"shared_input_snapshot"=>snapshot,"models"=>Dict{String,Any}[])
        for (mi,profile) in enumerate(profiles)
            op=collinear_lo_operator(xs,q2,q0,profile)
            snapshot==CollinsFit.numerical_snapshot((op.grid,soffer,initial,theta)) || error("cross-model dense inputs differ")
            evolved=apply_operator(op,initial);moments=vec(transpose(evolved)*weights)
            independent=[reference_integral(xs,evolved[:,f],lower,upper;rtol=p["reference_quad_rtol"],atol=p["reference_quad_atol_per_cell"]) for f in 1:2]
            check("$(profile.model_id) dense $resolution independent integration",moments,[v.value for v in independent],p["integration_max_scaled"])
            check("$(profile.model_id) dense $resolution quadrature estimate",[v.estimated_error/max(abs(v.value),p["relative_floor"]) for v in independent],[0.,0.],p["integration_max_scaled"];metric=:absolute)
            check("$(profile.model_id) dense $resolution endpoint",evolved[end,:],[0.,0.],p["endpoint_max_abs"];metric=:absolute)
            entry=Dict{String,Any}("model_id"=>profile.model_id,"evolved"=>evolved,"moments"=>moments,
                "independent_integration"=>[Dict("value"=>v.value,"estimated_error"=>v.estimated_error) for v in independent],
                "tau"=>evolution_time(profile,q2,q0),"signature"=>op.prepared_signature,"shared_input_snapshot"=>snapshot)
            if !isempty(dense);entry["refinement_error"]=scaled_error(moments,dense[end]["models"][mi]["moments"],p["relative_floor"]);end
            push!(row["models"],entry)
        end
        push!(dense,row);progress("dense",resolution,row);GC.gc();row
    end
    try
        for n in p["native_intervals"][1:3];native_run(n);end
        for n in p["dense_background_nodes"][1:3];dense_run(n);end
        ni=di=3
        while true
            cross=[scaled_error(native[end]["models"][mi]["moments"],dense[end]["models"][mi]["moments"],p["relative_floor"]) for mi in 1:3]
            extend_native=(any(r["refinement_error"]>p["refinement_max_scaled"] for r in native[end]["models"]) || any(>(p["cross_max_scaled"]),cross)) && ni<length(p["native_intervals"])
            extend_dense=(any(r["refinement_error"]>p["refinement_max_scaled"] for r in dense[end]["models"]) || any(>(p["cross_max_scaled"]),cross)) && di<length(p["dense_background_nodes"])
            extend_native || extend_dense || break
            if extend_native;ni+=1;native_run(p["native_intervals"][ni]);end
            if extend_dense;di+=1;dense_run(p["dense_background_nodes"][di]);end
        end
        derivatives=Dict{String,Any}[];dp=p["parameter_derivatives"];direction=Float64.(dp["direction"])
        for (mi,profile) in enumerate(profiles)
            id=profile.model_id;nr=native[end];row=nr["models"][mi]
            check("$id native final refinement",row["moments"],native[end-1]["models"][mi]["moments"],p["refinement_max_scaled"])
            check("$id dense final refinement",dense[end]["models"][mi]["moments"],dense[end-1]["models"][mi]["moments"],p["refinement_max_scaled"])
            check("$id independent evolution comparison",row["moments"],dense[end]["models"][mi]["moments"],p["cross_max_scaled"])
            xs,soffer,functional=nr["grid"],nr["soffer"],row["moment_functional"]
            objective=t->sum(weight*initial_transversity(KPSYParameters(t),pid,x,soffer[f,j])*functional[j]
                for (j,x) in enumerate(xs),(f,pid,weight) in ((1,2,.7),(2,1,-.4)))
            gradient=ForwardDiff.gradient(objective,theta);hessian=ForwardDiff.hessian(objective,theta);hvp=hessian*direction
            reference=initial_value_gradient(xs,soffer,theta,functional;bits=p["reference_bits"])
            check("$id full parameter gradient",gradient,Float64.(reference.gradient),dp["max_scaled"])
            check("$id Collins-only derivatives vanish",vcat(gradient[7:13],vec(hessian[7:13,:])),zeros(98),0.;metric=:absolute)
            derivative=Dict{String,Any}("model_id"=>id,"gradient"=>gradient,"hessian"=>hessian,"direction"=>direction,
                "hvp"=>hvp,"jvp"=>dot(gradient,direction),"independent_value_decimal"=>string(reference.value),
                "independent_gradient_decimal"=>string.(reference.gradient),"one_sided_trials"=>Dict{String,Any}[])
            for step in dp["fd_steps"]
                KPSYParameters(theta+2step*direction)
                trial=setprecision(p["reference_bits"]) do
                    h=BigFloat(step);t=BigFloat.(theta);v=BigFloat.(direction)
                    first=initial_value_gradient(xs,soffer,t+h*v,functional;bits=p["reference_bits"])
                    second=initial_value_gradient(xs,soffer,t+2h*v,functional;bits=p["reference_bits"])
                    jvp=(-3reference.value+4first.value-second.value)/(2h)
                    hfinite=(-3reference.gradient+4first.gradient-second.gradient)/(2h)
                    Dict("step"=>step,"jvp"=>Float64(jvp),"hvp"=>Float64.(hfinite),
                        "first_value_decimal"=>string(first.value),"second_value_decimal"=>string(second.value),
                        "first_gradient_decimal"=>string.(first.gradient),"second_gradient_decimal"=>string.(second.gradient))
                end
                check("$id one-sided JVP step=$step",[dot(gradient,direction)],[trial["jvp"]],dp["max_scaled"])
                check("$id one-sided HVP step=$step",hvp,trial["hvp"],dp["max_scaled"])
                push!(derivative["one_sided_trials"],trial)
            end
            push!(derivatives,derivative)
        end
        report["parameter_derivatives"]=derivatives
        stable=study_sources()==sources && all(bytes2hex(sha256(read(path)))==digest for (path,digest) in report["dependency_sha256"]) &&
            verify_coupling_report(couplingfile,nativefile)==coupling && verify_analytic_report(analyticfile,couplingfile,nativefile)==analytic &&
            CollinsFit.verify_parton_build(partonfile)==parton
        check("source/build/dependency stability",[stable ? 0. : 1.],[0.],0.;metric=:absolute)
        report["numerical_checks_passed"]=all(c["passed"] for c in checks)
        for (mi,profile) in enumerate(profiles),f in 1:2
            check("$(profile.model_id) paper raw center $(f==1 ? "u" : "d")",[native[end]["models"][mi]["moments"][f]],
                [d["reported_centers"][f]],d["raw_center_max_abs"];metric=:absolute)
        end
        report["raw_paper_centers_passed"]=all(c["passed"] for c in checks[end-5:end])
        report["coupling_only_sensitivity"]=[Dict("model_id"=>profiles[mi].model_id,"native_moments"=>native[end]["models"][mi]["moments"],
            "difference_from_existing_asymptotic"=>native[end]["models"][mi]["moments"]-native[end]["models"][1]["moments"],
            "dense_moments"=>dense[end]["models"][mi]["moments"]) for mi in 1:3]
        report["selected_native_intervals"]=native[end]["intervals"];report["selected_dense_background_nodes"]=dense[end]["background_nodes"]
        report["physical_partial_moments_evaluated"]=true
    finally;close(pdf);end
    report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    report
end

function main()
    length(ARGS)==6 || error("usage: check_coupling_moments.jl NATIVE_BUILD PARTON_BUILD PAPER COUPLING_REPORT ANALYTIC_REPORT OUTPUT_JSON")
    ispath(ARGS[6]) && error("refuse to overwrite earlier scientific evidence; use a new run directory")
    report=Dict{String,Any}("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "G2PERF_pass"=>false,"historical_coupling_recovered"=>false,"failed_paper_benchmark_replaced"=>false,
        "physical_partial_moments_evaluated"=>false,"numerical_checks_passed"=>false,"raw_paper_centers_passed"=>false,
        "checks"=>Dict{String,Any}[],"native_refinements"=>Dict{String,Any}[],"dense_refinements"=>Dict{String,Any}[],"julia_version"=>string(VERSION))
    try;sensitivity!(report,ARGS...);catch exception;report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace());
    finally;report["completed_at"]=string(now(UTC));writejson(ARGS[6],report);end
    println(JSON.json(Dict("status"=>report["status"],"checks"=>length(report["checks"]),"numerical_checks_passed"=>report["numerical_checks_passed"],
        "raw_paper_centers_passed"=>report["raw_paper_centers_passed"],"error"=>get(report,"error",nothing),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
end

if abspath(PROGRAM_FILE)==(@__FILE__);main();end
