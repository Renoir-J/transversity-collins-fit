module BoundedJointAudit
using CollinsFit, JSON, SHA, LinearAlgebra
const PROJECT=normpath(joinpath(@__DIR__,"../.."))
require(value,message)=value || error(message)
function close(a,b;tol=1e-11,floor=1e-8)
    x=vec(Float64.(a));y=vec(Float64.(b))
    require(length(x)==length(y) && all(isfinite,x) && all(isfinite,y),"invalid operands")
    require(maximum(abs.(x-y)./max.(abs.(y),floor))<=tol,"changed numerical operand")
end
"Exact replay of pinned NLopt2.11 BOBYQA initial rescaling/bound preparation, not a fit-coordinate transform."
function expected_first_evaluation(algorithm,start,steps,gh)
    algorithm==:LN_SBPLX && return Float64.(start)
    require(algorithm==:LN_BOBYQA && length(start)==length(steps)==13,"invalid initial mapping")
    scale=ones(13)
    if any(steps[i]!=steps[i-1] for i in 2:13)
        for i in 2:13;scale[i]=steps[i]/steps[1];end
    end
    require(all(s->isfinite(s) && s>0,scale),"unsupported native scale range")
    bounds=kpsy_parameter_bounds(gh);lo=collect(bounds.lower)./scale;hi=collect(bounds.upper)./scale
    x=Float64.(start)./scale;rho=abs(steps[1]/scale[1])
    for i in eachindex(x)
        sl=lo[i]-x[i];su=hi[i]-x[i]
        if sl>=-rho
            x[i]=sl>=0. ? lo[i] : lo[i]+rho
        elseif su<=rho
            x[i]=su<=0. ? hi[i] : hi[i]-rho
        end
    end
    [min(max(lo[i],x[i]+0.),hi[i])*scale[i] for i in eachindex(x)]
end
const EXPECTED=vcat(["numerics/bin/check_bounded_joint_fit.jl","numerics/bin/BoundedJointReference.jl",
    "numerics/bin/JointLikelihoodReference.jl","numerics/Project.toml","numerics/Manifest.toml"],
    ["numerics/src/"*name for name in readdir(joinpath(PROJECT,"numerics/src")) if endswith(name,".jl")],
    ["config/"*name for name in readdir(joinpath(PROJECT,"config")) if endswith(name,"-contract.json")])
function sources_current(report)
    require(Set(keys(report["sources"]))==Set(EXPECTED),"incomplete producer inventory")
    for path in EXPECTED
        require(report["sources"][path]==bytes2hex(sha256(read(joinpath(PROJECT,path)))),"stale source $path")
    end
end
include("BoundedJointReference.jl")
using .BoundedJointReference
const AUDIT_SOURCES=("numerics/bin/audit_bounded_joint_fit.jl","numerics/bin/BoundedJointAudit.jl",
    "numerics/bin/check_bounded_joint_cache.jl","numerics/test/bounded_report_controls.jl",
    "workflow/bounded-joint-fit.smk","scripts/check_cached_bounded_workflow.py")
audit_source_hashes()=Dict(path=>bytes2hex(sha256(read(joinpath(PROJECT,path)))) for path in AUDIT_SOURCES)

function audit_bounded_report(report)
    sources_current(report)
    initial_policy=JSON.parsefile(joinpath(PROJECT,"config/bounded-initial-roundoff-contract.json"))
    require(initial_policy["id"]=="KPSY-BOUNDED-INITIAL-IEEE-001" && initial_policy["scientific_acceptance_thresholds_changed"]===false,"changed initial arithmetic policy")
    contract=JSON.parsefile(joinpath(PROJECT,"config/bounded-joint-fit-contract.json"))
    require(report["contract"]==contract && report["contract_id"]==contract["id"],"changed contract")
    require(report["contract_sha256"]==bytes2hex(sha256(read(joinpath(PROJECT,"config/bounded-joint-fit-contract.json")))),"stale contract")
    require(!haskey(report,"error") && report["status"] in ("PASS","FAIL","INCOMPLETE"),"incomplete diagnostic report")
    require(report["scope"]=="synthetic_closure" && report["theory_scope"]=="reconstructed_theory","changed scope")
    for key in ("production_ready","physical_extraction","G2PERF_pass","experimental_data_admitted","real_data_fit_implemented","interval_coverage_verified")
        require(report[key]===false,"invalid admission flag $key")
    end
    budget=report["maxeval_per_start"]
    require(budget in contract["settings"]["resource_ladder_maxeval"],"changed resource policy")
    algorithm=Symbol(get(report,"algorithm","LN_SBPLX"))
    require(algorithm in (:LN_SBPLX,:LN_BOBYQA),"unknown solver")
    if haskey(report,"solver_comparison_contract")
        path=joinpath(PROJECT,"config/bounded-solver-comparison-contract.json")
        require(report["solver_comparison_contract"]==JSON.parsefile(path) &&
            report["solver_comparison_contract_sha256"]==bytes2hex(sha256(read(path))),"changed prospective solver comparison")
    else
        require(algorithm==:LN_SBPLX,"missing solver policy")
    end
    require(report["matrix_encoding"]=="outer arrays are Julia matrix columns","unknown matrix encoding")
    require(report["julia_version"]==string(VERSION),"different Julia runtime requires separate archival review")
    fixture=closure_fixture();l=fixture.likelihood;p=contract["validation"];fit=report["fit"]
    expected_solver=algorithm==:LN_SBPLX ? contract["id"] : "KPSY-BOUNDED-SOLVER-COMPARISON-001"
    require(fit["contract_id"]==contract["id"] && fit["solver_contract_id"]==expected_solver,"changed fitting contract")
    require(report["truth"]==fixture.truth && report["starts_requested"]==fixture.starts,"changed truth/starts")
    require(report["likelihood"]==JSON.parse(JSON.json(joint_likelihood_metadata(l)))==fit["likelihood"],"changed prepared model metadata")
    data=report["data"]
    require(data["point_ids"]==collect(l.point_ids) && data["nuisance_ids"]==collect(l.nuisance_ids),"changed data order")
    for (name,value) in (("observations",l.observations),("covariance",l.covariance),("additive",l.additive),("relative",l.relative),("nuisance_prior",l.nuisance_prior))
        require(data[name]==JSON.parse(JSON.json(value)),"changed synthetic data $name")
    end
    settings=BoundedJointSettings(algorithm=algorithm,maxeval=budget,trace=true,source="frozen synthetic 52-observation closure")
    require(fit["settings"]==JSON.parse(JSON.json(CollinsFit.bounded_settings_metadata(settings,l.model.gh))),"changed solver settings")
    require(fit["runtime"]==CollinsFit.bounded_runtime_metadata(),"stale solver wrapper/library")
    require(fit["sources_and_runtime_stable"]===true && fit["parameter_count"]==13 && fit["point_count"]==52,"unstable or wrong fitting dimensions")
    require(fit["scope"]=="synthetic_closure" && fit["theory_scope"]=="reconstructed_theory","changed child scope")
    for key in ("production_ready","physical_extraction","G2PERF_pass","experimental_data_admitted","real_data_fit_implemented","global_optimum_certified","interval_coverage_verified")
        require(fit[key]===false,"invalid child admission $key")
    end
    candidates=fit["starts"];require(length(candidates)==length(fixture.starts)==2,"missing requested start")
    trace_counts=Int[];statuses=String[]
    for (i,row) in enumerate(candidates)
        require(row["start"]==fixture.starts[i],"changed initial point")
        name=lowercase(row["algorithm_name"])
        identity=algorithm==:LN_SBPLX ? occursin("subplex",name) || occursin("sbplx",name) : occursin("bobyqa",name)
        require(identity,"changed native algorithm")
        require(row["postcheck_passed"]===true && row["postcheck_evaluations"]==1,"missing postcheck")
        require(row["evaluation_accounting_matches"]===true,"false native accounting")
        trace=row["trace"];n=length(trace)
        require(n>0 && n==row["callback_attempts"]==row["successful_evaluations"]==row["solver_evaluations"],"incomplete evaluation trace")
        expected_first=expected_first_evaluation(algorithm,fixture.starts[i],settings.initial_step,l.model.gh)
        require(reinterpret(UInt64,Float64.(trace[1]["theta"]))==reinterpret(UInt64,expected_first),"changed native initial evaluation")
        for (j,event) in enumerate(trace)
            require(event["evaluation"]==j,"changed trace order")
            CollinsFit.bounded_start(event["theta"],l.model.gh)
            require(event["chi2"] isa Real && isfinite(event["chi2"]) && event["chi2"]>=0,"nonfinite trace objective")
        end
        best=trace[argmin([e["chi2"] for e in trace])]
        require(row["retained_best_evaluation"]==Dict("theta"=>best["theta"],"chi2"=>best["chi2"]),"changed retained best")
        require(row["theta"]==best["theta"],"returned point absent from best trace")
        actual=joint_objective(l,Float64.(row["theta"]))
        close([row["chi2"],row["reported_chi2"],best["chi2"]],fill(actual,3))
        close([row["postcheck_difference"]],[abs(row["reported_chi2"]-row["chi2"])];tol=0.)
        stop=row["solver_stop"]
        status=stop in ("FTOL_REACHED","XTOL_REACHED") ? "PASS" : stop in ("MAXEVAL_REACHED","MAXTIME_REACHED") ? "INCOMPLETE" : "FAIL"
        require(row["status"]==status && row["optimizer_tolerance_stop"]==(status=="PASS"),"false convergence")
        stop=="MAXEVAL_REACHED" && require(n>=budget,"false budget exhaustion")
        push!(trace_counts,n);push!(statuses,status)
    end
    fit_status=any(==("FAIL"),statuses) ? "FAIL" : all(==("PASS"),statuses) ? "PASS" : "INCOMPLETE"
    require(fit["status"]==fit_status,"changed whole-fit status")
    best=argmin([r["chi2"] for r in candidates])
    require(fit["best_start_index"]==best && fit["best_candidate_status"]==statuses[best] && fit["theta"]==candidates[best]["theta"],"changed selected start")
    profile=joint_profile(l,Float64.(fit["theta"]))
    close([fit["chi2"]],[profile.chi2]);close(fit["predictions"],profile.predictions);close(fit["nuisance"],profile.nuisance)
    require(fit["report_profile_attempts"]==1,"missing final reporting profile")
    qr=setprecision(256) do
        augmented_reference(BigFloat.(profile.predictions),l.observations,l.covariance,l.additive,l.relative,l.nuisance_prior)
    end
    qr_error=abs(profile.chi2-Float64(qr.chi2))/max(abs(Float64(qr.chi2)),1e-8)
    flipped=copy(fixture.truth);flipped[[1,2,7,8]].*=-1
    distances=[maximum(abs.(fit["theta"]-truth)./Float64.(p["parameter_scales"])) for truth in (fixture.truth,flipped)]
    prediction_error=maximum(abs.(profile.predictions-l.observations));closure=report["closure"]
    close([closure["chi2"],closure["independent_QR_chi2"]],[profile.chi2,Float64(qr.chi2)])
    close([closure["max_absolute_prediction_error"]],[prediction_error]);close(closure["scaled_parameter_distances_to_both_sign_branches"],distances)
    require(closure["selected_sign_branch"]==argmin(distances) && closure["candidate_predictions_recomputed"]===true && closure["global_optimality_or_coverage_verified"]===false,"changed closure claim")
    information=closure_information(fixture);stored=report["local_information"]
    close(stored["singular_values"],information["singular_values"])
    for key in ("rank","relative_rank_threshold","full_local_rank","global_identifiability_verified","interval_coverage_verified","production_ready")
        require(stored[key]==information[key],"changed local information $key")
    end
    boundary=copy(fixture.truth);boundary[1]=1.;boundary[2]=-1.;boundary[3]=0.;boundary[5]=0.
    control=report["boundary_resource_control"];br=only(control["starts"])
    require(control["scope"]=="synthetic_closure" && control["theory_scope"]=="reconstructed_theory","changed boundary scope")
    for key in ("production_ready","physical_extraction","G2PERF_pass","experimental_data_admitted","real_data_fit_implemented","global_optimum_certified","interval_coverage_verified")
        require(control[key]===false,"invalid boundary admission $key")
    end
    boundary_settings=BoundedJointSettings(algorithm=algorithm,maxeval=1,source="deliberate physical-family boundary resource stop")
    require(control["settings"]==JSON.parse(JSON.json(CollinsFit.bounded_settings_metadata(boundary_settings,l.model.gh))) &&
        control["runtime"]==fit["runtime"] && control["sources_and_runtime_stable"]===true,"changed boundary settings/runtime")
    require(control["solver_contract_id"]==expected_solver && control["theta"]==boundary &&
        control["best_start_index"]==1 && control["report_profile_attempts"]==1,"changed boundary candidate")
    close(control["predictions"],predict_joint(l.model,boundary));close(control["nuisance"],zeros(2))
    require(control["status"]=="INCOMPLETE" && br["status"]=="INCOMPLETE" && br["solver_stop"]=="MAXEVAL_REACHED" && br["optimizer_tolerance_stop"]===false,"boundary resource stop promoted")
    require(br["start"]==br["theta"]==boundary && br["solver_evaluations"]==br["callback_attempts"]==br["successful_evaluations"]==1,"changed boundary accounting")
    require(length(br["trace"])==1 && br["trace"][1]["theta"]==boundary && br["trace"][1]["chi2"]==0. && control["chi2"]==0.,"changed boundary objective")
    require(all(isfinite,predict_joint(l.model,boundary)),"nonfinite process boundary")
    expected=Dict("52 explicit synthetic observations"=>length(l.point_ids)==52,
        "full local thirteen-coordinate information"=>information["full_local_rank"],
        "every requested start retained"=>length(candidates)==length(fixture.starts),
        "every start reaches configured tolerance stop"=>fit_status=="PASS",
        "reported sources and runtime stable"=>true,
        "independent augmented QR objective"=>qr_error<=1e-11,
        "noiseless closure chi2"=>profile.chi2<=p["closure_chi2_max"],
        "noiseless prediction recovery"=>prediction_error<=p["closure_max_absolute_prediction_error"],
        "all thirteen parameters recovered up to global sign"=>minimum(distances)<=p["closure_max_scaled_parameter_error"],
        "closed physical-family boundary stays finite and INCOMPLETE"=>true,
        "source files stable during complete job"=>true)
    checks=report["checks"]
    require(length(checks)==length(expected) && Set(row["name"] for row in checks)==Set(keys(expected)),"incomplete check ledger")
    require(all(row["passed"]===expected[row["name"]] for row in checks),"false check bit")
    expected_status=fit_status=="FAIL" || qr_error>1e-11 ? "FAIL" : fit_status=="INCOMPLETE" ? "INCOMPLETE" : all(values(expected)) ? "PASS" : "FAIL"
    require(report["status"]==expected_status,"false producer status")
    sources_current(report)
    Dict("status"=>report["status"],"review_status"=>"RETAINED_REPORT_ARITHMETIC_VERIFIED","producer_status"=>report["status"],
        "accepted_fitting_gate"=>report["status"]=="PASS","scope"=>"synthetic_closure_retained_report_review",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,
        "algorithm"=>String(algorithm),"source_bindings"=>length(EXPECTED),"trace_counts"=>trace_counts,"all_trace_parameter_domains_verified"=>true,
        "all_trace_objectives_replayed"=>false,"candidate_objectives_replayed"=>2,
        "selected_candidate_predictions_replayed"=>true,"independent_augmented_QR_scaled_error"=>qr_error,
        "independent_local_rank"=>information["rank"],"recomputed_checks"=>expected,
        "max_scaled_parameter_error_up_to_sign"=>minimum(distances),"max_absolute_prediction_error"=>prediction_error,
        "global_optimum_certified"=>false,"interval_coverage_verified"=>false)
end

end
