using CollinsFit, JSON, SHA, Dates, LinearAlgebra, Statistics
include("BoundedJointReference.jl")
using .BoundedJointReference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/bounded-joint-fit-contract.json")
const SOURCES=vcat([@__FILE__,joinpath(@__DIR__,"BoundedJointReference.jl"),joinpath(@__DIR__,"JointLikelihoodReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)),
    filter(p->endswith(p,"-contract.json"),readdir(joinpath(ROOT,"config");join=true)))
source_hashes()=Dict(replace(relpath(path,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(path))) for path in SOURCES)

function run_checks!(report,budget)
    contract=JSON.parsefile(CONTRACT);p=contract["validation"]
    budget in contract["settings"]["resource_ladder_maxeval"] || error("use an explicit frozen resource-ladder budget")
    before=source_hashes();report["sources"]=before;report["contract"]=contract
    report["contract_sha256"]=bytes2hex(sha256(read(CONTRACT)))
    println(JSON.json(Dict("stage"=>"preparing 52 synthetic observations","status"=>"IN_PROGRESS_NOT_A_GATE")));flush(stdout)
    fixture=closure_fixture();l=fixture.likelihood
    report["truth"]=fixture.truth;report["starts_requested"]=fixture.starts
    report["likelihood"]=joint_likelihood_metadata(l)
    report["data"]=Dict("point_ids"=>collect(l.point_ids),"observations"=>l.observations,"covariance"=>l.covariance,
        "additive"=>l.additive,"relative"=>l.relative,"nuisance_prior"=>l.nuisance_prior,"nuisance_ids"=>collect(l.nuisance_ids))
    information=closure_information(fixture);report["local_information"]=information
    checks=report["checks"]
    check(name,passed;details=nothing)=push!(checks,Dict("name"=>name,"passed"=>Bool(passed),"details"=>details))
    check("52 explicit synthetic observations",length(l.point_ids)==52)
    check("full local thirteen-coordinate information",information["full_local_rank"];details=information)
    information["full_local_rank"] || error("fixture does not certify identifiable thirteen-parameter recovery")
    joint_objective(l,fixture.truth)
    times=[@elapsed(joint_objective(l,fixture.truth)) for _ in 1:10]
    report["objective_timing_informational"]=Dict("samples"=>10,"median_seconds"=>median(times),
        "maximum_seconds"=>maximum(times),"production_performance_accepted"=>false)
    println(JSON.json(Dict("stage"=>"bounded multistart fit","observations"=>52,"parameters"=>13,
        "maxeval_per_start"=>budget,"objective_median_seconds_informational"=>median(times),"status"=>"IN_PROGRESS_NOT_A_GATE")));flush(stdout)
    settings=BoundedJointSettings(maxeval=budget,trace=true,source="frozen synthetic 52-observation closure")
    result=fit_joint(l,fixture.starts;settings=settings);report["fit"]=result
    check("every requested start retained",length(result["starts"])==length(fixture.starts))
    check("every start reaches configured tolerance stop",result["status"]=="PASS";details=[r["solver_stop"] for r in result["starts"]])
    check("reported sources and runtime stable",result["sources_and_runtime_stable"])
    hard_failure=result["status"]=="FAIL"
    if result["theta"]!==nothing && result["chi2"]!==nothing
        theta=Float64.(result["theta"]);profile=joint_profile(l,theta)
        independent=setprecision(256) do
            augmented_reference(BigFloat.(profile.predictions),l.observations,l.covariance,l.additive,l.relative,l.nuisance_prior)
        end
        scaled_difference=abs(profile.chi2-Float64(independent.chi2))/max(abs(Float64(independent.chi2)),1e-8)
        check("independent augmented QR objective",scaled_difference<=1e-11;
            details=Dict("actual"=>profile.chi2,"reference"=>Float64(independent.chi2),"scaled_error"=>scaled_difference,"limit"=>1e-11))
        hard_failure |= scaled_difference>1e-11
        flipped=copy(fixture.truth);flipped[[1,2,7,8]].*=-1
        scales=Float64.(p["parameter_scales"])
        distances=[maximum(abs.(theta-truth)./scales) for truth in (fixture.truth,flipped)]
        prediction_error=maximum(abs.(profile.predictions-l.observations))
        report["closure"]=Dict("chi2"=>profile.chi2,"max_absolute_prediction_error"=>prediction_error,
            "scaled_parameter_distances_to_both_sign_branches"=>distances,"selected_sign_branch"=>argmin(distances),
            "independent_QR_chi2"=>Float64(independent.chi2),"candidate_predictions_recomputed"=>true,
            "global_optimality_or_coverage_verified"=>false)
        check("noiseless closure chi2",profile.chi2<=p["closure_chi2_max"];details=p["closure_chi2_max"])
        check("noiseless prediction recovery",prediction_error<=p["closure_max_absolute_prediction_error"];details=p["closure_max_absolute_prediction_error"])
        check("all thirteen parameters recovered up to global sign",minimum(distances)<=p["closure_max_scaled_parameter_error"];details=p["closure_max_scaled_parameter_error"])
    else
        check("finite candidate required for closure",false)
    end
    # Evaluate actual process rows at unsupported-AD boundaries without disguising a short run as convergence.
    boundary=copy(fixture.truth);boundary[1]=1.;boundary[2]=-1.;boundary[3]=0.;boundary[5]=0.
    by=predict_joint(l.model,boundary)
    bl=JointGaussianLikelihood(l.model,l.point_ids,by,l.covariance;additive=l.additive,relative=l.relative,
        nuisance_prior=l.nuisance_prior,nuisance_ids=l.nuisance_ids,scope="synthetic_closure",source="explicit synthetic boundary path")
    boundary_result=fit_joint(bl,[boundary];settings=BoundedJointSettings(maxeval=1,source="deliberate physical-family boundary resource stop"))
    report["boundary_resource_control"]=boundary_result
    check("closed physical-family boundary stays finite and INCOMPLETE",boundary_result["status"]=="INCOMPLETE" &&
        boundary_result["chi2"]!==nothing && isfinite(boundary_result["chi2"]) && boundary_result["starts"][1]["start"]==boundary)
    stable=source_hashes()==before;check("source files stable during complete job",stable)
    report["status"]=hard_failure || !stable ? "FAIL" : result["status"]=="INCOMPLETE" ? "INCOMPLETE" :
        all(check["passed"] for check in checks) ? "PASS" : "FAIL"
end

function main()
    length(ARGS)==2 || error("usage: check_bounded_joint_fit.jl FROZEN_MAXEVAL OUTPUT")
    budget=parse(Int,ARGS[1]);output=ARGS[2];ispath(output) && error("use a fresh fitting report path")
    report=Dict{String,Any}("status"=>"FAIL","scope"=>"synthetic_closure","theory_scope"=>"reconstructed_theory",
        "contract_id"=>"KPSY-BOUNDED-JOINT-FIT-001","production_ready"=>false,"physical_extraction"=>false,
        "G2PERF_pass"=>false,"experimental_data_admitted"=>false,"real_data_fit_implemented"=>false,
        "interval_coverage_verified"=>false,"maxeval_per_start"=>budget,"checks"=>Dict{String,Any}[],"julia_version"=>string(VERSION))
    try
        run_checks!(report,budget)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(output,report)
    end
    failed=[entry["name"] for entry in report["checks"] if !entry["passed"]]
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "error"=>get(report,"error",nothing),"failed"=>failed)))
    report["status"]=="PASS" || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
