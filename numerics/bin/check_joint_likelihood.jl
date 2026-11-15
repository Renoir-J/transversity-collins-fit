using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff
include("JointLikelihoodReference.jl")
include("PaperMoments.jl")
using .JointLikelihoodReference
using .PaperMoments: compare_check

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/joint-likelihood-contract.json")
const INPUT_CONTRACT=joinpath(ROOT,"config/joint-pinned-input-contract.json")
const SOURCES=vcat([@__FILE__,joinpath(@__DIR__,"JointLikelihoodReference.jl"),joinpath(@__DIR__,"PaperMoments.jl"),
    CONTRACT,INPUT_CONTRACT,joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    [joinpath(ROOT,"config",n*"-contract.json") for n in ("initial","evolution","hoppet","parton","polarized-process",
        "leptonic-observable","leptonic-bin","sidis-convention","star-evolved","star-born")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in SOURCES)

function checks!(report,kind,nativefile,partonfile)
    c=JSON.parsefile(CONTRACT);ip=JSON.parsefile(INPUT_CONTRACT);p=c["validation_policy"]
    hashes=source_hashes();report["sources"]=hashes
    report["contract_id"]=c["id"];report["contract_sha256"]=bytes2hex(sha256(read(CONTRACT)))
    report["input_contract_sha256"]=bytes2hex(sha256(read(INPUT_CONTRACT)))
    report["limitations"]=vcat(c["limitations"],ip["limitations"])
    report["policy"]=p;report["matrix_encoding"]="outer arrays are Julia matrix columns"
    dependencies=Dict{String,String}()
    if kind=="inputs"
        native=CollinsFit.verified_hoppet_build(nativefile);parton=CollinsFit.verify_parton_build(partonfile)
        report["native_build"]=native;report["parton_build"]=parton
        for path in (nativefile,partonfile);dependencies[path]=bytes2hex(sha256(read(path)));end
        report["input_preparation_policy"]=ip
    end
    report["dependency_sha256"]=dependencies
    println(JSON.json(Dict("stage"=>"preparing joint model","kind"=>kind,"status"=>"IN_PROGRESS_NOT_A_GATE")));flush(stdout)
    model=kind=="analytic" ? fixture_model() : JointLikelihoodReference.pinned_fixture_model(nativefile,partonfile)
    fixture=fixture_likelihood(model);l=fixture.likelihood;theta=fixture.theta
    report["model"]=joint_model_metadata(model);report["likelihood"]=joint_likelihood_metadata(l)
    report["parameters"]=theta
    report["data"]=Dict("point_ids"=>collect(l.point_ids),"observations"=>l.observations,"covariance"=>l.covariance,
        "additive"=>l.additive,"relative"=>l.relative,"nuisance_prior"=>l.nuisance_prior,"nuisance_ids"=>collect(l.nuisance_ids))
    checks=report["checks"]
    check(name,a,b,tol=p["algebra_max_scaled"];metric=:scaled)=push!(checks,compare_check(name,a,b,tol;metric=metric,floor=p["relative_floor"]))
    profile=joint_profile(l,theta)
    reference=augmented_reference(profile.predictions,l.observations,l.covariance,l.additive,l.relative,l.nuisance_prior)
    check("independent augmented QR profile",[profile.chi2],[reference.chi2])
    check("independent augmented QR nuisance",profile.nuisance,reference.nuisance)
    flipped=copy(theta);flipped[[1,2,7,8]].*=-1
    check("global transversity Collins sign ambiguity preserved",predict_joint(model,flipped),profile.predictions)
    check("global sign ambiguity in joint profile",[joint_objective(l,flipped)],[profile.chi2])
    design=l.additive+reshape(profile.predictions,:,1).*l.relative
    effective=l.covariance+design*l.nuisance_prior*transpose(design)
    check("effective covariance identity",[profile.chi2],[dot(profile.residual,effective\profile.residual)])
    check("original-coordinate nuisance stationarity",transpose(design)*(l.covariance\profile.shifted_residual),l.nuisance_prior\profile.nuisance)
    report["profile"]=Dict("chi2"=>profile.chi2,"predictions"=>profile.predictions,"nuisance"=>profile.nuisance,
        "residual"=>profile.residual,"shifted_residual"=>profile.shifted_residual,
        "independent_QR_chi2"=>reference.chi2,"independent_QR_nuisance"=>reference.nuisance)
    zero=JointGaussianLikelihood(model,l.point_ids,l.observations,l.covariance;scope="synthetic_closure",source="synthetic no-nuisance control")
    check("zero nuisance GLS reduction",[joint_objective(zero,theta)],[dot(profile.residual,l.covariance\profile.residual)])
    U=[1. .2;-.1 1.3];Ui=inv(U)
    rotated=JointGaussianLikelihood(model,l.point_ids,l.observations,l.covariance;additive=l.additive*U,relative=l.relative*U,
        nuisance_prior=Matrix(Symmetric(Ui*l.nuisance_prior*transpose(Ui))),nuisance_ids=["rotated 1","rotated 2"],
        scope="synthetic_closure",source="explicit synthetic nuisance reparameterization")
    check("nuisance-coordinate invariant profile",[joint_objective(rotated,theta)],[profile.chi2])
    check("nuisance-coordinate invariant physical pulls",U*joint_profile(rotated,theta).nuisance,profile.nuisance)
    indices=[7,1,3];sub=joint_subset(l,[l.point_ids[i] for i in indices];source="synthetic reordered marginal subset")
    subref=augmented_reference(profile.predictions[indices],l.observations[indices],l.covariance[indices,indices],
        l.additive[indices,:],l.relative[indices,:],l.nuisance_prior)
    check("reordered marginal subset profile",[joint_objective(sub,theta)],[subref.chi2])
    check("subset retains covariance",sub.covariance,l.covariance[indices,indices],0.;metric=:absolute)
    check("subset retains correlated nuisance prior",sub.nuisance_prior,l.nuisance_prior,0.;metric=:absolute)
    keep=findall(o->o.process!=:SIA,model.observations);holdout=joint_holdout(l,:SIA;source="synthetic process holdout")
    holdref=augmented_reference(profile.predictions[keep],l.observations[keep],l.covariance[keep,keep],l.additive[keep,:],l.relative[keep,:],l.nuisance_prior)
    check("process holdout marginal profile",[joint_objective(holdout,theta)],[holdref.chi2])
    data_design=l.additive+reshape(l.observations,:,1).*l.relative
    wrong_data=dot(profile.residual,(l.covariance+data_design*l.nuisance_prior*transpose(data_design))\profile.residual)
    wrong_diag=augmented_reference(profile.predictions,l.observations,Matrix(Diagonal(diag(l.covariance))),l.additive,l.relative,l.nuisance_prior).chi2
    marginal_logdet=logdet(Symmetric(effective))-logdet(Symmetric(l.covariance))
    report["counterfactual_values"]=Dict("profile"=>profile.chi2,"data_scaled"=>wrong_data,"diagonalized_C"=>wrong_diag,
        "marginal_relative_logdet"=>marginal_logdet)
    for (name,difference) in (("data-scaled nuisance is different",abs(wrong_data-profile.chi2)),
        ("discarded correlations are different",abs(wrong_diag-profile.chi2)),("marginal logdet is not profile chi2",abs(marginal_logdet)))
        check(name,[difference>p["relative_floor"] ? 0. : 1.],[0.],0.;metric=:absolute)
    end
    controls=Dict{String,Any}[]
    function reject(name,f)
        rejected=false;message=""
        try;f();catch exception;rejected=true;message=sprint(showerror,exception);end
        push!(controls,Dict("name"=>name,"rejected"=>rejected,"error"=>message))
    end
    for scope in ("production","experimental_candidate","historical_benchmark_candidate")
        reject("reject $scope",()->JointGaussianLikelihood(model,l.point_ids,l.observations,l.covariance;scope=scope,source="not admitted"))
    end
    reject("wrong point order",()->JointGaussianLikelihood(model,reverse(l.point_ids),l.observations,l.covariance;scope="synthetic_closure",source="wrong order"))
    reject("unsupported lattice recipe",()->JointObservation("lattice",(kind=:truncated_moment,);source="unsupported"))
    reject("duplicate source phase",()->JointObservation("phase",model.observations[1].plan;phase=:sin_phi_sum,source="wrong"))
    for field in (:observations,:covariance,:additive,:relative,:nuisance_prior,:nuisance_factor)
        damaged=deepcopy(l);a=getfield(damaged,field);a[1]=nextfloat(a[1])
        reject("mutated $field",()->joint_objective(damaged,theta))
    end
    report["negative_controls"]=controls
    check("all scope identity and mutation controls reject",[all(r["rejected"] for r in controls) ? 0. : 1.],[0.],0.;metric=:absolute)
    println(JSON.json(Dict("stage"=>"independent 13-parameter finite differences","kind"=>kind,"status"=>"IN_PROGRESS_NOT_A_GATE")));flush(stdout)
    direction=Float64.(ip["derivative_direction"])
    derivative=derivative_reference(l,theta,direction,p["finite_difference_steps"];bits=p["reference_bits"])
    gradient=ForwardDiff.gradient(t->joint_objective(l,t),theta)
    hessian=ForwardDiff.hessian(t->joint_objective(l,t),theta);hvp=hessian*direction
    check("complete 13-coordinate gradient",gradient,derivative.gradient,p["derivative_max_scaled"])
    check("complete 13-coordinate HVP",hvp,derivative.hvp,p["derivative_max_scaled"])
    buffer=zeros(13);value=joint_value_and_gradient!(buffer,theta,l)
    check("in-place gradient agreement",buffer,gradient,0.;metric=:absolute)
    check("in-place value agreement",[value],[profile.chi2],0.;metric=:absolute)
    report["derivatives"]=Dict("gradient"=>gradient,"hessian"=>hessian,"direction"=>direction,"hvp"=>hvp,
        "reference_gradient"=>derivative.gradient,"reference_hvp"=>derivative.hvp,"finite_difference_trials"=>derivative.trials)
    stable=source_hashes()==hashes && all(bytes2hex(sha256(read(path)))==sha for (path,sha) in dependencies)
    if kind=="inputs"
        stable &= CollinsFit.verified_hoppet_build(nativefile)==report["native_build"] && CollinsFit.verify_parton_build(partonfile)==report["parton_build"]
    end
    check("sources and build dependencies stable",[stable ? 0. : 1.],[0.],0.;metric=:absolute)
    report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    report
end

function main()
    kind=isempty(ARGS) ? "" : ARGS[1]
    (kind=="analytic" && length(ARGS)==2) || (kind=="inputs" && length(ARGS)==4) ||
        error("usage: check_joint_likelihood.jl analytic OUTPUT | inputs NATIVE_BUILD PARTON_BUILD OUTPUT")
    output=last(ARGS);ispath(output) && error("use a fresh report path; earlier evidence is retained")
    report=Dict{String,Any}("status"=>"FAIL","scope"=>"synthetic_closure","theory_scope"=>"reconstructed_theory",
        "input_kind"=>(kind=="inputs" ? "pinned_CT10_DSSV08_DSS14_native_rows" : "explicitly_synthetic_source_rows"),
        "kind"=>kind,"production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,
        "real_data_fit_implemented"=>false,"experimental_data_admitted"=>false,"checks"=>Dict{String,Any}[],"julia_version"=>string(VERSION))
    try;checks!(report,kind,kind=="inputs" ? ARGS[2] : nothing,kind=="inputs" ? ARGS[3] : nothing)
    catch exception;report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally;report["completed_at"]=string(now(UTC));writejson(output,report);end
    println(JSON.json(Dict("status"=>report["status"],"kind"=>kind,"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "error"=>get(report,"error",nothing),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
end

if abspath(PROGRAM_FILE)==(@__FILE__);main();end
