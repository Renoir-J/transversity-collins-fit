const BOUNDED_JOINT_CONTRACT=normpath(joinpath(@__DIR__,"../../config/bounded-joint-fit-contract.json"))
const BOUNDED_SOLVER_COMPARISON_CONTRACT=normpath(joinpath(@__DIR__,"../../config/bounded-solver-comparison-contract.json"))

struct BoundedJointSettings
    algorithm::Symbol
    initial_step::NTuple{13,Float64}
    xtol_abs::NTuple{13,Float64}
    ftol_abs::Float64
    maxeval::Int
    maxtime::Float64
    trace::Bool
    source::String
    function BoundedJointSettings(;algorithm=:LN_SBPLX,initial_step=nothing,xtol_step_fraction=nothing,ftol_abs=nothing,
            maxeval=nothing,maxtime=0.,trace=true,source)
        algorithm isa Symbol && algorithm in (:LN_SBPLX,:LN_BOBYQA) || error("explicit supported bounded algorithm required")
        policy=algorithm==:LN_SBPLX ? JSON.parsefile(BOUNDED_JOINT_CONTRACT)["settings"] :
            JSON.parsefile(BOUNDED_SOLVER_COMPARISON_CONTRACT)["algorithms"]["LN_BOBYQA"]
        steps=initial_step===nothing ? policy["initial_step"] : initial_step
        steps isa Union{AbstractVector,Tuple} && length(steps)==13 || error("13 explicit initial steps required")
        owned=Tuple(fixed_scalar(v,"optimizer initial step") for v in steps)
        all(>(0),owned) || error("positive initial steps required")
        fraction=fixed_scalar(xtol_step_fraction===nothing ? policy["xtol_step_fraction"] : xtol_step_fraction,"optimizer xtol fraction")
        tolerance=fixed_scalar(ftol_abs===nothing ? policy["ftol_abs"] : ftol_abs,"optimizer ftol")
        time=fixed_scalar(maxtime,"optimizer resource time")
        fraction>0 && tolerance>0 && time>=0 || error("positive tolerances and nonnegative resource time required")
        budget=maxeval===nothing ? policy["default_maxeval"] : maxeval
        budget isa Integer && !(budget isa Bool) && 1<=budget<=typemax(Cint) || error("positive native-compatible evaluation budget required")
        trace isa Bool || error("explicit Bool trace setting required")
        absolute=Tuple(v*fraction for v in owned)
        all(v->isfinite(v) && v>0,absolute) || error("finite positive absolute coordinate tolerances required")
        new(algorithm,owned,absolute,tolerance,Int(budget),time,trace,joint_text(source,"optimizer settings source"))
    end
end

"The complete direct binary64 domain, including closed zero-shape/Nh boundaries."
function kpsy_parameter_bounds(gh)
    width=fixed_scalar(gh,"fixed fit gh");width>0 || error("positive gh required")
    (lower=(-1.,-1.,0.,0.,0.,0.,-Inf,-Inf,0.,0.,0.,0.,-Inf),
     upper=(1.,1.,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,prevfloat(width)))
end

function bounded_start(start,gh)
    start isa AbstractVector && length(start)==13 || error("explicit thirteen-coordinate start required")
    owned=[fixed_scalar(v,"fixed optimizer start") for v in start]
    KPSYParameters(owned;gh=gh)
    owned
end

function bounded_settings_metadata(s::BoundedJointSettings,gh)
    bounds=kpsy_parameter_bounds(gh)
    encode(v)=isinf(v) ? (v<0 ? "-Inf" : "Inf") : v
    Dict("algorithm"=>String(s.algorithm),
        "algorithm_contract_id"=>s.algorithm==:LN_SBPLX ? "KPSY-BOUNDED-JOINT-FIT-001" : "KPSY-BOUNDED-SOLVER-COMPARISON-001",
        "initial_step"=>collect(s.initial_step),"xtol_abs"=>collect(s.xtol_abs),
        "ftol_abs"=>s.ftol_abs,"ftol_rel"=>0.,"xtol_rel"=>0.,"maxeval"=>s.maxeval,"maxtime"=>s.maxtime,
        "trace"=>s.trace,"source"=>s.source,"parameter_names"=>string.(INITIAL_PARAMETER_NAMES),
        "lower_bounds"=>encode.(collect(bounds.lower)),"upper_bounds"=>encode.(collect(bounds.upper)),
        "coordinates"=>"direct physical binary64; no surrogate bounds or post-hoc clipping")
end

function bounded_runtime_metadata()
    entry=pathof(NLopt);package_root=dirname(dirname(entry))
    paths=[joinpath(directory,name) for (directory,_,names) in walkdir(joinpath(package_root,"src"))
        for name in names if endswith(name,".jl")]
    library=realpath(Libdl.dlpath(NLopt.libnlopt))
    Dict("wrapper_version"=>string(Base.pkgversion(NLopt)),"native_version"=>string(NLopt.version()),
        "library"=>library,"library_sha256"=>bytes2hex(sha256(read(library))),
        "wrapper_sources"=>Dict(path=>bytes2hex(sha256(read(path))) for path in paths),
        "manifest_sha256"=>bytes2hex(sha256(read(normpath(joinpath(@__DIR__,"../Manifest.toml"))))),
        "full_system_runtime_revalidated"=>false)
end

"One owned native optimizer. The callback is Julia and never requests a gradient."
function bounded_joint_start(f,start,s::BoundedJointSettings,gh)
    x0=bounded_start(start,gh);bounds=kpsy_parameter_bounds(gh)
    opt=NLopt.Opt(s.algorithm,13)
    algorithm_name=NLopt.algorithm_name(opt)
    name=lowercase(algorithm_name)
    identity=s.algorithm==:LN_SBPLX ? occursin("subplex",name) || occursin("sbplx",name) : occursin("bobyqa",name)
    identity || error("native algorithm identity mismatch")
    NLopt.lower_bounds!(opt,collect(bounds.lower));NLopt.upper_bounds!(opt,collect(bounds.upper))
    NLopt.initial_step!(opt,collect(s.initial_step));NLopt.xtol_abs!(opt,collect(s.xtol_abs))
    NLopt.xtol_rel!(opt,0.);NLopt.ftol_abs!(opt,s.ftol_abs);NLopt.ftol_rel!(opt,0.)
    NLopt.maxeval!(opt,s.maxeval);NLopt.maxtime!(opt,s.maxtime)
    evaluations=Ref(0);successful=Ref(0);trace=Dict{String,Any}[]
    best_value=Ref(Inf);best_theta=Ref{Union{Nothing,Vector{Float64}}}(nothing)
    function callback(x,gradient)
        evaluations[]+=1
        isempty(gradient) || error("bounded derivative-free solver requested an unsupported gradient")
        own=bounded_start(x,gh);before=copy(own)
        value=f(own)
        own==before || error("objective mutated parameter coordinates")
        value=fixed_scalar(value,"finite objective; penalties and nonfinite fallbacks are forbidden")
        successful[]+=1
        s.trace && push!(trace,Dict("evaluation"=>evaluations[],"theta"=>own,"chi2"=>value))
        if value<best_value[];best_value[]=value;best_theta[]=copy(own);end
        value
    end
    NLopt.min_objective!(opt,callback)
    result=Dict{String,Any}("status"=>"FAIL","start"=>x0,"algorithm_name"=>algorithm_name,
        "solver_stop"=>"NOT_RETURNED","optimizer_tolerance_stop"=>false,"trace"=>trace,
        "theta"=>nothing,"chi2"=>nothing,"postcheck_passed"=>false,"postcheck_evaluations"=>0)
    try
        value,theta,stop=NLopt.optimize(opt,x0)
        result["solver_stop"]=string(stop)
        if stop in (:MAXEVAL_REACHED,:MAXTIME_REACHED) && !isfinite(value) && successful[]==0
            result["status"]="INCOMPLETE"
            result["candidate_absent_reason"]="resource stop before any finite objective evaluation"
        else
        owned=bounded_start(theta,gh)
        isfinite(value) || error("nonfinite returned optimizer value")
        result["postcheck_evaluations"]=1
        rechecked=fixed_scalar(f(copy(owned)),"postchecked joint objective")
        difference=abs(value-rechecked)
        agrees=difference<=1e-11*max(abs(rechecked),1e-8)
        agrees || error("returned objective does not match independently re-evaluated candidate")
        result["theta"]=owned;result["chi2"]=rechecked;result["reported_chi2"]=value
        result["postcheck_difference"]=difference;result["postcheck_passed"]=true
        result["optimizer_tolerance_stop"]=stop in (:FTOL_REACHED,:XTOL_REACHED)
        result["status"]=result["optimizer_tolerance_stop"] ? "PASS" :
            stop in (:MAXEVAL_REACHED,:MAXTIME_REACHED) ? "INCOMPLETE" : "FAIL"
        end
    catch exception
        result["status"]="FAIL";result["error"]=sprint(showerror,exception,catch_backtrace())
    end
    result["solver_evaluations"]=NLopt.numevals(opt)
    result["callback_attempts"]=evaluations[];result["successful_evaluations"]=successful[]
    result["evaluation_accounting_matches"]=NLopt.numevals(opt)==evaluations[]
    result["retained_best_evaluation"]=best_theta[]===nothing ? nothing : Dict("theta"=>best_theta[],"chi2"=>best_value[])
    result["evaluation_accounting_matches"] || (result["status"]="FAIL")
    result
end

"Complete multistart job on prepared thirteen-parameter predictions; synthetic data only."
function fit_joint(likelihood::JointGaussianLikelihood,starts;settings::BoundedJointSettings)
    joint_likelihood_valid(likelihood)
    likelihood.scope=="synthetic_closure" || error("real-data fitting admission is not implemented")
    starts isa Union{AbstractVector,Tuple} && !isempty(starts) || error("nonempty explicit start vectors required")
    owned=[bounded_start(start,likelihood.model.gh) for start in starts]
    before=likelihood.signature;runtime=bounded_runtime_metadata();metadata=joint_likelihood_metadata(likelihood)
    candidates=Dict{String,Any}[]
    for (index,start) in enumerate(owned)
        candidate=bounded_joint_start(t->joint_objective(likelihood,t),start,settings,likelihood.model.gh)
        push!(candidates,candidate)
        @info "bounded joint start completed" start_index=index algorithm=String(settings.algorithm) status=candidate["status"] stop=candidate["solver_stop"] evaluations=candidate["solver_evaluations"]
    end
    stable=false;failure=nothing
    try
        joint_likelihood_valid(likelihood)
        stable=likelihood.signature==before && bounded_runtime_metadata()==runtime
    catch exception
        failure=sprint(showerror,exception,catch_backtrace())
    end
    status=(!stable || any(r->r["status"]=="FAIL",candidates)) ? "FAIL" :
        all(r->r["status"]=="PASS",candidates) ? "PASS" : "INCOMPLETE"
    finite_candidates=findall(r->r["postcheck_passed"]===true,candidates)
    best=isempty(finite_candidates) ? nothing : finite_candidates[argmin([candidates[i]["chi2"] for i in finite_candidates])]
    profile=nothing;report_attempts=0
    if stable && best!==nothing
        try
            report_attempts=1;profile=joint_profile(likelihood,candidates[best]["theta"])
        catch exception
            status="FAIL";failure=sprint(showerror,exception,catch_backtrace())
        end
    end
    result=Dict{String,Any}("status"=>status,"scope"=>"synthetic_closure","theory_scope"=>"reconstructed_theory",
        "contract_id"=>"KPSY-BOUNDED-JOINT-FIT-001","production_ready"=>false,"physical_extraction"=>false,
        "solver_contract_id"=>settings.algorithm==:LN_SBPLX ? "KPSY-BOUNDED-JOINT-FIT-001" : "KPSY-BOUNDED-SOLVER-COMPARISON-001",
        "G2PERF_pass"=>false,"real_data_fit_implemented"=>false,"experimental_data_admitted"=>false,
        "global_optimum_certified"=>false,"interval_coverage_verified"=>false,
        "settings"=>bounded_settings_metadata(settings,likelihood.model.gh),"likelihood"=>metadata,
        "runtime"=>runtime,"sources_and_runtime_stable"=>stable,"starts"=>candidates,
        "parameter_count"=>13,"point_count"=>length(likelihood.point_ids),"best_start_index"=>best,
        "best_candidate_status"=>best===nothing ? nothing : candidates[best]["status"],
        "report_profile_attempts"=>report_attempts,
        "theta"=>best===nothing ? nothing : candidates[best]["theta"],"chi2"=>profile===nothing ? nothing : profile.chi2,
        "predictions"=>profile===nothing ? nothing : profile.predictions,"nuisance"=>profile===nothing ? nothing : profile.nuisance,
        "limitations"=>["Tolerance stops do not certify global optimality or scientific admission",
            "Every requested start is retained; any failed/incomplete start prevents overall PASS",
            "Explicit synthetic observations only; no real-data loader, coverage, full-x lattice or new extraction",
            "Global transversity/Collins sign ambiguity and all historical/covariance/acceptance limitations remain"])
    failure===nothing || (result["source_error"]=failure)
    result
end
