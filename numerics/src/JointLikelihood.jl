"Prepared observation families; unsupported observations are never dropped."
const JointChildPlan=Union{LeptonicObservablePlan,LeptonicBinPlan,StarPointPlan,StarBinPlan}

joint_child_model(p::LeptonicObservablePlan)=p.context.leptonic.model
joint_child_model(p::LeptonicBinPlan)=joint_child_model(first(p.nodes))
joint_child_model(p::StarPointPlan)=p.fragmentation.model
joint_child_model(p::StarBinPlan)=joint_child_model(first(p.points))
joint_child_process(p::LeptonicObservablePlan)=p.context.leptonic.process
joint_child_process(p::LeptonicBinPlan)=joint_child_process(first(p.nodes))
function joint_child_process(p::StarPointPlan)
    energy=first(p.born).kinematics.sqrt_s
    energy==200. ? :STAR_200 : energy==510. ? :STAR_510 : error("joint STAR process requires explicit 200/510 GeV geometry")
end
joint_child_process(p::StarBinPlan)=joint_child_process(first(p.points))
joint_child_valid(p::LeptonicObservablePlan)=observable_valid(p)
joint_child_valid(p::LeptonicBinPlan)=bin_valid(p)
joint_child_valid(p::StarPointPlan)=star_point_valid(p)
function joint_child_valid(p::StarBinPlan)
    foreach(star_point_valid,p.points)
    star_bin_signature(p)==p.signature || error("mutated joint STAR bin")
end
joint_common_settings(m)=(m.bmax,m.Q02,m.g2,m.gq,m.gh,m.coupling.lambda_GeV,m.coupling.nf)
joint_text(value,label)=value isa AbstractString && !isempty(strip(value)) ? String(value) : error("nonblank $label required")

"One signed, source-identified point/bin in a fixed common-parameter model."
struct JointObservation
    id::String
    plan::JointChildPlan
    phase::Union{Nothing,Symbol}
    process::Symbol
    source::String
    signature::String
    function JointObservation(id,plan;phase=nothing,source)
        plan isa JointChildPlan || error("unsupported joint observation; lattice requires a separately validated full-x plan")
        identifier=joint_text(id,"point identifier");origin=joint_text(source,"observation source")
        joint_child_valid(plan)
        if plan isa StarPointPlan
            phase in (:sin_phiS_minus_phiH,:sin_phiH_minus_phiS) || error("explicit STAR point phase required")
        else
            phase===nothing || error("phase already belongs to the prepared observable/bin; no second correction")
        end
        process=joint_child_process(plan)
        signature=numerical_snapshot((identifier,plan.signature,phase,process,origin))
        new(identifier,plan,phase,process,origin,signature)
    end
end

function joint_observation_valid(o::JointObservation)
    joint_child_valid(o.plan)
    numerical_snapshot((o.id,o.plan.signature,o.phase,o.process,o.source))==o.signature &&
        joint_child_process(o.plan)==o.process || error("mutated joint observation identity or phase")
end

struct JointPredictionModel
    observations::Tuple{Vararg{JointObservation}}
    point_ids::Tuple{Vararg{String}}
    settings::Tuple
    gh::Float64
    source::String
    signature::String
    function JointPredictionModel(observations;source)
        observations isa Union{AbstractVector,Tuple} && !isempty(observations) &&
            all(o->o isa JointObservation,observations) || error("nonempty prepared joint observations required")
        origin=joint_text(source,"joint model source");owned=Tuple(observations)
        foreach(joint_observation_valid,owned)
        ids=Tuple(o.id for o in owned)
        length(unique(ids))==length(ids) || error("duplicate joint point identifier")
        settings=joint_common_settings(joint_child_model(first(owned).plan))
        all(o->joint_common_settings(joint_child_model(o.plan))==settings,owned) || error("joint observations mix physical model settings")
        signature=numerical_snapshot((Tuple(o.signature for o in owned),ids,settings,origin))
        new(owned,ids,settings,settings[5],origin,signature)
    end
end

function joint_model_valid(model::JointPredictionModel)
    foreach(joint_observation_valid,model.observations)
    numerical_snapshot((Tuple(o.signature for o in model.observations),model.point_ids,model.settings,model.source))==model.signature ||
        error("mutated joint model")
end

function predict_joint(model::JointPredictionModel,p::KPSYParameters)
    joint_model_valid(model)
    p.gh==model.gh || error("joint parameter/width mismatch")
    [o.plan isa StarPointPlan ? predict_star(o.plan,p;phase=o.phase) :
        o.plan isa StarBinPlan ? predict_star(o.plan,p) :
        o.plan isa LeptonicBinPlan ? predict_bin(o.plan,p) : predict_observable(o.plan,p)
        for o in model.observations]
end
predict_joint(model::JointPredictionModel,theta::AbstractVector)=predict_joint(model,KPSYParameters(theta;gh=model.gh))

function joint_model_metadata(model::JointPredictionModel)
    joint_model_valid(model)
    Dict{String,Any}("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "G2PERF_pass"=>false,"contract_id"=>"KPSY-JOINT-PROFILE-001","point_ids"=>collect(model.point_ids),
        "processes"=>string.([o.process for o in model.observations]),"source"=>model.source,
        "common_model_settings"=>model.settings,"signature"=>model.signature)
end

function joint_fixed_matrix(value,label)
    value isa AbstractMatrix || error("matrix required for $label")
    [fixed_scalar(value[i,j],label) for i in axes(value,1),j in axes(value,2)]
end
function joint_positive_covariance(value,n,label)
    C=joint_fixed_matrix(value,label)
    size(C)==(n,n) && C==transpose(C) || error("matching symmetric $label required; no symmetrization")
    n==0 || isposdef(Symmetric(C)) || error("positive definite $label required; no clipping")
    C
end

"""Profiled Gaussian observations for joint-model validation.

The current public constructor admits explicit synthetic closure data only.
It never promotes historical/experimental candidates or production inputs.
"""
struct JointGaussianLikelihood
    model::JointPredictionModel
    point_ids::Tuple{Vararg{String}}
    observations::Vector{Float64}
    covariance::Matrix{Float64}
    additive::Matrix{Float64}
    relative::Matrix{Float64}
    nuisance_prior::Matrix{Float64}
    nuisance_factor::Matrix{Float64}
    nuisance_ids::Tuple{Vararg{String}}
    scope::String
    source::String
    signature::String
    function JointGaussianLikelihood(model::JointPredictionModel,point_ids,observations,covariance;
            additive=zeros(length(model.point_ids),0),relative=zeros(length(model.point_ids),0),
            nuisance_prior=zeros(0,0),nuisance_ids=String[],scope,source)
        scope isa AbstractString && scope=="synthetic_closure" ||
            error("joint likelihood validation currently requires explicit synthetic_closure observations; candidate/production data are not admitted")
        joint_model_valid(model);origin=joint_text(source,"observation data source")
        point_ids isa Union{AbstractVector,Tuple} || error("ordered point identifiers required")
        ids=Tuple(joint_text(id,"data point identifier") for id in point_ids)
        ids==model.point_ids || error("joint data/prediction point order mismatch")
        n=length(ids)
        observations isa AbstractVector && length(observations)==n || error("matching fixed observations required")
        y=[fixed_scalar(v,"fixed observation") for v in observations]
        C=joint_positive_covariance(covariance,n,"observation covariance")
        A=joint_fixed_matrix(additive,"additive nuisance design");R=joint_fixed_matrix(relative,"relative nuisance design")
        size(A,1)==size(R,1)==n && size(A,2)==size(R,2) || error("matching ordered nuisance columns required")
        labels=Tuple(joint_text(id,"nuisance identifier") for id in nuisance_ids)
        m=size(A,2)
        length(labels)==length(unique(labels))==m || error("unique nuisance identifiers matching columns required")
        S=joint_positive_covariance(nuisance_prior,m,"nuisance covariance")
        L=m==0 ? zeros(0,0) : Matrix(cholesky(Symmetric(S)).L)
        signature=numerical_snapshot((model.signature,ids,y,C,A,R,S,L,labels,scope,origin))
        new(model,ids,y,C,A,R,S,L,labels,String(scope),origin,signature)
    end
end

function joint_likelihood_valid(likelihood::JointGaussianLikelihood)
    joint_model_valid(likelihood.model)
    numerical_snapshot((likelihood.model.signature,likelihood.point_ids,likelihood.observations,
        likelihood.covariance,likelihood.additive,likelihood.relative,likelihood.nuisance_prior,
        likelihood.nuisance_factor,likelihood.nuisance_ids,likelihood.scope,likelihood.source))==likelihood.signature ||
        error("mutated joint observations, covariance, nuisance settings or source")
end

"Exact Gaussian nuisance profiling; B depends on the predicted theory, not measured y."
function joint_profile(likelihood::JointGaussianLikelihood,parameters)
    joint_likelihood_valid(likelihood)
    prediction=predict_joint(likelihood.model,parameters)
    all(isfinite,prediction) || error("nonfinite joint predictions")
    residual=likelihood.observations-prediction
    design=likelihood.additive+reshape(prediction,:,1).*likelihood.relative
    B=design*likelihood.nuisance_factor
    chi2,pulls=profile_nuisance(residual,likelihood.covariance,B)
    nuisance=likelihood.nuisance_factor*pulls
    (chi2=chi2,predictions=prediction,nuisance=nuisance,standardized_pulls=pulls,
        residual=residual,shifted_residual=residual-design*nuisance)
end
joint_objective(likelihood::JointGaussianLikelihood,parameters)=joint_profile(likelihood,parameters).chi2
function joint_value_and_gradient!(gradient,theta,likelihood::JointGaussianLikelihood)
    length(gradient)==length(theta)==13 || error("13 gradient coordinates required")
    ForwardDiff.gradient!(gradient,t->joint_objective(likelihood,t),theta)
    joint_objective(likelihood,theta)
end

function joint_likelihood_metadata(likelihood::JointGaussianLikelihood)
    joint_likelihood_valid(likelihood)
    Dict{String,Any}("scope"=>likelihood.scope,"theory_scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,
        "contract_id"=>"KPSY-JOINT-PROFILE-001","source"=>likelihood.source,
        "point_ids"=>collect(likelihood.point_ids),"nuisance_ids"=>collect(likelihood.nuisance_ids),
        "likelihood"=>"Gaussian nuisance-profile chi2; not marginalized posterior",
        "relative_design"=>"theory-scaled","signature"=>likelihood.signature,
        "model"=>joint_model_metadata(likelihood.model))
end

"Marginal data subset: retain its complete covariance and the shared nuisance prior."
function joint_subset(likelihood::JointGaussianLikelihood,point_ids;source)
    joint_likelihood_valid(likelihood)
    point_ids isa Union{AbstractVector,Tuple} && !isempty(point_ids) || error("nonempty explicit point selection required")
    ids=Tuple(joint_text(id,"selected point identifier") for id in point_ids)
    length(ids)==length(unique(ids)) || error("duplicate selected point")
    indices=[begin
        i=findfirst(==(id),likelihood.point_ids)
        i===nothing && error("unknown selected joint point")
        i
    end for id in ids]
    model=JointPredictionModel([likelihood.model.observations[i] for i in indices];source=source)
    JointGaussianLikelihood(model,ids,likelihood.observations[indices],likelihood.covariance[indices,indices];
        additive=likelihood.additive[indices,:],relative=likelihood.relative[indices,:],
        nuisance_prior=likelihood.nuisance_prior,nuisance_ids=likelihood.nuisance_ids,
        scope=likelihood.scope,source=source)
end

function joint_holdout(likelihood::JointGaussianLikelihood,process::Symbol;source)
    joint_likelihood_valid(likelihood)
    any(o->o.process==process,likelihood.model.observations) || error("requested holdout process is absent")
    ids=[o.id for o in likelihood.model.observations if o.process!=process]
    joint_subset(likelihood,ids;source=source)
end
