include("../bin/JointLikelihoodReference.jl")
using .JointLikelihoodReference

@testset "mixed reconstructed observations preserve identity, signs and source scope" begin
    model=fixture_model();fixture=fixture_likelihood(model);likelihood=fixture.likelihood;theta=fixture.theta
    @test length(model.observations)==8
    @test Set(o.process for o in model.observations)==Set((:SIDIS,:SIA,:STAR_200,:STAR_510))
    expected=[o.plan isa StarPointPlan ? predict_star(o.plan,theta;phase=o.phase) :
        o.plan isa StarBinPlan ? predict_star(o.plan,theta) : o.plan isa LeptonicBinPlan ? predict_bin(o.plan,theta) :
        predict_observable(o.plan,theta) for o in model.observations]
    @test predict_joint(model,theta)==expected
    flipped=copy(theta);flipped[[1,2,7,8]].*=-1
    @test predict_joint(model,flipped) ≈ expected rtol=1e-12 atol=1e-16
    @test joint_objective(likelihood,flipped) ≈ joint_objective(likelihood,theta) rtol=1e-12
    @test any(<(0),expected) && any(>(0),expected)
    @test predict_joint(model,KPSYParameters(theta))==expected
    metadata=joint_likelihood_metadata(likelihood)
    @test metadata["scope"]=="synthetic_closure" && metadata["theory_scope"]=="reconstructed_theory"
    @test metadata["production_ready"]===false && metadata["physical_extraction"]===false
    metadata["point_ids"][1]="changed copy"
    @test model.point_ids[1]=="synthetic:1"
    @test_throws ErrorException JointPredictionModel([model.observations[1],model.observations[1]];source="synthetic duplicate")
    @test_throws ErrorException JointObservation("lattice",(kind=:truncated_moment,);source="unsupported full-x fixture")
    @test_throws ErrorException JointObservation("phase",model.observations[6].plan;source="missing STAR phase")
    @test_throws ErrorException JointObservation("phase",model.observations[7].plan;phase=:sin_phiS_minus_phiH,source="duplicate bin phase")
    @test_throws ErrorException JointObservation("phase",model.observations[1].plan;phase=:sin_phi_sum,source="duplicate leptonic phase")
    @test_throws ErrorException JointObservation(" ",model.observations[1].plan;source="synthetic")
    different=JointLikelihoodReference.synthetic_leptonic(:SIDIS;model=KPSYCSS(g2=.7))
    different_node=sidis_observable(different;definition=:structure_function_ratio,source="synthetic changed physical model")
    @test_throws ErrorException JointPredictionModel([model.observations[1],JointObservation("different",different_node;source="synthetic")];source="mixed model")
    @test_throws ErrorException predict_joint(model,KPSYParameters(theta;gh=.05))
    @test_throws ErrorException predict_joint(model,[theta[1:12]...])
end

@testset "profile chi2 has independent QR/covariance identities and theory-scaled errors" begin
    f=fixture_likelihood();l=f.likelihood;t=f.theta;p=joint_profile(l,t)
    ref=augmented_reference(p.predictions,l.observations,l.covariance,l.additive,l.relative,l.nuisance_prior)
    @test p.chi2 ≈ ref.chi2 rtol=1e-11 atol=1e-14
    @test p.nuisance ≈ ref.nuisance rtol=1e-11 atol=1e-14
    design=l.additive+reshape(p.predictions,:,1).*l.relative
    effective=l.covariance+design*l.nuisance_prior*transpose(design)
    @test p.chi2 ≈ dot(p.residual,effective\p.residual) rtol=1e-11
    @test p.shifted_residual==p.residual-design*p.nuisance
    # Direct stationary equation in the original, correlated nuisance coordinates.
    @test transpose(design)*(l.covariance\p.shifted_residual) ≈ l.nuisance_prior\p.nuisance rtol=1e-10 atol=1e-13
    zero=JointGaussianLikelihood(l.model,l.point_ids,l.observations,l.covariance;scope="synthetic_closure",source="synthetic no nuisance")
    @test joint_objective(zero,t) ≈ dot(p.residual,l.covariance\p.residual) rtol=1e-12
    @test isempty(joint_profile(zero,t).nuisance)
    data_scaled=l.additive+reshape(l.observations,:,1).*l.relative
    wrong=dot(p.residual,(l.covariance+data_scaled*l.nuisance_prior*transpose(data_scaled))\p.residual)
    @test abs(p.chi2-wrong)>1e-8
    diagonal=JointGaussianLikelihood(l.model,l.point_ids,l.observations,Matrix(Diagonal(diag(l.covariance)));
        additive=l.additive,relative=l.relative,nuisance_prior=l.nuisance_prior,nuisance_ids=l.nuisance_ids,
        scope="synthetic_closure",source="synthetic deliberately discarded correlations")
    @test abs(p.chi2-joint_objective(diagonal,t))>1e-8
    # An invertible nuisance-coordinate change leaves the likelihood invariant.
    U=[1. .2;-.1 1.3];Ui=inv(U)
    rotated=JointGaussianLikelihood(l.model,l.point_ids,l.observations,l.covariance;
        additive=l.additive*U,relative=l.relative*U,nuisance_prior=Matrix(Symmetric(Ui*l.nuisance_prior*transpose(Ui))),
        nuisance_ids=["new synthetic coordinate 1","new synthetic coordinate 2"],scope="synthetic_closure",source="synthetic reparameterization")
    @test joint_objective(rotated,t) ≈ p.chi2 rtol=1e-11
    @test U*joint_profile(rotated,t).nuisance ≈ p.nuisance rtol=1e-11 atol=1e-14
end

@testset "subsets retain shared priors and cross-process covariance" begin
    f=fixture_likelihood();l=f.likelihood;t=f.theta
    ids=[l.point_ids[7],l.point_ids[1],l.point_ids[3]];indices=[7,1,3]
    sub=joint_subset(l,ids;source="synthetic reordered marginal subset")
    @test sub.point_ids==Tuple(ids)
    @test sub.covariance==l.covariance[indices,indices]
    @test sub.additive==l.additive[indices,:] && sub.relative==l.relative[indices,:]
    @test sub.nuisance_prior==l.nuisance_prior && sub.nuisance_ids==l.nuisance_ids
    @test predict_joint(sub.model,t)==predict_joint(l.model,t)[indices]
    ref=augmented_reference(predict_joint(l.model,t)[indices],l.observations[indices],l.covariance[indices,indices],
        l.additive[indices,:],l.relative[indices,:],l.nuisance_prior)
    @test joint_objective(sub,t) ≈ ref.chi2 rtol=1e-11
    holdout=joint_holdout(l,:SIA;source="synthetic SIA holdout")
    @test all(o->o.process!=:SIA,holdout.model.observations)
    @test length(holdout.point_ids)==5
    @test_throws ErrorException joint_holdout(l,:lattice;source="absent")
    @test_throws ErrorException joint_subset(l,[l.point_ids[1],l.point_ids[1]];source="duplicate")
    @test_throws ErrorException joint_subset(l,["absent"];source="unknown")
    @test_throws ErrorException joint_subset(l,String[];source="empty")
end

@testset "joint data fail closed on candidate scope, malformed inputs and mutation" begin
    f=fixture_likelihood();l=f.likelihood;t=f.theta
    make(y=l.observations,C=l.covariance;kwargs...)=JointGaussianLikelihood(l.model,l.point_ids,y,C;source="synthetic guard",kwargs...)
    for scope in ("experimental_candidate","historical_benchmark_candidate","production","reconstructed_theory",true,1,nothing)
        @test_throws ErrorException make(;scope=scope)
    end
    @test_throws UndefKeywordError make()
    @test_throws ErrorException JointGaussianLikelihood(l.model,reverse(l.point_ids),l.observations,l.covariance;scope="synthetic_closure",source="wrong order")
    @test_throws ErrorException make(l.observations[1:end-1];scope="synthetic_closure")
    @test_throws ErrorException make(fill(NaN,length(l.point_ids));scope="synthetic_closure")
    values=Any[l.observations...];values[1]=true
    @test_throws ErrorException make(values;scope="synthetic_closure")
    asym=copy(l.covariance);asym[1,2]+=.001
    @test_throws ErrorException make(l.observations,asym;scope="synthetic_closure")
    @test_throws ErrorException make(l.observations,zeros(size(l.covariance));scope="synthetic_closure")
    @test_throws ErrorException make(;scope="synthetic_closure",additive=l.additive)
    @test_throws ErrorException make(;scope="synthetic_closure",additive=l.additive,relative=l.relative,nuisance_ids=l.nuisance_ids,nuisance_prior=zeros(2,2))
    @test_throws ErrorException make(;scope="synthetic_closure",additive=l.additive,relative=l.relative,nuisance_ids=["same","same"],nuisance_prior=l.nuisance_prior)
    @test_throws ErrorException ForwardDiff.derivative(a->make(l.observations,l.covariance.*a;scope="synthetic_closure").covariance[1,1],1.)
    for field in (:observations,:covariance,:additive,:relative,:nuisance_prior,:nuisance_factor)
        owned=deepcopy(l);array=getfield(owned,field);array[1]=nextfloat(array[1])
        @test_throws ErrorException joint_objective(owned,t)
    end
    y=copy(l.observations);C=copy(l.covariance)
    owned=make(y,C;scope="synthetic_closure");before=joint_objective(owned,t)
    y[1]+=.1;C[1,1]+=.1
    @test joint_objective(owned,t)==before
    child=deepcopy(l);child.model.observations[1].plan.provenance["production_ready"]=true
    @test_throws ErrorException joint_objective(child,t)
end

@testset "joint likelihood full 13-coordinate AD and independent high-precision finite differences" begin
    f=fixture_likelihood();l=f.likelihood;t=f.theta
    direction=[.2,.1,.1,-.2,.02,.1,.2,-.1,.1,-.1,.1,-.2,.01]
    contract=JSON.parsefile(joinpath(@__DIR__,"../../config/joint-likelihood-contract.json"));p=contract["validation_policy"]
    reference=derivative_reference(l,t,direction,p["finite_difference_steps"];bits=p["reference_bits"])
    objective=x->joint_objective(l,x)
    gradient=ForwardDiff.gradient(objective,t);hvp=ForwardDiff.hessian(objective,t)*direction
    scaled(a,b)=maximum(abs.(a-b)./max.(abs.(b),p["relative_floor"]))
    @test scaled(gradient,reference.gradient)<=p["derivative_max_scaled"]
    @test scaled(hvp,reference.hvp)<=p["derivative_max_scaled"]
    buffer=zeros(13)
    @test joint_value_and_gradient!(buffer,t,l)==objective(t)
    @test buffer==gradient
    @test all(isfinite,gradient) && all(isfinite,hvp)
    boundary=copy(t);boundary[2]=-1.
    @test isfinite(objective(boundary)) && all(isfinite,ForwardDiff.gradient(objective,boundary))
end
