module BoundedJointReference
using CollinsFit, JSON, LinearAlgebra, ForwardDiff
include("JointLikelihoodReference.jl")
using .JointLikelihoodReference: augmented_reference
export closure_fixture, closure_information, augmented_reference

const CONTRACT=normpath(joinpath(@__DIR__,"../../config/bounded-joint-fit-contract.json"))

"Varied kinematics with explicitly synthetic fixed source rows, not native evolution."
function leptonic_fixture(process,x,z1,z2,Q2,P)
    model=KPSYCSS();radial=RadialRule([0.,.25,1.,4.,8.];order=4)
    plan=leptonic_plan(process,model,radial,Q2,P,z1;z2=process==:SIA ? z2 : nothing)
    first=[(1+.07i)*(1+.03b) for i in 1:10,b in radial.b]
    second=[(1+.11i)*(1+.05b) for i in 1:10,b in radial.b]
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    queries=process==:SIDIS ? [x] : [z1,z2]
    op=JointLikelihoodReference.synthetic_matching(process,kind,queries,radial,model)
    samples=PolarizedInitialSamples(kind,op.grid,
        [(1+.3i)*fraction*(1-fraction)^3 for i in 1:(process==:SIDIS ? 2 : 4),fraction in op.grid];input_kind=:analytic)
    op2=process==:SIDIS ? JointLikelihoodReference.synthetic_matching(process,:KPSY_Hhat3_homogeneous,[z1],radial,model) : nothing
    samples2=op2===nothing ? nothing : PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,
        [(1+.3i)*fraction*(1-fraction)^3 for i in 1:4,fraction in op2.grid];input_kind=:analytic)
    polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=samples2,
        target=process==:SIDIS ? NucleonTarget(:proton) : nothing,charges=process==:SIDIS ? (1,0) : (1,-1))
end

function closure_fixture(;boundary=false)
    policy=JSON.parsefile(CONTRACT)["validation"];geometry=policy["geometry"]
    observations=JointObservation[]
    convention=SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:DNN_removed,
        source="explicit synthetic closure phase and normalization")
    for x in geometry["SIDIS_x"],z in geometry["SIDIS_z"],charge in geometry["SIDIS_charges"]
        i=length(observations)+1;k=mod1(i,2)
        p=leptonic_fixture(:SIDIS,x,z,nothing,geometry["SIDIS_Q2"][k],geometry["SIDIS_P"][k])
        charge==1 || (p=with_pion_charges(p,charge))
        node=sidis_observable(p,convention)
        push!(observations,JointObservation("closure:SIDIS:$i",node;source="explicit synthetic x/z/charge closure point"))
    end
    angular=SIAMeanOfRatio(.6;source="synthetic angular mean",averaging_measure="fixed synthetic measure")
    for z1 in geometry["SIA_z1"],z2 in geometry["SIA_z2"],ratio in geometry["SIA_ratios"]
        i=length(observations)+1;k=mod1(i,2)
        p=leptonic_fixture(:SIA,nothing,z1,z2,geometry["SIA_Q2"][k],geometry["SIA_P"][k])
        node=sia_observable(p;ratio=Symbol(ratio),angular=angular,source="explicit synthetic SIA closure point")
        push!(observations,JointObservation("closure:SIA:$i",node;source="explicit synthetic z1/z2/ratio closure point"))
    end
    for energy in (200.,510.)
        p=JointLikelihoodReference.synthetic_star(energy)
        for momentum in (.2,.4)
            i=length(observations)+1
            push!(observations,JointObservation("closure:STAR:$i",with_momentum(p,momentum);
                phase=:sin_phiS_minus_phiH,source="explicit synthetic STAR closure point"))
        end
    end
    model=JointPredictionModel(observations;source="52-observation optimizer fixture; explicitly synthetic source rows and geometry")
    truth=Float64.(policy["closure_truth"])
    if boundary
        truth[1]=1.;truth[2]=-1.;truth[3]=0.;truth[5]=0.
    end
    y=predict_joint(model,truth);n=length(y)
    sigma=[.002+.00001i for i in 1:n]
    C=[sigma[i]*sigma[j]*.2^abs(i-j) for i in 1:n,j in 1:n]
    A=[.0001sin(i*j) for i in 1:n,j in 1:2];R=[.01cos(.3i+j) for i in 1:n,j in 1:2]
    likelihood=JointGaussianLikelihood(model,model.point_ids,y,C;additive=A,relative=R,
        nuisance_prior=[1. .35;.35 1.4],nuisance_ids=["synthetic closure shared 1","synthetic closure shared 2"],
        scope="synthetic_closure",source="noiseless truth with explicitly synthetic correlated covariance and theory-scaled nuisances")
    (likelihood=likelihood,truth=truth,starts=[Float64.(row) for row in policy["closure_starts"]],policy=policy)
end

"Local information at the interior truth, never global identifiability or interval coverage."
function closure_information(fixture)
    l=fixture.likelihood;theta=fixture.truth;p=predict_joint(l.model,theta)
    J=ForwardDiff.jacobian(t->predict_joint(l.model,t),theta)
    D=l.additive+reshape(p,:,1).*l.relative
    V=l.covariance+D*l.nuisance_prior*transpose(D)
    L=cholesky(Symmetric(V)).L
    scaled=(L\J).*reshape(Float64.(fixture.policy["parameter_scales"]),1,:)
    singular=svdvals(scaled);rank=count(>(first(singular)*1e-10),singular)
    Dict("scope"=>"synthetic_closure_local_information","point_count"=>length(p),"parameter_count"=>13,
        "rank"=>rank,"singular_values"=>singular,"relative_rank_threshold"=>1e-10,
        "full_local_rank"=>length(p)>13 && rank==13,"global_identifiability_verified"=>false,
        "interval_coverage_verified"=>false,"production_ready"=>false)
end
end
