module JointLikelihoodReference
using CollinsFit, LinearAlgebra
export fixture_model, fixture_likelihood, augmented_reference, derivative_reference

# Explicit synthetic rows: assembly/inference fixtures, never native evolution evidence.
function synthetic_matching(process,kind,queries,radial,model)
    grid=sort(unique(vcat(.05,.1,queries,.75,1.)))
    qs=css_mu_b2.(Ref(model),radial.b);blocks=CollinsFit.HoppetPointBlock[]
    for query in queries
        rows=zeros(length(qs),length(grid));rows[:,findfirst(==(query),grid)].=1 .+ .01 .*radial.b
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),rows))
    end
    data=(grid,collect(queries),qs,model.Q02,evolution_time.(Ref(model.coupling),qs,model.Q02),kind,
        5,32,:synthetic_fixture,blocks,sum(sizeof(b.rows) for b in blocks),0.)
    alpha=alpha_s.(Ref(model.coupling),qs)
    HoppetMatchingPlan(data...,process,alpha,CollinsFit.hoppet_point_snapshot(data...;process=process,alpha=alpha))
end

function synthetic_leptonic(process;model=KPSYCSS())
    radial=RadialRule([0.,.25,1.,4.,8.];order=4)
    plan=leptonic_plan(process,model,radial,process==:SIDIS ? 4. : 110.,.2,.3;z2=process==:SIA ? .5 : nothing)
    first=[(1+.07i)*(1+.03b) for i in 1:10,b in radial.b]
    second=[(1+.11i)*(1+.05b) for i in 1:10,b in radial.b]
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    op=synthetic_matching(process,kind,process==:SIDIS ? [.15] : [.3,.5],radial,model)
    samples=PolarizedInitialSamples(kind,op.grid,[(1+.3i)*x*(1-x)^3 for i in 1:(process==:SIDIS ? 2 : 4),x in op.grid];input_kind=:analytic)
    op2=process==:SIDIS ? synthetic_matching(process,:KPSY_Hhat3_homogeneous,[.3],radial,model) : nothing
    samples2=op2===nothing ? nothing : PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,
        [(1+.3i)*x*(1-x)^3 for i in 1:4,x in op2.grid];input_kind=:analytic)
    polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=samples2,
        target=process==:SIDIS ? NucleonTarget(:proton) : nothing,charges=process==:SIDIS ? (1,0) : (1,-1))
end

function synthetic_star(energy)
    model=KPSYCSS();radial=RadialRule([0.,.1,.5,1.,3.,8.];order=4);pt=20.
    ids=(CollinsFit.INITIAL_FLAVORS...,21);xmin=jet_born_xmin(energy,pt,.3);xa=xmin .+(1-xmin).*[.15,.4,.75]
    born=[JetBornPlan(ids,JetBornKinematics(energy,pt,.3,x);fragmentation=:quark_only) for x in xa]
    H=synthetic_matching(:STAR_C11,:KPSY_Hhat3_homogeneous,[.3],radial,model)
    Hs=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,H.grid,[(1+.2i)*(1-x)^3 for i in 1:4,x in H.grid];input_kind=:analytic)
    D=[pid==21 ? 0. : (1+.02i)*(1+.01b) for (i,pid) in enumerate(ids),b in radial.b]
    ff=star_fragmentation_plan(model,radial,pt^2,.3,.25,ids,D,H,Hs;fragmentation=:quark_only,gluon_model=:none,charge=1)
    grid=sort(unique(vcat(.05,xa,1.)));blocks=CollinsFit.HoppetPointBlock[]
    for x in xa
        row=zeros(1,length(grid));row[1,findfirst(==(x),grid)]=1+.02x
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),row))
    end
    data=(grid,xa,[pt^2],model.Q02,evolution_time.(Ref(model.coupling),[pt^2],model.Q02),
        :transversity,5,32,:synthetic_fixture,blocks,sum(sizeof(b.rows) for b in blocks),0.)
    h=HoppetPointPlan(data...,CollinsFit.hoppet_point_snapshot(data...))
    hs=PolarizedInitialSamples(:transversity,grid,[(1+.3i)*(1-x)^2 for i in 1:2,x in grid];input_kind=:analytic)
    A=[(1+.1i)*(1-x)^3 for i in 1:11,x in xa]
    B=[(1+.07i)*(1-b.kinematics.x_b)^3 for i in 1:11,b in born]
    star_point_plan(born,[.2,.3,.4],fill(alpha_s(model.coupling,pt^2),3),A,B,h,hs,ff;
        jet_radius=.6,axis=:standard,scale_prescription=:historical_pT_C11,measure_source="synthetic three-node dxa")
end

function assemble_fixture_model(sidis,sia,star,star510;source)
    convention=SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:DNN_removed,source="synthetic source phase")
    plus=sidis_observable(sidis,convention)
    minus=sidis_observable(with_pion_charges(sidis,-1),SIDISAsymmetryConvention(
        phase=:sin_phi_sum_minus_pi,normalization=:DNN_removed,source="synthetic shifted phase"))
    angular=SIAMeanOfRatio(.6;source="synthetic angular mean",averaging_measure="fixed synthetic measure")
    ul=sia_observable(sia;ratio=:UL,angular=angular,source="synthetic UL")
    uc=sia_observable(sia;ratio=:UC,angular=angular,momentum=MomentumRule(.05,.35;order=8),source="synthetic UC integral")
    sidisbin=leptonic_bin([plus,sidis_observable(with_momentum(sidis,.35),convention)],[.7,1.1];
        source="synthetic SIDIS bin",measure="two nodes",weight_definition="synthetic cross-section weights",acceptance="supplied synthetic weights")
    siabin=leptonic_bin([ul,sia_observable(with_momentum(sia,.35);ratio=:UL,angular=angular,source="synthetic UL node")],[.7 1.;1.1 .6];
        source="synthetic SIA bin",measure="two nodes",weight_definition="synthetic U/L cross-section weights",acceptance="supplied synthetic weights")
    starbin=star_bin([star510,with_momentum(star510,.4)],[.5,.7];phase=:sin_phiH_minus_phiS,measure_source="synthetic STAR bin")
    children=(plus,minus,ul,uc,sidisbin,star,starbin,siabin)
    observations=[JointObservation("synthetic:$i",p;phase=p isa StarPointPlan ? :sin_phiS_minus_phiH : nothing,source="synthetic observation $i") for (i,p) in enumerate(children)]
    JointPredictionModel(observations;source=source)
end

fixture_model()=assemble_fixture_model(synthetic_leptonic(:SIDIS),synthetic_leptonic(:SIA),
    synthetic_star(200.),synthetic_star(510.);source="mixed reconstructed process assembly with synthetic source rows")

function pinned_fixture_model(nativefile,partonfile)
    receipt=CollinsFit.verify_parton_build(partonfile);pdf=CT10Input(partonfile;policy=:reference)
    try
        hel=load_dssv08(receipt["input_root"]);ff=load_dss14(receipt["input_root"])
        model=KPSYCSS();radial=RadialRule([0.,.25,1.,4.,8.];order=4)
        sidis=prepare_polarized_sidis(pdf,hel,ff,nativefile,model,radial,4.,.2,.15,.3,1;
            matching_order=8,intervals=1024,range_policy=:reference)
        sia=prepare_polarized_sia(ff,nativefile,model,radial,110.,.2,.3,.5,1,-1;
            matching_order=8,intervals=1024,range_policy=:reference)
        stars=[prepare_star_point(pdf,hel,ff,nativefile,model,radial,energy,20.,.3,.3,.25,1;
            fragmentation=:quark_only,gluon_model=:none,jet_radius=.6,axis=:standard,
            scale_prescription=:historical_pT_C11,outer_order=4,matching_order=8,intervals=1024,
            range_policy=:reference) for energy in (200.,510.)]
        assemble_fixture_model(sidis,sia,stars...;source="pinned CT10/DSSV08/DSS14/native rows; synthetic kinematics, observations and bin weights")
    finally
        close(pdf)
    end
end

function fixture_likelihood(model=fixture_model())
    theta=collect(kpsy_published_parameters().values);theta[2]=-.8
    for i in (3,4,5,6,9,10,11,12);theta[i]+=.2;end
    t=predict_joint(model,theta);n=length(t)
    y=t+[.01*cos(i) for i in 1:n]
    sigma=[.02+.001i for i in 1:n]
    C=[sigma[i]*sigma[j]*.3^abs(i-j) for i in 1:n,j in 1:n]
    A=[.001*sin(i*j) for i in 1:n,j in 1:2]
    R=[.02*cos(.3i+j) for i in 1:n,j in 1:2]
    likelihood=JointGaussianLikelihood(model,model.point_ids,y,C;additive=A,relative=R,
        nuisance_prior=[1. .35;.35 1.4],nuisance_ids=["synthetic shared 1","synthetic shared 2"],
        scope="synthetic_closure",source="explicit synthetic covariance and source nuisances")
    (likelihood=likelihood,theta=theta)
end

"Independent augmented least squares; no profiled implementation or Woodbury solve."
function augmented_reference(prediction,observed,C,A,R,S)
    T=eltype(prediction);n=length(prediction);m=size(A,2)
    L=cholesky(Symmetric(Matrix{T}(C))).L
    residual=T.(observed)-prediction
    rhs=L\residual
    m==0 && return (chi2=sum(abs2,rhs),nuisance=T[])
    prior=cholesky(Symmetric(Matrix{T}(S))).L
    design=Matrix{T}(A)+reshape(prediction,:,1).*Matrix{T}(R)
    augmented=vcat(L\(design*prior),Matrix{T}(I,m,m));target=vcat(rhs,zeros(T,m))
    pulls=qr(augmented)\target
    (chi2=sum(abs2,target-augmented*pulls),nuisance=prior*pulls)
end

function derivative_reference(likelihood,theta,direction,steps;bits=256)
    setprecision(bits) do
        t=BigFloat.(theta);v=BigFloat.(direction)
        oracle=x->augmented_reference(predict_joint(likelihood.model,x),likelihood.observations,likelihood.covariance,
            likelihood.additive,likelihood.relative,likelihood.nuisance_prior).chi2
        rows=map(steps) do step
            h=BigFloat(step);gradient=BigFloat[];hvp=BigFloat[]
            gradient_values=Vector{BigFloat}[];hvp_values=Vector{BigFloat}[]
            for i in eachindex(t)
                e=zeros(BigFloat,length(t));e[i]=1
                g=[oracle(t+h*e),oracle(t-h*e)]
                q=[oracle(t+h*e+h*v),oracle(t+h*e-h*v),oracle(t-h*e+h*v),oracle(t-h*e-h*v)]
                push!(gradient_values,g);push!(hvp_values,q)
                push!(gradient,(g[1]-g[2])/(2h))
                push!(hvp,(q[1]-q[2]-q[3]+q[4])/(4h^2))
            end
            (step=step,gradient=gradient,hvp=hvp,gradient_values=gradient_values,hvp_values=hvp_values)
        end
        length(rows)==2 && rows[1].step==2rows[2].step || error("two declared halved finite-difference steps required")
        (gradient=Float64.((4rows[2].gradient-rows[1].gradient)/3),
            hvp=Float64.((4rows[2].hvp-rows[1].hvp)/3),
            trials=[Dict("step"=>r.step,"gradient_decimal"=>string.(r.gradient),"hvp_decimal"=>string.(r.hvp),
                "gradient_function_values_decimal"=>[string.(v) for v in r.gradient_values],
                "hvp_function_values_decimal"=>[string.(v) for v in r.hvp_values]) for r in rows])
    end
end
end
