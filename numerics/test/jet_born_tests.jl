using Test, CollinsFit, ForwardDiff, QuadGK
include(joinpath(@__DIR__,"../bin/JetBornReference.jl"))
using .JetBornReference

@testset "Born hard factors against PDG angular cross sections and CHL spin transfer" begin
    channels=((:qqprime_t,:qqprime),(:qqprime_u,:zero),(:qq_same,:qq_same),
        (:qbarq_same_t,:qbarq_same),(:qbarq_same_u,:qbarq_opposite),(:qbarq_other,:zero),
        (:qbarq_gluon,:zero),(:qg_quark,:qg),(:qg_gluon,:zero),(:gq_quark,:zero),
        (:gq_gluon,:zero),(:gg_quark,:zero),(:gg_gluon,:zero))
    for (uk,tk) in channels,y in (1//32,1//8,1//4,1//2,3//4,7//8,31//32),scale in (1,17)
        s=big(scale)//1;t=-s*y;u=-s*(1-y)
        U=CollinsFit.jet_born_hard(uk,s,t,u);T=CollinsFit.jet_born_transfer(tk,s,t,u)
        reference=JetBornReference.reference_hard(uk,tk,y)
        @test U==reference.UU
        @test T==reference.UT
        @test U>0 && abs(T)<=U
    end
    @test CollinsFit.jet_born_transfer(:qq_same,4//1,-3//1,-1//1)==0
    @test CollinsFit.jet_born_transfer(:qq_same,8//1,-7//1,-1//1)<0
    @test CollinsFit.jet_born_transfer(:qbarq_opposite,2//1,-1//1,-1//1)==-8//27
    @test_throws ErrorException CollinsFit.jet_born_hard(:unknown,2.,-1.,-1.)
    @test_throws ErrorException CollinsFit.jet_born_transfer(:unknown,2.,-1.,-1.)
end

@testset "Physical jet delta reduction and tagged flavor routing" begin
    kin=JetBornKinematics(200.,20.,0.4,0.4)
    @test 0<kin.x_min<kin.x_a<1 && 0<kin.x_b<1
    @test isapprox(kin.s+kin.t+kin.u,0.;atol=1e-11)
    @test kin.x_b≈-kin.x_a*kin.T/(kin.x_a*kin.S+kin.U)
    @test jet_born_prefactor(kin,0.2)≈0.2^2/(kin.S*kin.x_a*kin.x_b*(kin.x_a*kin.S+kin.U))
    larger=JetBornKinematics(600.,60.,0.4,0.4)
    @test larger.x_b≈kin.x_b
    @test jet_born_prefactor(larger,0.2)≈jet_born_prefactor(kin,0.2)/3^4
    pids=[21,2,-2,1,-1,3,-3]
    full=JetBornPlan(pids,kin;fragmentation=:quark_and_gluon)
    quarks=JetBornPlan(pids,kin;fragmentation=:quark_only)
    @test length(full.terms)==127
    @test length(quarks.terms)==108
    @test count(t->t.transfer_kind!=:zero,full.terms)==48
    @test all(t->full.pids[t.a]!=21 && full.pids[t.c]!=21 && abs(full.pids[t.a])==abs(full.pids[t.c]),filter(t->t.transfer_kind!=:zero,full.terms))
    pids[2]=4
    @test full.pids==(21,2,-2,1,-1,3,-3)
    @test_throws Exception setfield!(full,:fragmentation,:quark_only)
    five=JetBornPlan([21,1,-1,2,-2,3,-3,4,-4,5,-5],kin;fragmentation=:quark_and_gluon)
    @test length(five.terms)==331
    @test count(t->t.transfer_kind!=:zero,five.terms)==120
    for bad in ([21,2],[21,2,-2,2],[0,2,-2],Any[21,true,-1],[21,6,-6],Int[])
        @test_throws ErrorException JetBornPlan(bad,kin;fragmentation=:quark_only)
    end
    @test_throws UndefKeywordError JetBornPlan([21,2,-2],kin)
    @test_throws ErrorException JetBornPlan([21,2,-2],kin;fragmentation=:automatic)
    for args in ((200.,20.,0.4,0.01),(200.,20.,0.4,1.),(200.,100.,0.,0.9),(200.,-20.,0.,0.4),(Inf,20.,0.,0.4),(200.,20.,Inf,0.4))
        @test_throws ErrorException JetBornKinematics(args...)
    end
    @test_throws ErrorException jet_born_prefactor(kin,0.)
    @test_throws ErrorException jet_born_prefactor(kin,true)
end

@testset "Born one-channel limits, phases and source-weighted integration" begin
    ids=[21,2,-2,1,-1,3,-3]
    first=JetBornPlan(ids,JetBornKinematics(200.,20.,0.,0.3);fragmentation=:quark_only)
    second=JetBornPlan(ids,JetBornKinematics(200.,30.,0.,0.5);fragmentation=:quark_only)
    n=length(ids);fA=zeros(n);fB=zeros(n);D=zeros(n);h=zeros(n);C=zeros(n)
    fA[2]=h[2]=D[2]=C[2]=fB[1]=1.
    value=jet_born_structure(first,fA,fB,D,h,C;alpha_s=0.2)
    y=-first.kinematics.t/first.kinematics.s
    @test jet_born_asymmetry(value;phase=:sin_phiS_minus_phiH)≈2(1-y)/(1+(1-y)^2)
    @test jet_born_asymmetry(value;phase=:sin_phiH_minus_phiS)≈-value.UT/value.UU
    doubled=jet_born_structure(first,fA,fB,D,h,C;alpha_s=0.4)
    @test doubled.UU≈4value.UU && doubled.UT≈4value.UT
    @test jet_born_asymmetry(doubled;phase=:sin_phiS_minus_phiH)≈value.UT/value.UU
    @test_throws ErrorException jet_born_asymmetry(value;phase=:DNN_removed)
    @test_throws ErrorException jet_born_asymmetry((UU=0.,UT=1.);phase=:sin_phiS_minus_phiH)
    @test_throws ErrorException jet_born_asymmetry((UU=-1.,UT=1.);phase=:sin_phiS_minus_phiH)
    @test_throws ErrorException jet_born_contractions(first,fA[2:end],fB,D,h,C)
    broken=copy(C);broken[2]=NaN
    @test_throws ErrorException jet_born_contractions(first,fA,fB,D,h,broken)
    arrays=map(v->hcat(v,v),(fA,fB,D,h,C))
    integrated=jet_born_integrals([first,second],[1.,2.],[0.2,0.17],arrays...;measure_source="synthetic fixed-node cross-section control")
    v2=jet_born_structure(second,fA,fB,D,h,C;alpha_s=0.17)
    @test integrated.UU≈value.UU+2v2.UU
    @test integrated.UT≈value.UT+2v2.UT
    @test !isapprox(integrated.UT/integrated.UU,(value.UT/value.UU+2v2.UT/v2.UU)/3;rtol=1e-12)
    @test_throws ErrorException jet_born_integrals([first,second],[1.,-1.],[0.2,0.17],arrays...;measure_source="synthetic")
    @test_throws ErrorException jet_born_integrals([first,second],[0.,0.],[0.2,0.17],arrays...;measure_source="synthetic")
    @test_throws ErrorException jet_born_integrals([first,second],[1.,2.],[0.2,0.17],arrays...;measure_source="")
    incompatible=JetBornPlan(ids,JetBornKinematics(510.,30.,0.,0.5);fragmentation=:quark_only)
    @test_throws ErrorException jet_born_integrals([first,incompatible],[1.,2.],[0.2,0.17],arrays...;measure_source="synthetic")
    # Channels with zero transversity transfer remain zero, independently of h/C values.
    hzero=zeros(n);hzero[1]=1.
    @test jet_born_contractions(first,ones(n),ones(n),ones(n),hzero,ones(n)).UT==0
    Cg=zeros(n);Cg[1]=1.
    @test jet_born_contractions(first,ones(n),ones(n),ones(n),ones(n),Cg).UT==0
end

@testset "Collins mass/z bridge and distinct differential versus integrated Gaussian controls" begin
    H=-0.07;z=0.3;width=0.16
    for mass in (0.134,0.13957,0.497)
        first_moment=trento_first_moment(H,z,mass)
        @test jet_collins_fourier_moment(H,z)≈mass*first_moment
        gaussian(j)=j/(z*mass)*(2z^2*mass^2*first_moment)*exp(-j^2/width)/(pi*width^2)
        expected(j)=-H*j*exp(-j^2/width)/(pi*width^2)
        @test all(isapprox(gaussian(j),expected(j);rtol=1e-12) for j in (0.,0.2,0.7,1.3))
        integral=quadgk(j->2pi*j*gaussian(j),0.,Inf;rtol=1e-12)[1]
        @test integral≈-sqrt(pi)*H/(2sqrt(width))
        @test gaussian(0.)==0 && integral!=0
    end
    @test ForwardDiff.derivative(h->jet_collins_fourier_moment(h,z),H)≈-1/(2z)
    @test_throws ErrorException jet_collins_fourier_moment(H,0.)
    @test_throws ErrorException jet_collins_fourier_moment(NaN,z)
end

@testset "Reference arithmetic rejects executable code" begin
    environment=Dict{Symbol,Any}(:s=>2//1,:t=>-1//1,:u=>-1//1)
    @test JetBornReference.arithmetic_expression("4*(s^2+u^2)/(9*t^2)",environment)==20//9
    for text in ("run(`echo forbidden`)","eval(1)","begin 1;2 end","s=1","Base.sin(s)","s^(1/2)")
        @test_throws Exception JetBornReference.arithmetic_expression(text,environment)
    end
end
