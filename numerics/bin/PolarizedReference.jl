"Independent Mellin/adaptive reference for a contracted analytic parameter family."
module PolarizedReference
using QuadGK, SpecialFunctions, ForwardDiff
include("MatchedReference.jl")
include("LeptonicReference.jl")
using .MatchedReference
using .LeptonicReference: reference_sudakov
export analytic_reference, analytic_amplitude

const IDS=(-5,-4,-3,-2,-1,1,2,3,4,5)
const FF_BASES=(1.1,.9,.4,.3)
swap(q)=abs(q)==1 ? 2sign(q) : abs(q)==2 ? sign(q) : q
charge_value(f,q,c)=c==1 ? f(q) : c==-1 ? f(-q) : (f(q)+f(-q))/2

function analytic_amplitude(process,charge1,charge2,weights,a)
    # a contains Nh_u, Nh_d, Nc_u, Nc_d, gc. Shapes are fixed at
    # (a_u,a_d,b_u,b_d,alpha_u,alpha_d,beta_u,beta_d)=(1,1,2,2,1,1,0,0).
    hp(q)=q==2 ? (27/4)*1.3a[1] : q==1 ? (27/4)*.8a[2] : zero(a[1])
    function H(q,c)
        abs(q)>=4 && return zero(a[3])
        q==3 && return FF_BASES[3]*a[4]
        q==-3 && return FF_BASES[4]*a[4]
        fav,unfav=FF_BASES[1]*a[3],FF_BASES[2]*a[4]
        c==0 && return (fav+unfav)/2
        (c==1 ? q in (2,-1) : q in (-2,1)) ? fav : unfav
    end
    sum((abs(q) in (2,4) ? 4/9 : 1/9)*
        (process==:SIDIS ? (weights[1]*hp(q)+weights[2]*hp(swap(q)))*H(q,charge1) : H(q,charge1)*H(-q,charge2)) for q in IDS)
end

function analytic_reference(model,process,Q2,P,x,z1,z2,c1,c2,target,amplitudes,B;moments=false,momentum_interval=nothing)
    cache=Dict{Float64,Tuple{Float64,Float64}}()
    # Fubini: adaptively integrate P*dP for each Bessel kernel before the
    # independent adaptive b integral. No production momentum/radial weights.
    momentum_cache=Dict{Tuple{Int,Float64},Float64}()
    bounds=if momentum_interval===nothing
        nothing
    else
        length(momentum_interval)==2 || error("two momentum bounds required")
        low,high=Float64.(momentum_interval)
        isfinite(low) && isfinite(high) && 0<=low<high || error("finite ordered nonnegative momentum bounds required")
        (low,high)
    end
    function kernel(n,b)
        bounds===nothing && return besselj(n,P*b/z1)
        get!(momentum_cache,(n,b)) do
            quadgk(p->p*besselj(n,p*b/z1),bounds...;rtol=1e-12,atol=1e-14)[1]
        end
    end
    b0,b1=11-2model.coupling.nf/3,102-38model.coupling.nf/3
    function legs(b)
        get!(cache,b) do
            mu2=4exp(-2Float64(Base.MathConstants.eulergamma))*(inv(b^2)+inv(model.bmax^2))
            L=log(mu2/model.coupling.lambda_GeV^2)
            alpha=4pi/(b0*L)*(1-b1*log(L)/(b0^2*L))
            tau=quadgk(log(model.Q02),log(mu2);rtol=1e-12,atol=1e-13) do t
                LL=t-2log(model.coupling.lambda_GeV)
                2/(b0*LL)*(1-b1*log(LL)/(b0^2*LL))
            end[1]
            first=matched_polynomial(process==:SIDIS ? x : z1,tau,alpha,process,
                process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous)
            second=matched_polynomial(process==:SIDIS ? z1 : z2,tau,alpha,process,:KPSY_Hhat3_homogeneous)
            first,second
        end
    end
    A(q)=1+.1q;D(q)=.8-.05q
    flavor_uu=sum((abs(q) in (2,4) ? 4/9 : 1/9)*
        (process==:SIDIS ? (target.unpolarized[1]*A(q)+target.unpolarized[2]*A(swap(q)))*charge_value(D,q,c1) :
            charge_value(A,q,c1)*charge_value(D,-q,c2)) for q in IDS)
    weights=process==:SIDIS ? target.polarized : (0.,0.)
    amplitude=analytic_amplitude(process,c1,c2,weights,amplitudes)
    shape=process==:SIDIS ? inv(z1^2) : inv(z1^2)+inv(z2^2)
    common(b)=model.g2*log(sqrt(1+(b/model.bmax)^2))*log(sqrt(Q2/model.Q02))+
        (process==:SIDIS ? model.gq*b^2 : 0.)+reference_sudakov(model,b,Q2)
    UU=flavor_uu*quadgk(b->b*kernel(0,b)*exp(-common(b)-model.gh*shape*b^2)/(2pi*z1^2),
        0.,.001,.1,1.,4.,B;rtol=1e-11,atol=1e-13)[1]
    n=process==:SIDIS ? 1 : 2
    prefactor=process==:SIDIS ? -inv(2z1^3) : inv(4z1^3*z2)
    function integral(power)
        quadgk(0.,.001,.1,1.,4.,B;rtol=1e-11,atol=1e-13) do b
            M1,M2=legs(b)
            prefactor*b^(n+1)*kernel(n,b)*M1*M2*(shape*b^2)^power*
                exp(-common(b)-(model.gh-amplitudes[5])*shape*b^2)/(2pi)
        end[1]
    end
    I=integral(0)
    result=(UU=UU,Collins=amplitude*I)
    moments || return result
    I1,I2=integral(1),integral(2)
    f(a)=analytic_amplitude(process,c1,c2,weights,a)
    g=ForwardDiff.gradient(f,amplitudes);H=ForwardDiff.hessian(f,amplitudes)
    J=g*I;J[5]=amplitude*I1
    Hess=H*I
    Hess[1:4,5]=g[1:4]*I1;Hess[5,1:4]=g[1:4]*I1;Hess[5,5]=amplitude*I2
    (values=result,gradient=J,hessian=Hess)
end
end
