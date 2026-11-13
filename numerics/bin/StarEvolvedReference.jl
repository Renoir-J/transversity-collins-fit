"Scalar reference assembly; shared fixed native rows and Born coefficients are explicit."
module StarEvolvedReference
using LinearAlgebra, QuadGK, SpecialFunctions
include("LeptonicReference.jl")
include("MatchedReference.jl")
export scalar_reference, polynomial_mellin, reference_single_sudakov, reference_matching,
       analytic_reference, analytic_pdf, analytic_ff, analytic_ff_bases

analytic_pdf(pid,x)=(pid==21 ? 3. : 1+.02pid)*(1-x)^3*(1+.1x)
analytic_ff(pid,z)=(pid==21 ? 2.3 : .8+.025pid)*(1-z)^3
const analytic_ff_bases=(1.1,.9,.4,.3)

function reference_single_sudakov(model,b,Q2;gluon=false)
    lo=log(4exp(-2Float64(Base.MathConstants.eulergamma))*(inv(b^2)+inv(model.bmax^2)))
    hi=log(Q2);nf=model.coupling.nf;b0=11-2nf/3;b1=102-38nf/3
    function integrand(t)
        L=t-2log(model.coupling.lambda_GeV);a=4/(b0*L)*(1-b1*log(L)/(b0^2*L))
        CF,CA=4/3,3.;color=gluon ? CA : CF
        A=color*a+color/2*(CA*(67/18-pi^2/6)-5nf/9)*a^2
        B=gluon ? -b0*a/2 : -3CF*a/2
        (A*(hi-t)+B)/2
    end
    quadgk(integrand,lo,hi;rtol=1e-12,atol=1e-13)[1]
end

"Independent integration in v rather than prepared log(y) matching weights."
function reference_matching(z,alpha,parton,quark,gluon;singlet=nothing)
    if parton==:quark
        f=v->((4/3)*((1-v)/2+(1+v^2)*log(v)/(1-v))*quark(z/v)+
            (4/3)*(v/2+(2/v-2+v)*log(v))*gluon(z/v))/v
        quark(z)+alpha/pi*quadgk(f,z,1.;rtol=1e-12,atol=1e-13)[1]
    elseif parton==:gluon
        singlet===nothing && error("gluon reference requires the complete quark singlet")
        f=v->(6*(v/(1-v)+(1-v)/v+v*(1-v))*log(v)*gluon(z/v)+
            .5*(v*(1-v)+(v^2+(1-v)^2)*log(v))*singlet(z/v))/v
        gluon(z)+alpha/pi*quadgk(f,z,1.;rtol=1e-12,atol=1e-13)[1]
    else
        error("reference parton must be quark/gluon")
    end
end

"No CollinsFit initial functions, tensor contraction or apply_operator call."
function initial_columns(samples,t)
    n=length(samples.grid);columns=size(samples.bases,1);out=zeros(eltype(t),n,columns)
    for j in 1:n-1,k in 1:columns
        x=samples.grid[j]
        if columns==2
            N,a,b=t[k],t[k+2],t[k+4]
            norm=exp((a+b)*log(a+b)-(iszero(a) ? zero(a) : a*log(a))-(iszero(b) ? zero(b) : b*log(b)))
            out[j,k]=N*x^a*(1-x)^b*norm*samples.bases[k,j]
        else
            i=k==1 ? 7 : 8
            out[j,k]=t[i]*x^t[i+2]*(1-x)^t[i+4]*samples.bases[k,j]
        end
    end
    out
end
function scalar_rows(operator,initial)
    out=zeros(eltype(initial),length(operator.queries),size(initial,2),length(operator.Q2))
    for (i,block) in enumerate(operator.blocks),s in eachindex(operator.Q2),c in axes(initial,2),j in eachindex(block.indices)
        out[i,c,s]+=block.rows[s,j]*initial[block.indices[j],c]
    end
    out
end
function collin_value(H,pid,charge,j)
    abs(pid)>=4 && return zero(H[1,1,1])
    abs(pid)==3 && return H[1,pid==3 ? 3 : 4,j]
    charge==0 && return (H[1,1,j]+H[1,2,j])/2
    H[1,(charge==1 ? pid in (2,-1) : pid in (-2,1)) ? 1 : 2,j]
end

"Prepare only fixed scalar-reference weights; retain native rows as shared evidence."
function scalar_reference(p;momentum_interval=nothing)
    f=p.fragmentation;m=f.model;r=f.radial;z=f.z;ids=f.pids
    common=[m.g2/2*log(hypot(1,b/m.bmax))*log(sqrt(f.Q2/m.Q02)) for b in r.b]
    q=[reference_single_sudakov(m,b,f.Q2)+common[j] for (j,b) in enumerate(r.b)]
    g=[reference_single_sudakov(m,b,f.Q2;gluon=true)+(9/4)*common[j] for (j,b) in enumerate(r.b)]
    shape=[(b/z)^2 for b in r.b]
    function kernel(n,b)
        momentum_interval===nothing && return besselj(n,f.jT*b/z)
        quadgk(jt->jt*besselj(n,jt*b/z),momentum_interval...;rtol=1e-12,atol=1e-14)[1]
    end
    J0=[kernel(0,b) for b in r.b];J1=[kernel(1,b) for b in r.b]
    D=[sum(r.db[j]*b*J0[j]*f.matched[i,j]*
        exp(-(pid==21 ? g[j] : q[j])-m.gh*shape[j])/(2pi*z^2) for (j,b) in enumerate(r.b)) for (i,pid) in enumerate(ids)]
    Cweights=[-r.db[j]*b^2*J1[j]*exp(-q[j])/(4pi*z^3) for (j,b) in enumerate(r.b)]
    function evaluate(t)
        h=scalar_rows(p.operator,initial_columns(p.samples,t))
        H=scalar_rows(f.operator,initial_columns(f.samples,t))
        C=[sum(Cweights[j]*collin_value(H,pid,f.charge,j)*exp(-(m.gh-t[13])*shape[j]) for j in eachindex(r.b)) for pid in ids]
        U,T=zero(t[1]),zero(t[1])
        for (j,born) in enumerate(p.born)
            k=born.kinematics;pref=p.weights[j]*p.couplings[j]^2/(k.S*k.x_a*k.x_b*(k.x_a*k.S+k.U))
            for term in born.terms
                U+=pref*term.unpolarized*p.fA[term.a,j]*p.fB[term.b,j]*D[term.c]
                pid=ids[term.a];pid in (2,1) || continue
                T+=pref*term.transfer*h[j,pid==2 ? 1 : 2,1]*p.fB[term.b,j]*C[term.c]
            end
        end
        (UU=U,UT=T)
    end
    evaluate
end

"Independent Mellin solution for x^a(1-x)^b, with beam or C11-Hhat3 matching."
function polynomial_mellin(x,tau,alpha;collins=false,a=1,b=3)
    0<x<1 && a>=0 && b isa Integer && b>=1 || error("interior beta-polynomial reference required")
    direction=cis(3pi/4)
    function integrand(t)
        N=1.5+t*direction
        digamma=MatchedReference.CollinearReference.contour_digamma(N+1)
        PN=(4/3)*(1.5-2*(digamma+Float64(Base.MathConstants.eulergamma)))
        moment=factorial(b)/prod(N+a+j for j in 0:b)
        regular=collins ? -(8/3)*MatchedReference.contour_trigamma(N+1) : zero(N)
        imag(direction*exp(-N*log(x)+tau*PN)*moment*(1+alpha/pi*regular))/pi
    end
    quadgk(integrand,0.,1.,5.,20.,Inf;rtol=1e-11,atol=1e-13)[1]
end

"Adaptive Fourier/Mellin and recoil-rapidity oracle for the declared polynomial family."
function analytic_reference(model,sqrt_s,pt,eta,z,jT,charge,fragmentation,t,hard;B=24.)
    t[3:6]==[1.,1.,2.,2.] && t[9:12]==[1.,1.,0.,0.] || error("reference requires the declared polynomial shapes")
    ids=(-5,-4,-3,-2,-1,1,2,3,4,5,21);qids=filter(!=(21),ids)
    nf=model.coupling.nf;lambda=model.coupling.lambda_GeV;b0=11-2nf/3;b1=102-38nf/3
    function alpha(q2)
        L=log(q2/lambda^2)
        4pi/(b0*L)*(1-b1*log(L)/(b0^2*L))
    end
    function time(q2)
        quadgk(v->alpha(exp(v))/(2pi),log(model.Q02),log(q2);rtol=1e-12,atol=1e-13)[1]
    end
    d0=Dict(pid=>analytic_ff(pid,z) for pid in ids)
    dq=Dict(pid=>reference_matching(z,1.,:quark,y->analytic_ff(pid,y),y->analytic_ff(21,y))-d0[pid] for pid in qids)
    dg=reference_matching(z,1.,:gluon,y->0.,y->analytic_ff(21,y);
        singlet=y->sum(analytic_ff(pid,y) for pid in qids))-d0[21]
    # One scalar common Collins transform; the four initial FF amplitudes are
    # attached only after independent evolution/matching of z*(1-z)^3.
    function radial(b)
        q2=4exp(-2Float64(Base.MathConstants.eulergamma))*(inv(b^2)+inv(model.bmax^2))
        a=alpha(q2);tau=time(q2)
        common=model.g2/2*log(hypot(1,b/model.bmax))*log(sqrt(pt^2/model.Q02))
        eq=reference_single_sudakov(model,b,pt^2)+common
        eg=reference_single_sudakov(model,b,pt^2;gluon=true)+(9/4)*common
        shape=(b/z)^2
        U=[pid==21 && fragmentation==:quark_only ? 0. :
            b*besselj(0,jT*b/z)*(d0[pid]+a*(pid==21 ? dg : dq[pid]))*
            exp(-(pid==21 ? eg : eq)-model.gh*shape)/(2pi*z^2) for pid in ids]
        H=polynomial_mellin(z,tau,a;collins=true)
        C=-b^2*besselj(1,jT*b/z)*H*exp(-eq-(model.gh-t[13])*shape)/(4pi*z^3)
        vcat(U,C)
    end
    transforms,radial_error=quadgk(radial,0.,.001,.02,.1,.4,1.,2.,4.,8.,B;rtol=1e-10,atol=1e-12)
    D=transforms[1:11];I=transforms[12]
    function C(pid)
        abs(pid)>=4 && return 0.
        pid==3 && return I*t[8]*analytic_ff_bases[3]
        pid==-3 && return I*t[8]*analytic_ff_bases[4]
        fav=t[7]*analytic_ff_bases[1];unfav=t[8]*analytic_ff_bases[2]
        charge==0 && return I*(fav+unfav)/2
        I*((charge==1 ? pid in (2,-1) : pid in (-2,1)) ? fav : unfav)
    end
    Cs=C.(ids);tau=time(pt^2);c=pt/sqrt_s
    low=-log(inv(c)-exp(-eta));high=log(inv(c)-exp(eta))
    function recoil(y)
        xa=c*(exp(eta)+exp(y));xb=c*(exp(-eta)+exp(-y))
        s=xa*xb*sqrt_s^2;tt=-xa*sqrt_s*pt*exp(-eta);u=-xb*sqrt_s*pt*exp(eta)
        h=polynomial_mellin(xa,tau,0.)
        U,T=0.,0.
        for (ia,a) in enumerate(ids),(ib,b) in enumerate(ids),(ic,cid) in enumerate(ids)
            cid==21 && fragmentation==:quark_only && continue
            HU,HT=hard(a,b,cid,s,tt,u)
            U+=HU*analytic_pdf(a,xa)*analytic_pdf(b,xb)*D[ic]
            ha=a==2 ? (27/4)*1.3t[1]*h : a==1 ? (27/4)*.8t[2]*h : 0.
            T+=HT*ha*analytic_pdf(b,xb)*Cs[ic]
        end
        [U,T]./(xa*xb)
    end
    value,recoil_error=quadgk(recoil,low,high;rtol=1e-10,atol=1e-12)
    factor=alpha(pt^2)^2/sqrt_s^4
    (values=(UU=factor*value[1],UT=factor*value[2]),
        radial_estimated_error=radial_error,recoil_estimated_error=recoil_error,
        radial_integrals=transforms,recoil_integrals=value,
        shared="Born hard coefficient definitions only; no prepared CSS/matching/Fourier/native evolution rows in the oracle")
end
end
