"Massless observed-parton kinematics after integrating the unobserved momentum fraction."
struct JetBornKinematics
    sqrt_s::Float64
    pT::Float64
    eta::Float64
    x_a::Float64
    x_b::Float64
    x_min::Float64
    S::Float64
    T::Float64
    U::Float64
    s::Float64
    t::Float64
    u::Float64
    function JetBornKinematics(sqrt_s,pT,eta,x_a)
        rs=positive_fixed(sqrt_s,"jet sqrt_s");pt=positive_fixed(pT,"jet pT")
        rapidity=fixed_scalar(eta,"jet rapidity");xa=fixed_scalar(x_a,"jet x_a")
        xmin=jet_born_xmin(rs,pt,rapidity)
        xmin<xa<1 || error("jet x_a must lie strictly inside its physical integration interval")
        S=rs^2;T=-rs*pt*exp(-rapidity);U=-rs*pt*exp(rapidity)
        A=xa*S+U;xb=-xa*T/A
        s=xa*xb*S;t=xa*T;u=xb*U
        all(isfinite,(S,T,U,A,xb,s,t,u)) && S>0 && A>0 && 0<xb<1 && s>0 && t<0 && u<0 ||
            error("invalid massless jet kinematics")
        new(rs,pt,rapidity,xa,xb,xmin,S,T,U,s,t,u)
    end
end

function jet_born_xmin(sqrt_s,pT,eta)
    rs=positive_fixed(sqrt_s,"jet sqrt_s");pt=positive_fixed(pT,"jet pT")
    rapidity=fixed_scalar(eta,"jet rapidity");xT=2pt/rs
    denominator=2-xT*exp(-rapidity)
    isfinite(denominator) && denominator>0 || error("jet has no physical x interval")
    xmin=xT*exp(rapidity)/denominator
    isfinite(xmin) && 0<xmin<1 || error("jet has no physical x interval")
    xmin
end

function positive_fixed(value,label)
    x=fixed_scalar(value,label)
    x>0 || error("$label must be positive")
    x
end

# H_U is dimensionless: d sigma/dt = pi*alpha_s^2/s^2 * H_U.
# Keep rational constants exact for independent Rational/BigFloat references.
function jet_born_hard(kind,s,t,u)
    kind==:qqprime_t && return 4*(s^2+u^2)/(9t^2)
    kind==:qqprime_u && return 4*(s^2+t^2)/(9u^2)
    kind==:qq_same && return 4*((s^2+u^2)/t^2+(s^2+t^2)/u^2)/9-8s^2/(27t*u)
    kind==:qbarq_same_t && return 4*((s^2+u^2)/t^2+(t^2+u^2)/s^2)/9-8u^2/(27s*t)
    kind==:qbarq_same_u && return jet_born_hard(:qbarq_same_t,s,u,t)
    kind==:qbarq_other && return 4*(t^2+u^2)/(9s^2)
    kind==:qbarq_gluon && return 32*(t/u+u/t)/27-8*(t^2+u^2)/(3s^2)
    kind==:qg_quark && return (s^2+u^2)/t^2-4*(s/u+u/s)/9
    kind==:qg_gluon && return jet_born_hard(:qg_quark,s,u,t)
    kind==:gq_quark && return jet_born_hard(:qg_quark,s,u,t)
    kind==:gq_gluon && return jet_born_hard(:qg_quark,s,t,u)
    kind==:gg_quark && return (t/u+u/t)/6-3*(t^2+u^2)/(8s^2)
    kind==:gg_gluon && return 9*(3-t*u/s^2-s*u/t^2-s*t/u^2)/2
    error("unsupported observed-parton Born channel")
end

function jet_born_transfer(kind,s,t,u)
    kind==:zero && return zero(s+t+u)
    kind==:qg && return 8//9-2s*u/t^2
    kind==:qq_same && return 8*(-s*u/t^2+s/(3t))/9
    kind==:qbarq_same && return 8*(-s*u/t^2+u/(3t))/9
    kind==:qbarq_opposite && return -(8//27)*one(s+t+u)
    kind==:qqprime && return -8s*u/(9t^2)
    error("unsupported transverse-transfer Born channel")
end

function jet_born_channel(a,b,c)
    if a==21 && b==21
        return (c==21 ? :gg_gluon : :gg_quark,:zero)
    elseif b==21
        c==a && return (:qg_quark,:qg)
        c==21 && return (:qg_gluon,:zero)
    elseif a==21
        c==b && return (:gq_quark,:zero)
        c==21 && return (:gq_gluon,:zero)
    elseif a==b
        c==a && return (:qq_same,:qq_same)
    elseif a == -b
        c==a && return (:qbarq_same_t,:qbarq_same)
        c==b && return (:qbarq_same_u,:qbarq_opposite)
        c==21 && return (:qbarq_gluon,:zero)
        return (:qbarq_other,:zero)
    else
        c==a && return (:qqprime_t,:qqprime)
        c==b && return (:qqprime_u,:zero)
    end
    nothing
end

struct JetBornTerm
    a::Int
    b::Int
    c::Int
    unpolarized_kind::Symbol
    transfer_kind::Symbol
    unpolarized::Float64
    transfer::Float64
end

"Immutable fixed hard coefficients; this contains neither evolution nor experimental admission."
struct JetBornPlan
    kinematics::JetBornKinematics
    pids::Tuple{Vararg{Int}}
    fragmentation::Symbol
    terms::Tuple{Vararg{JetBornTerm}}
    function JetBornPlan(pids,kinematics::JetBornKinematics;fragmentation)
        fragmentation in (:quark_only,:quark_and_gluon) || error("explicit jet fragmentation species required")
        all(p->p isa Integer && !(p isa Bool) && (p==21 || 1<=abs(p)<=5),pids) || error("invalid jet PDG flavor")
        ids=Tuple(Int(p) for p in pids)
        length(ids)>=3 && length(unique(ids))==length(ids) && count(==(21),ids)==1 &&
            all(p==21 || -p in ids for p in ids) || error("jet basis requires gluon21 and complete quark/antiquark pairs")
        s,t,u=kinematics.s,kinematics.t,kinematics.u
        terms=JetBornTerm[]
        for a in eachindex(ids),b in eachindex(ids),c in eachindex(ids)
            fragmentation==:quark_only && ids[c]==21 && continue
            kinds=jet_born_channel(ids[a],ids[b],ids[c])
            kinds===nothing && continue
            uk,tk=kinds
            U=Float64(jet_born_hard(uk,s,t,u));T=Float64(jet_born_transfer(tk,s,t,u))
            isfinite(U) && U>0 && isfinite(T) || error("invalid fixed Born hard coefficient")
            push!(terms,JetBornTerm(a,b,c,uk,tk,U,T))
        end
        new(kinematics,ids,fragmentation,Tuple(terms))
    end
end

"Common alpha_s^2/S times the x_b delta Jacobian and both dx/x PDF factors."
function jet_born_prefactor(kinematics::JetBornKinematics,alpha_s)
    coupling=positive_fixed(alpha_s,"jet alpha_s")
    k=kinematics
    value=coupling^2/(k.S*k.x_a*k.x_b*(k.x_a*k.S+k.U))
    isfinite(value) && value>0 || error("invalid jet Born prefactor")
    value
end

"D is the TMD FF; C already equals jperp/(z*M_h)*H1perp in the Trento convention."
function jet_born_contractions(plan::JetBornPlan,fA,fB,D,hA,C)
    n=length(plan.pids)
    all(v->v isa AbstractVector && length(v)==n && all(isfinite,v),(fA,fB,D,hA,C)) ||
        error("jet responses must be finite flavor vectors in the declared PDG order")
    U=zero(fA[1]*fB[1]*D[1]);T=zero(hA[1]*fB[1]*C[1])
    for term in plan.terms
        U+=term.unpolarized*fA[term.a]*fB[term.b]*D[term.c]
        term.transfer_kind==:zero && continue
        T+=term.transfer*hA[term.a]*fB[term.b]*C[term.c]
    end
    (UU=U,UT=T)
end

function jet_born_structure(plan::JetBornPlan,fA,fB,D,hA,C;alpha_s)
    result=jet_born_contractions(plan,fA,fB,D,hA,C)
    prefactor=jet_born_prefactor(plan.kinematics,alpha_s)
    (UU=prefactor*result.UU,UT=prefactor*result.UT)
end

"Integrate supplied fixed cross-section/measure weights before forming an asymmetry."
function jet_born_integrals(plans,weights,couplings,fA,fB,D,hA,C;measure_source)
    measure_source isa AbstractString && !isempty(strip(measure_source)) || error("explicit jet measure source required")
    n=length(plans)
    n>0 && length(weights)==length(couplings)==n && all(p->p isa JetBornPlan,plans) || error("invalid jet node collection")
    first=plans[1];nf=length(first.pids)
    all(p->p.pids==first.pids && p.fragmentation==first.fragmentation && p.kinematics.sqrt_s==first.kinematics.sqrt_s,plans) ||
        error("jet nodes must share the declared flavor/model/beam-energy contract")
    all(v->v isa AbstractMatrix && size(v)==(nf,n) && all(isfinite,v),(fA,fB,D,hA,C)) ||
        error("jet responses must be flavor-by-node matrices")
    w=[fixed_scalar(value,"jet quadrature/measure weight") for value in weights]
    all(>=(0),w) && sum(w)>0 || error("jet measure weights must be nonnegative with positive sum")
    factors=[jet_born_prefactor(plans[i].kinematics,couplings[i])*w[i] for i in 1:n]
    U=zero(fA[1,1]*fB[1,1]*D[1,1]);T=zero(hA[1,1]*fB[1,1]*C[1,1])
    for i in 1:n
        value=jet_born_contractions(plans[i],view(fA,:,i),view(fB,:,i),view(D,:,i),view(hA,:,i),view(C,:,i))
        U+=factors[i]*value.UU;T+=factors[i]*value.UT
    end
    (UU=U,UT=T)
end

function jet_born_asymmetry(value;phase)
    phase in (:sin_phiS_minus_phiH,:sin_phiH_minus_phiS) || error("explicit supported jet azimuthal phase required")
    isfinite(value.UU) && value.UU>0 && isfinite(value.UT) || error("jet asymmetry requires a positive integrated unpolarized yield")
    (phase==:sin_phiS_minus_phiH ? 1 : -1)*value.UT/value.UU
end

"The dimensionful moment in the STAR J1 transform is M_h*H1perp(1), not Hhat3."
function jet_collins_fourier_moment(Hhat3,z)
    zz=css_z(z)
    isfinite(Hhat3) || error("non-finite Collins moment")
    -Hhat3/(2zz)
end
