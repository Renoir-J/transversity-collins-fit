"Independent scalar assembly and adaptive momentum reference; no observable implementation calls."
module ObservableReference
using CollinsFit, QuadGK, SpecialFunctions
using ..PolarizedReference: analytic_reference
using ..PolarizedReference.LeptonicReference: reference_sudakov
export prepare_generic_reference, mellin_observable_reference

const IDS=(-5,-4,-3,-2,-1,1,2,3,4,5)
const PAIRS=((1,1),(1,-1),(-1,1),(-1,-1))
swap(q)=abs(q)==1 ? 2sign(q) : abs(q)==2 ? sign(q) : q
charge2(q)=abs(q) in (2,4) ? 4/9 : 1/9
charge_value(f,q,c)=c==1 ? f(q) : c==-1 ? f(-q) : (f(q)+f(-q))/2
angular_y(y)=2(1-y)/(1+(1-y)^2)

function reference_weights(context,interval)
    plan=context.leptonic;z=plan.z1;n=plan.process==:SIDIS ? 1 : 2
    prefactor=plan.process==:SIDIS ? -inv(2z^3) : inv(4z^3*plan.z2)
    kernel(k,b)=interval===nothing ? besselj(k,plan.P*b/z) :
        quadgk(p->p*besselj(k,p*b/z),interval...;rtol=1e-12,atol=1e-14)[1]
    uu=[db*b*kernel(0,b)/(2pi*z^2) for (b,db) in zip(plan.radial.b,plan.radial.db)]
    pol=[prefactor*db*b^(n+1)*kernel(n,b)/(2pi) for (b,db) in zip(plan.radial.b,plan.radial.db)]
    uu,pol
end

"""
Prepare an independent fixed-row derivative reference. Native rows and fixed
source tables are shared evidence, not an independent evolution solution.
Scalar initial functions, charge loops, CSS and adaptive P integration are
separate from the optimized observable assembly. BigFloat-backed row products
use generic Julia arithmetic rather than Float64 tangent batching.
"""
function prepare_generic_reference(context;momentum_interval=nothing,y=.4,angular=.6)
    plan=context.leptonic;model=plan.model
    weights,polweights=reference_weights(context,momentum_interval)
    shape=[b^2*(plan.process==:SIDIS ? inv(plan.z1^2) : inv(plan.z1^2)+inv(plan.z2^2)) for b in plan.radial.b]
    fixed=[reference_sudakov(model,b,plan.Q2)+
        model.g2*log(sqrt(1+(b/model.bmax)^2))*log(sqrt(plan.Q2/model.Q02))+
        (plan.process==:SIDIS ? model.gq*b^2 : 0.) for b in plan.radial.b]
    damping=exp.(-fixed.-model.gh.*shape)
    index=Dict(q=>i for (i,q) in enumerate(plan.flavors))
    function denominator(c1,c2)
        sum(eachindex(plan.radial.b)) do j
            A(q)=context.unpolarized_first[index[q],j]
            D(q)=context.unpolarized_second[index[q],j]
            value=sum(charge2(q)*(plan.process==:SIDIS ?
                (context.target.unpolarized[1]*A(q)+context.target.unpolarized[2]*A(swap(q)))*charge_value(D,q,c1) :
                charge_value(A,q,c1)*charge_value(D,-q,c2)) for q in IDS)
            weights[j]*damping[j]*value
        end
    end
    denominators=plan.process==:SIDIS ? [denominator(context.charges...)] : [denominator(c...) for c in PAIRS]
    all(isfinite,denominators) && all(>(0),denominators) || error("positive reference charge denominators required")
    function initial_columns(samples,p)
        if samples.kind==:transversity
            [initial_transversity(p,pid,x,samples.bases[col,j]) for (j,x) in enumerate(samples.grid),(col,pid) in enumerate((2,1))]
        else
            [initial_collins(p,pid,1,x,samples.bases[:,j]) for (j,x) in enumerate(samples.grid),pid in (2,1,3,-3)]
        end
    end
    function prediction(theta)
        p=KPSYParameters(theta;gh=model.gh)
        first=apply_operator(context.first_operator,initial_columns(context.first_samples,p))
        second=plan.process==:SIDIS ? apply_operator(context.second_operator,initial_columns(context.second_samples,p)) : nothing
        damp=exp.(-fixed.-(model.gh-theta[13]).*shape)
        function H(array,query,q,c,j)
            abs(q)>=4 && return zero(eltype(array))
            abs(q)==3 && return array[query,q==3 ? 3 : 4,j]
            c==0 && return (array[query,1,j]+array[query,2,j])/2
            array[query,(c==1 ? q in (2,-1) : q in (-2,1)) ? 1 : 2,j]
        end
        function polarized(c1,c2)
            sum(eachindex(plan.radial.b)) do j
                hp(q)=q==2 ? first[1,1,j] : q==1 ? first[1,2,j] : zero(eltype(first))
                value=sum(charge2(q)*(plan.process==:SIDIS ?
                    (context.target.polarized[1]*hp(q)+context.target.polarized[2]*hp(swap(q)))*H(second,1,q,c1,j) :
                    H(first,1,q,c1,j)*H(first,2,-q,c2,j)) for q in IDS)
                polweights[j]*damp[j]*value
            end
        end
        if plan.process==:SIDIS
            ratio=polarized(context.charges...)/denominators[1]
            [ratio,angular_y(y)*ratio]
        else
            values=[polarized(c...) for c in PAIRS]
            U,L=values[2]+values[3],values[1]+values[4]
            DU,DL=denominators[2]+denominators[3],denominators[1]+denominators[4]
            angular.*[U/DU-L/DL,U/DU-(U+L)/(DU+DL)]
        end
    end
    (predict=prediction,denominators=denominators,uu_weights=weights,collins_weights=polweights,
     scope="generic_fixed_row_scalar_assembly_and_adaptive_P_reference")
end

"Mellin/adaptive values and amplitude/gc moments combined directly by paper equations."
function mellin_observable_reference(context,theta;momentum_interval=nothing,y=.4,angular=.6)
    plan=context.leptonic;indices=[1,2,7,8,13];a=theta[indices]
    reference(c1,c2)=analytic_reference(plan.model,plan.process,plan.Q2,plan.P,.1,plan.z1,plan.z2,
        c1,c2,context.target,a,last(plan.radial.edges);moments=true,momentum_interval=momentum_interval)
    if plan.process==:SIDIS
        r=reference(context.charges...);D=r.values.UU
        D>0 || error("positive Mellin SIDIS denominator required")
        factors=[1.,angular_y(y)]
        (values=factors.*r.values.Collins./D,
         gradient=reduce(vcat,[permutedims(f.*r.gradient./D) for f in factors]),
         hessian=cat([f.*r.hessian./D for f in factors]...;dims=3),
         charge_moments=[r])
    else
        r=[reference(c...) for c in PAIRS]
        DU=r[2].values.UU+r[3].values.UU;DL=r[1].values.UU+r[4].values.UU
        DU>0 && DL>0 || error("positive Mellin SIA category denominators required")
        # Fixed UU: each normalized observable is a linear combination of
        # the independent numerator values, gradients and Hessians.
        coeffs=[angular.*[-inv(DL),inv(DU),inv(DU),-inv(DL)],
                angular.*[-inv(DU+DL),inv(DU)-inv(DU+DL),inv(DU)-inv(DU+DL),-inv(DU+DL)]]
        (values=[sum(c[i]*r[i].values.Collins for i in 1:4) for c in coeffs],
         gradient=reduce(vcat,[permutedims(sum(c[i]*r[i].gradient for i in 1:4)) for c in coeffs]),
         hessian=cat([sum(c[i]*r[i].hessian for i in 1:4) for c in coeffs]...;dims=3),
         charge_moments=r)
    end
end
end
