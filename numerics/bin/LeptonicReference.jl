"""Independent variable transforms and analytic oracles; no CollinsFit numerical calls."""
module LeptonicReference
using QuadGK, SpecialFunctions
export reference_match, reference_sudakov, reference_structure, gaussian_hankel, angular_bessel

function reference_match(process,channel,x,alpha,quark,gluon)
    setprecision(256) do
        CF,TR = BigFloat(4)/3,BigFloat(1)/2
        low,a = BigFloat(x),BigFloat(alpha)/BigFloat(pi)
        delta = process == :SIDIS ? -2CF : CF*(BigFloat(pi)^2-8)/4
        result = (1+a*delta)*quark(low)
        (low==1 || channel==:transversity || alpha==0) && return Float64(result)
        function integrand(v)
            lr = log(v)/(1-v)
            if channel == :pdf
                rq,rg = CF*(1-v)/2,TR*v*(1-v)
            elseif channel == :ff
                rq = CF*((1-v)/2+(1+v^2)*lr)
                rg = CF*(v/2+(2/v-2+v)*log(v))
            else
                rq,rg = 2CF*v*lr,zero(v)
            end
            (rq*quark(low/v)+rg*gluon(low/v))/v
        end
        Float64(result+a*quadgk(integrand,low,one(low);rtol=big"1e-28",atol=big"1e-32")[1])
    end
end

function reference_sudakov(model,b,Q2)
    lower = log(4exp(-2Float64(Base.MathConstants.eulergamma))*(1/b^2+1/model.bmax^2))
    upper = log(Q2)
    nf,Lambda = model.coupling.nf,model.coupling.lambda_GeV
    b0,b1 = 11-2nf/3,102-38nf/3
    function integrand(t)
        L = t-2log(Lambda)
        a = 4/(b0*L)*(1-b1*log(L)/(b0^2*L))
        A = (4/3)*a+(2/3)*(3*(67/18-pi^2/6)-5nf/9)*a^2
        A*(upper-t)-2a
    end
    quadgk(integrand,lower,upper;rtol=1e-12,atol=1e-13)[1]
end

gaussian_hankel(n,k,A) = k^n*exp(-k^2/(4A))/(2^(n+1)*A^(n+1)*2pi)

function angular_bessel(n,t)
    # Real angular projections of exp(i*t*cos(phi)); no Bessel function call.
    f = n==0 ? (phi -> cos(t*cos(phi))) : n==1 ?
        (phi -> cos(phi)*sin(t*cos(phi))) : (phi -> -cos(2phi)*cos(t*cos(phi)))
    quadgk(f,0.,2pi;rtol=1e-13,atol=1e-13)[1]/(2pi)
end

function reference_structure(process,model,Q2,P,z1,z2,ids,first,second,pol_first,pol_second,gc,B)
    # Scalar dictionaries keyed by physical PDG flavor; use no precomputed plan.
    function integrand(b,polarized)
        A,D = polarized ? (pol_first,pol_second) : (first,second)
        amplitude = sum((abs(q) in (2,4) ? 4/9 : 1/9)*A(q,b)*D(process==:SIA ? -q : q,b) for q in ids)
        common = model.g2*log(sqrt(1+(b/model.bmax)^2))*log(sqrt(Q2/model.Q02))
        width = polarized ? model.gh-gc : model.gh
        gaussian = process==:SIDIS ? model.gq+width/z1^2 : width*(1/z1^2+1/z2^2)
        n = !polarized ? 0 : process==:SIDIS ? 1 : 2
        prefactor = !polarized ? inv(z1^2) : process==:SIDIS ? -inv(2z1^3) : inv(4z1^3*z2)
        prefactor*amplitude*b^(n+1)*besselj(n,P*b/z1)*
            exp(-reference_sudakov(model,b,Q2)-common-gaussian*b^2)/(2pi)
    end
    (UU=quadgk(b -> integrand(b,false),0.,.001,.1,1.,4.,B;rtol=1e-11,atol=1e-13)[1],
     Collins=quadgk(b -> integrand(b,true),0.,.001,.1,1.,4.,B;rtol=1e-11,atol=1e-13)[1])
end
end
