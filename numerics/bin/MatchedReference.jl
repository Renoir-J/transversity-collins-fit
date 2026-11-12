module MatchedReference

using QuadGK
include("CollinearReference.jl")
using .CollinearReference: contour_digamma
export contour_trigamma, matched_polynomial

"Stable continuation for the independent Collins Mellin moment."
function contour_trigamma(z::ComplexF64)
    imag(z)<0 && return conj(contour_trigamma(conj(z)))
    if real(z)<.5
        q=exp(2pi*im*z)
        return -4pi^2*q/(q-1)^2-contour_trigamma(1-z)
    end
    value=zero(z)
    while real(z)<24
        value+=inv(z^2)
        z+=1
    end
    value+inv(z)+inv(2z^2)+inv(6z^3)-inv(30z^5)+inv(42z^7)-inv(30z^9)+
        5/(66z^11)-691/(2730z^13)
end

"Analytic beta-polynomial moment, evolved and matched by Mellin inversion."
function matched_polynomial(x,tau,alpha,process,quantity; a=1,b=3,angle=3pi/4)
    0<x<1 || error("interior reference fraction required")
    ((process==:SIDIS && quantity in (:transversity,:KPSY_Hhat3_homogeneous)) ||
        (process==:SIA && quantity==:KPSY_Hhat3_homogeneous)) || error("unsupported matching reference")
    a>=0 && b isa Integer && b>=1 || error("analytic beta-polynomial input required")
    direction=cis(angle)
    # Independent constants: do not call the implementation's coefficient helper.
    delta=process==:SIDIS ? -8/3 : (pi^2-8)/3
    function integrand(t)
        N=1.5+t*direction
        PN=(4/3)*(1.5-2*(contour_digamma(N+1)+Float64(Base.MathConstants.eulergamma)))
        hN=factorial(b)/prod(N+a+j for j in 0:b)
        regular=quantity==:transversity ? zero(N) : -(8/3)*contour_trigamma(N+1)
        coefficient=1+(alpha/pi)*(delta+regular)
        imag(direction*exp(-N*log(x)+tau*PN)*hN*coefficient)/pi
    end
    quadgk(integrand,0,1,5,20,Inf;rtol=1e-11,atol=1e-13)[1]
end

end
