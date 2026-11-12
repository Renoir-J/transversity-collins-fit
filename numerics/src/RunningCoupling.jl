"A fixed-nf asymptotic two-loop coupling, distinct from an exact two-loop RG solution."
struct AsymptoticTwoLoopCoupling
    lambda_GeV::Float64
    nf::Int
    function AsymptoticTwoLoopCoupling(lambda_GeV, nf)
        lambda = fixed_scalar(lambda_GeV, "Lambda")
        lambda > 0 || error("positive Lambda required")
        nf isa Integer && 3 <= nf <= 6 || error("fixed integer nf from three to six required")
        new(lambda, nf)
    end
end

AsymptoticTwoLoopCoupling(; lambda_GeV=0.225, nf=5) = AsymptoticTwoLoopCoupling(lambda_GeV, nf)

function fixed_scalar(value, name)
    value isa Union{AbstractFloat,Integer} && !(value isa Bool) ||
        error("fixed numeric $name required; setting derivatives are not implemented")
    converted = Float64(value)
    isfinite(converted) || error("finite $name required")
    converted
end

function coupling_log(coupling::AsymptoticTwoLoopCoupling, Q2)
    q2 = fixed_scalar(Q2, "Q2")
    q2 > 0 || error("positive Q2 required")
    L = log(q2) - 2log(coupling.lambda_GeV)
    L > 0 || error("Q2 is outside the coupling Landau-pole domain")
    L
end

"The specified 1/L - beta1/beta0^2*log(L)/L^2 expansion; no flavor matching."
function alpha_s(coupling::AsymptoticTwoLoopCoupling, Q2)
    L = coupling_log(coupling, Q2)
    beta0, beta1 = 11 - 2coupling.nf/3, 102 - 38coupling.nf/3
    alpha = 4pi/(beta0*L) * (1 - beta1/beta0^2 * log(L)/L)
    isfinite(alpha) && alpha > 0 || error("invalid asymptotic two-loop coupling")
    alpha
end

"Exact integral of this explicitly truncated coupling, integral alpha_s/(2*pi) d ln Q2."
function evolution_time(coupling::AsymptoticTwoLoopCoupling, Q2, Q02)
    L0 = coupling_log(coupling, Q02)
    coupling_log(coupling, Q2)
    beta0, beta1 = 11 - 2coupling.nf/3, 102 - 38coupling.nf/3
    # log1p retains precision for nearby scales and avoids an overflowing ratio.
    delta = abs(Q2-Q02) < Q02/2 ? log1p((Q2-Q02)/Q02) : log(Q2)-log(Q02)
    ratio = delta/L0
    logratio = log1p(ratio)
    correction = (logratio - (log(L0)+1)*ratio)/(L0*(1+ratio))
    tau = (2/beta0) * (logratio + beta1/beta0^2 * correction)
    isfinite(tau) || error("nonfinite integrated coupling")
    tau
end
