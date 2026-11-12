"Reconstructed KPSY Sudakov settings; no Wilson coefficients, PDFs/FFs or process fit."
struct KPSYCSS
    coupling::AsymptoticTwoLoopCoupling
    bmax::Float64
    Q02::Float64
    g2::Float64
    gq::Float64
    gh::Float64
    nodes::Vector{Float64}
    weights::Vector{Float64}
end

function KPSYCSS(; coupling=AsymptoticTwoLoopCoupling(), bmax=1.5, Q02=2.4,
                 g2=0.84, gq=0.106, gh=0.042, order=64)
    bmax, Q02, g2, gq, gh = [fixed_scalar(v, n) for (v,n) in
        ((bmax,"bmax"),(Q02,"Q02"),(g2,"g2"),(gq,"gq"),(gh,"gh"))]
    bmax > 0 && bmax*coupling.lambda_GeV < 1 || error("bmax must be positive and below 1/Lambda")
    coupling_log(coupling, Q02)
    g2 >= 0 && gq >= 0 && gh > 0 || error("invalid fixed nonperturbative widths")
    order isa Integer && 2 <= order <= 512 || error("quadrature order must be an integer from 2 to 512")
    k = collect(1:order-1)
    eig = eigen(SymTridiagonal(zeros(order), k ./ sqrt.(4 .* k.^2 .- 1)))
    KPSYCSS(coupling, bmax, Q02, g2, gq, gh, eig.values, 2 .* eig.vectors[1,:].^2)
end

function css_bstar(model::KPSYCSS, b)
    b = fixed_scalar(b, "b")
    b > 0 || error("positive b required; b=0 is not an evaluation node")
    inv(hypot(inv(b), inv(model.bmax)))
end

function css_mu_b2(model::KPSYCSS, b)
    value = (2exp(-Float64(Base.MathConstants.eulergamma))/css_bstar(model,b))^2
    isfinite(value) || error("nonfinite mu_b squared")
    value
end

"KPSY Eq.27: A through alpha_s^2 and B through alpha_s, with alpha_s/pi expansion."
function css_coefficients(alpha, nf)
    CF, CA, TR = 4/3, 3.0, 0.5
    A2 = CF/2 * (CA*(67/18-pi^2/6) - (10/9)*TR*nf)
    a = alpha/pi
    (CF*a + A2*a^2, -1.5CF*a)
end

function css_integral(lower2, upper2, coefficients, nodes, weights)
    lower, upper = log(lower2), log(upper2)
    half = (upper-lower)/2
    center = (upper+lower)/2
    total = sum(eachindex(nodes)) do i
        t = center + half*nodes[i]
        A, B = coefficients(exp(t))
        weights[i] * (A*(upper-t) + B)
    end
    value = half * total
    isfinite(value) || error("nonfinite perturbative Sudakov exponent")
    value
end

"Paper-defined signed limits; no cap mu_b<=Q and no clamp S_pert>=0 is inserted."
function sudakov_perturbative(model::KPSYCSS, b, Q2)
    coupling_log(model.coupling, Q2)
    lower2 = css_mu_b2(model,b)
    coefficients(q2) = css_coefficients(alpha_s(model.coupling,q2), model.coupling.nf)
    css_integral(lower2, Q2, coefficients, model.nodes, model.weights)
end

function css_np_common(model, b, Q2)
    css_bstar(model,b)
    coupling_log(model.coupling,Q2)
    model.g2 * 0.5log1p((b/model.bmax)^2) * 0.5log(Q2/model.Q02)
end

function css_z(z)
    value = fixed_scalar(z, "fragmentation z")
    0 < value <= 1 || error("fragmentation z must lie in (0,1]")
    value
end

function collins_width(model, gc)
    gc isa Real && isfinite(gc) && model.gh-gc > 0 ||
        error("finite gc with gh-gc>0 required for a convergent Collins Gaussian")
    model.gh-gc
end

finite_np(value) = isfinite(value) ? value : error("nonfinite nonperturbative exponent")

sidis_snp(model::KPSYCSS,b,Q2,z) = finite_np(css_np_common(model,b,Q2) +
    (model.gq+model.gh/css_z(z)^2)*b^2)

sidis_collins_snp(model::KPSYCSS,b,Q2,z,gc) = finite_np(css_np_common(model,b,Q2) +
    (model.gq+collins_width(model,gc)/css_z(z)^2)*b^2)

sia_snp(model::KPSYCSS,b,Q2,z1,z2) = finite_np(css_np_common(model,b,Q2) +
    model.gh*(inv(css_z(z1)^2)+inv(css_z(z2)^2))*b^2)

sia_collins_snp(model::KPSYCSS,b,Q2,z1,z2,gc) = finite_np(css_np_common(model,b,Q2) +
    collins_width(model,gc)*(inv(css_z(z1)^2)+inv(css_z(z2)^2))*b^2)
