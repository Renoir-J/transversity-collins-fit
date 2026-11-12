"""Fixed-grid, fixed-nf LO collinear evolution reconstructed from KPSY15 Eq. 69–71.

The input is h1(x), or explicitly Hhat^(3)(z) in the KPSY homogeneous approximation.
It is not x*h1, a Trento Collins first moment, or a complete TMD evolution operator.
"""
struct CollinearLOOperator
    grid::Vector{Float64}
    matrix::Matrix{Float64}
    quantity::Symbol
    Q2::Float64
    Q02::Float64
    alpha0::Float64
    nf::Int
    coupling_model::Symbol
    lambda_GeV::Union{Nothing,Float64}
end

function collinear_grid(grid)
    grid isa AbstractVector || error("one-dimensional grid required")
    all(v -> v isa Union{AbstractFloat,Integer}, grid) || error("fixed numeric grid required; grid derivatives are not implemented")
    x = Float64.(collect(grid))
    length(x) >= 3 && all(isfinite, x) && 0 < first(x) && last(x) == 1 &&
        all(diff(x) .> 0) || error("grid must increase strictly from positive xmin to exactly one")
    x
end

"""Generator P, with d h/d ln(Q²) = alpha_s/(2*pi) P h and h(1)=0.

Piecewise-linear h(y) is integrated analytically after the plus subtraction:
P*h = CF*[2*int_x^1 (h(x/z)-h(x))/(1-z) dz + (3/2+2*log(1-x))*h(x)].
The endpoint row is the zero limit on the imposed h(1)=0 subspace.
"""
function collinear_lo_generator(grid)
    x = collinear_grid(grid)
    n = length(x)
    P = zeros(n, n)
    for i in 1:n-1
        xi = x[i]
        P[i, i] = 1.5 + 2log1p(-xi)
        # On the first cell, h(y)-h(xi) cancels y-xi exactly.
        weight = xi * log1p((x[i+1] - xi) / xi) / (x[i+1] - xi)
        P[i, i] -= 2weight
        P[i, i+1] += 2weight
        for j in i+1:n-1
            low, high = x[j], x[j+1]
            width = high - low
            ly = log1p(width / low)
            lr = log1p(width / (low - xi))
            P[i, i] += 2(ly - lr)
            P[i, j] += 2(-high * ly + (high - xi) * lr) / width
            P[i, j+1] += 2(low * ly + (xi - low) * lr) / width
        end
    end
    P .*= 4 / 3
    all(isfinite, P) || error("nonfinite collinear generator")
    P
end

"""Build a reusable evolution matrix. Grid, scales, alpha0 and nf are fixed inputs.

The scalar alpha0 constructor uses one-loop running with fixed integer nf.
An explicit AsymptoticTwoLoopCoupling selects the separate truncated two-loop
definition; neither constructor implements thresholds or process matching.
"""
function collinear_lo_operator(grid, Q2, Q02, alpha0; nf=3, quantity=:transversity)
    quantity in (:transversity, :KPSY_Hhat3_homogeneous) || error("unsupported collinear quantity/convention")
    nf isa Integer && 3 <= nf <= 6 || error("nf must be a fixed integer from three to six")
    all(v -> v isa Union{AbstractFloat,Integer}, (Q2, Q02, alpha0)) ||
        error("fixed numeric settings required; operator-setting derivatives are not implemented")
    all(isfinite, (Q2, Q02, alpha0)) || error("finite fixed coupling inputs required")
    x = collinear_grid(grid)
    alpha = alpha_lo(Q2, Q02, alpha0; nf=nf)
    tau = -2log(alpha / alpha0) / (11 - 2nf / 3)
    U = Matrix(exp(UpperTriangular(tau * collinear_lo_generator(x))))
    all(isfinite, U) || error("nonfinite collinear evolution operator")
    CollinearLOOperator(x, U, quantity, Q2, Q02, alpha0, nf, :one_loop, nothing)
end

"LO splitting with the explicitly specified asymptotic two-loop running coupling."
function collinear_lo_operator(grid, Q2, Q02, coupling::AsymptoticTwoLoopCoupling; quantity=:transversity)
    quantity in (:transversity, :KPSY_Hhat3_homogeneous) || error("unsupported collinear quantity/convention")
    x = collinear_grid(grid)
    tau = evolution_time(coupling, Q2, Q02)
    U = Matrix(exp(UpperTriangular(tau * collinear_lo_generator(x))))
    all(isfinite, U) || error("nonfinite collinear evolution operator")
    CollinearLOOperator(x, U, quantity, Q2, Q02, alpha_s(coupling,Q02), coupling.nf,
                        :asymptotic_two_loop, coupling.lambda_GeV)
end

function apply_operator(operator::CollinearLOOperator, initial::AbstractVector)
    length(initial) == length(operator.grid) || error("initial distribution/grid size mismatch")
    all(isfinite, initial) && iszero(last(initial)) || error("finite initial distribution with h(1)=0 required")
    operator.matrix * initial
end

function operator_adjoint(operator::CollinearLOOperator, cotangent::AbstractVector)
    length(cotangent) == length(operator.grid) && all(isfinite, cotangent) || error("invalid collinear cotangent")
    result = transpose(operator.matrix) * cotangent
    # The last ordinate is constrained to zero and is not a fitted coordinate.
    result[end] = zero(eltype(result))
    result
end
