module CollinearReference

using QuadGK, SpecialFunctions
export validation_grid, interpolate_grid, numerical_moment, contour_digamma,
       mellin_polynomial, mellin_polynomial_tau, max_scaled

function validation_grid(n, xmin)
    nlog = (n - 1) ÷ 4 + 1
    low = exp.(range(log(xmin), log(0.02); length=nlog))
    high = 1 .- range(sqrt(0.98), 0; length=n-nlog+1).^2
    low[1], high[1], high[end] = xmin, 0.02, 1.0
    vcat(low[1:end-1], high)
end

function interpolate_grid(grid, values, x)
    first(grid) <= x <= last(grid) || error("no reference extrapolation")
    i = min(searchsortedlast(grid, x), length(grid)-1)
    fraction = (x - grid[i]) / (grid[i+1] - grid[i])
    (1-fraction) * values[i] + fraction * values[i+1]
end

function numerical_moment(grid, values, n)
    sum(quadgk(x -> x^(n-1) * (values[i] + (values[i+1]-values[i]) *
            (x-grid[i]) / (grid[i+1]-grid[i])), grid[i], grid[i+1]; atol=1e-15)[1]
        for i in 1:length(grid)-1)
end

function contour_digamma(z)
    # Stable reflection in the upper half-plane avoids trigonometric overflow.
    imag(z) < 0 && return conj(contour_digamma(conj(z)))
    if real(z) < 0.5
        q = exp(2pi * im * z)
        return digamma(1-z) - pi * im * (q+1)/(q-1)
    end
    digamma(z)
end

function mellin_polynomial_tau(x, tau; angle=3pi/4)
    0 < x < 1 || error("interior x required for the Mellin reference")
    direction = cis(angle)
    function integrand(t)
        N = 1.5 + t * direction
        initial_moment = 6 / prod(N + j for j in 1:4)
        PN = (4/3) * (1.5 - 2 * (contour_digamma(N+1) + Float64(Base.MathConstants.eulergamma)))
        imag(direction * exp(-N * log(x) + tau * PN) * initial_moment) / pi
    end
    quadgk(integrand, 0, 1, 5, 20, Inf; rtol=1e-11, atol=1e-13)[1]
end

function mellin_polynomial(x, Q2, Q02, alpha0; nf=3, angle=3pi/4)
    beta0 = 11 - 2nf / 3
    tau = (2/beta0) * log1p(alpha0 * beta0 * log(Q2/Q02)/(4pi))
    mellin_polynomial_tau(x, tau; angle=angle)
end

function max_scaled(actual, reference)
    scale = max.(abs.(reference), 0.01maximum(abs.(reference)))
    maximum(abs.(actual-reference) ./ max.(scale, 1e-15))
end

end
