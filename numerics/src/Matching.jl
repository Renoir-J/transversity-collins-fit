"Fixed Float64 quadrature settings; fitted input values retain generic Julia arithmetic."
const GAUSS_RULE_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const GAUSS_RULE_LOCK = ReentrantLock()

function gauss_rule(order)
    order isa Integer && !(order isa Bool) && 2 <= order <= 512 || error("Gauss order must be 2..512")
    lock(GAUSS_RULE_LOCK) do
        nodes,weights = get!(GAUSS_RULE_CACHE,Int(order)) do
            # Prepare accurate symmetric roots/weights once, then round to the
            # declared Float64 numerical path. Eigenvector weights can retain
            # larger errors in strongly cancelling Fourier integrals.
            setprecision(192) do
                n = Int(order)
                roots,ws = zeros(n),zeros(n)
                function polynomial(z)
                    p,previous = one(z),zero(z)
                    for k in 1:n
                        p,previous = ((2k-1)*z*p-(k-1)*previous)/k,p
                    end
                    p,n*(z*p-previous)/(z^2-1)
                end
                for k in 1:cld(n,2)
                    z = cos(BigFloat(pi)*(k-big"0.25")/(n+big"0.5"))
                    converged = false
                    for _ in 1:20
                        p,derivative = polynomial(z)
                        change = p/derivative
                        z -= change
                        if abs(change) <= 16eps(one(z))
                            converged = true
                            break
                        end
                    end
                    converged || error("Legendre root did not converge")
                    _,derivative = polynomial(z)
                    w = 2/((1-z^2)*derivative^2)
                    roots[k],roots[n+1-k] = -Float64(z),Float64(z)
                    ws[k]=ws[n+1-k]=Float64(w)
                end
                roots,ws
            end
        end
        copy(nodes),copy(weights)
    end
end

function matching_channel(process,channel)
    ((process == :SIDIS && channel in (:pdf,:transversity,:ff,:collins)) ||
     (process == :SIA && channel in (:ff,:collins))) || error("unsupported process/matching channel")
end

function matching_delta(process,channel)
    matching_channel(process,channel)
    process == :SIDIS ? -8/3 : (pi^2-8)/3
end

function matching_regular(channel,v)
    0 < v <= 1 || error("matching fraction must be in (0,1]")
    L = v == 1 ? -one(v) : (v > .5 ? log1p(v-1) : log(v))/(1-v)
    if channel == :pdf
        return ((2/3)*(1-v),.5v*(1-v))
    elseif channel == :ff
        return ((4/3)*((1-v)/2+(1+v^2)*L),
                (4/3)*(v/2+(1+(1-v)^2)*log(v)/v))
    elseif channel == :collins
        return ((8/3)*v*L,zero(v))
    elseif channel == :transversity
        return (zero(v),zero(v))
    end
    error("unsupported matching channel")
end

struct MatchingRule
    process::Symbol
    channel::Symbol
    x::Float64
    alpha::Float64
    coordinates::Vector{Float64} # x first, then convolution queries y
    quark_weights::Vector{Float64}
    gluon_weights::Vector{Float64}
end

"One unweighted matching leg; quadrature is in log(y), with explicit knot splitting."
function matching_rule(process,channel,x,alpha;order=64,knots=Float64[])
    matching_channel(process,channel)
    x,a = css_z(x),fixed_scalar(alpha,"matching alpha_s")
    a >= 0 || error("matching alpha_s must be nonnegative")
    nodes,weights = gauss_rule(order)
    grid = [fixed_scalar(v,"matching knot") for v in knots]
    all(v -> 0 < v <= 1,grid) && all(diff(grid).>0) || error("matching knots must increase in (0,1]")
    query = [x]
    wq,wg = [1+(a/pi)*matching_delta(process,channel)],[0.]
    if a > 0 && channel != :transversity && x < 1
        edges = vcat(x,filter(y -> x < y < 1,grid),1.)
        for i in 1:length(edges)-1
            lo,hi = log(edges[i]),log(edges[i+1])
            half = (hi-lo)/2
            for k in eachindex(nodes)
                y = exp((lo+hi)/2+half*nodes[k])
                rq,rg = matching_regular(channel,x/y)
                push!(query,y)
                push!(wq,(a/pi)*half*weights[k]*rq)
                push!(wg,(a/pi)*half*weights[k]*rg)
            end
        end
    end
    MatchingRule(process,channel,x,a,query,wq,wg)
end

"Apply fixed matching weights to generic (possibly AD-valued) parton inputs."
function match_leg(rule::MatchingRule,quark;gluon=nothing)
    length(quark)==length(rule.coordinates) && all(isfinite,quark) || error("invalid matching quark values")
    value = dot(rule.quark_weights,quark)
    if isnothing(gluon)
        all(iszero,rule.gluon_weights) || error("matching requires explicit gluon values")
    else
        rule.channel in (:pdf,:ff) || error("no gluon matching for transversity/Collins")
        length(gluon)==length(quark) && all(isfinite,gluon) || error("invalid matching gluon values")
        value += dot(rule.gluon_weights,gluon)
    end
    isfinite(value) || error("nonfinite matching result")
    value
end
