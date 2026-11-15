module PaperMoments

using CollinsFit, LinearAlgebra, QuadGK
export moment_weights, reference_integral, direct_initial, rounding_vertices,
       compare_check, scaled_error

"Exact cell integrals of a linear interpolant in x, including partial boundary cells."
function moment_weights(grid, lower, upper)
    x=CollinsFit.collinear_grid(grid)
    isfinite(lower) && isfinite(upper) && first(x)<=lower<upper<=last(x) ||
        error("finite ordered moment bounds inside the grid required; no extrapolation")
    weights=zeros(length(x))
    for j in 1:length(x)-1
        left,right=max(lower,x[j]),min(upper,x[j+1])
        left<right || continue
        scale=(right-left)/(2(x[j+1]-x[j]))
        weights[j]+=scale*((x[j+1]-left)+(x[j+1]-right))
        weights[j+1]+=scale*((left-x[j])+(right-x[j]))
    end
    weights
end

"Independent log-x quadrature of each clipped interpolated cell, including dx=x*dlog(x)."
function reference_integral(grid,values,lower,upper;rtol=1e-12,atol=1e-15)
    x=CollinsFit.collinear_grid(grid)
    length(x)==length(values) && all(isfinite,values) || error("finite reference values required")
    isfinite(lower) && isfinite(upper) && first(x)<=lower<upper<=last(x) || error("invalid reference bounds")
    result=0.;estimate=0.
    for j in 1:length(x)-1
        left,right=max(lower,x[j]),min(upper,x[j+1])
        left<right || continue
        integrand(y)=begin
            v=exp(y)
            (values[j]+(values[j+1]-values[j])*(v-x[j])/(x[j+1]-x[j]))*v
        end
        value,error=quadgk(integrand,log(left),log(right);rtol=rtol,atol=atol)
        result+=value;estimate+=error
    end
    (value=result,estimated_error=estimate)
end

"Independent direct powers in BigFloat, not the production log-shape routine."
function direct_initial(grid,soffer,theta;bits=256)
    size(soffer)==(2,length(grid)) || error("u,d Soffer rows required")
    all(isfinite,soffer) && all(>=(0),soffer) || error("nonnegative finite Soffer inputs required")
    length(theta)==13 || error("13 physical initial parameters required")
    setprecision(bits) do
        out=zeros(length(grid),2)
        for f in 1:2
            n,a,b=BigFloat(theta[f]),BigFloat(theta[f+2]),BigFloat(theta[f+4])
            normalization=(a+b)^(a+b)/(a^a*b^b)
            for (j,value) in enumerate(grid)
                value==1 && continue
                x=BigFloat(value)
                out[j,f]=Float64(n*x^a*(1-x)^b*normalization*BigFloat(soffer[f,j]))
            end
        end
        out
    end
end

"All 64 printing-box vertices; a sampled sensitivity diagnostic, not interval extrema."
function rounding_vertices(theta,bounds)
    length(theta)==13 && length(bounds)==6 || error("six printed transversity intervals required")
    all(b -> length(b)==2 && all(isfinite,b) && b[1]<=b[2],bounds) || error("invalid printing intervals")
    [begin
        t=copy(theta)
        for j in 1:6
            t[j]=bounds[j][1+Int((mask>>(j-1))&1)]
        end
        t
    end for mask in 0:63]
end

scaled_error(a,b,floor)=maximum(abs.(a.-b)./max.(abs.(b),floor))

"Retain actual/reference operands so a status string cannot stand in for a check."
function compare_check(name,actual,reference,tolerance;metric=:scaled,floor=1e-8)
    a,b=Float64.(vec(collect(actual))),Float64.(vec(collect(reference)))
    length(a)==length(b)>0 || error("nonempty matching check operands required")
    isfinite(tolerance) && tolerance>=0 && isfinite(floor) && floor>0 || error("invalid fixed tolerance")
    difference=metric==:scaled ? scaled_error(a,b,floor) : metric==:absolute ? maximum(abs.(a-b)) :
        Base.error("choose scaled or absolute check")
    Dict("name"=>name,"actual"=>a,"reference"=>b,"error"=>difference,"tolerance"=>tolerance,
         "metric"=>string(metric),"relative_floor"=>floor,
         "passed"=>all(isfinite,a) && all(isfinite,b) && isfinite(difference) && difference<=tolerance)
end

end
