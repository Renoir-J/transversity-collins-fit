"Stored inward derivative references at fixed rows; no independent evolution claim."
module HistoricalSIAIntegralReference
using ForwardDiff
export inward_stencil, mathematical_controls, derivative_reference, rows, matrix_rows

rows(A)=collect.(eachrow(A))
matrix_rows(A)=reduce(vcat,permutedims.(A))
scaled(A,B,floor)=maximum(abs.(Float64.(A).-Float64.(B))./max.(abs.(Float64.(B)),floor))
stencil_samples(samples,h)=(-25samples[1]+48samples[2]-36samples[3]+16samples[4]-3samples[5])/(12h)
inward_stencil(f,theta,direction,h)=stencil_samples([f(theta+k*h*direction) for k in 0:4],h)

function mathematical_controls(policy)
    setprecision(policy["reference_bits"]) do
        x=big"0.01";theta=[x];direction=[one(x)];records=Dict{String,Any}[]
        for degree in 0:5,h in parse.(BigFloat,string.(policy["steps"]))
            polynomial=t->t[1]^degree
            samples=[polynomial(theta+k*h*direction) for k in 0:4]
            actual=stencil_samples(samples,h);exact=degree==0 ? zero(x) : degree*x^(degree-1)
            bias=degree==5 ? -24h^4 : zero(h)
            push!(records,Dict("function"=>"monomial","degree"=>degree,"x"=>string(x),"step"=>string(h),
                "samples"=>string.(samples),"estimate"=>string(actual),"analytic_derivative"=>string(exact),
                "known_stencil_bias"=>string(bias),"passed"=>abs(actual-exact-bias)<big"1e-60"))
        end
        exponential=t->exp(7t[1]);errors=BigFloat[]
        for h in parse.(BigFloat,string.(policy["steps"]))
            samples=[exponential(theta+k*h*direction) for k in 0:4]
            actual=stencil_samples(samples,h);exact=7exp(7x);error=actual-exact;push!(errors,error)
            push!(records,Dict("function"=>"exp7","x"=>string(x),"step"=>string(h),"samples"=>string.(samples),
                "estimate"=>string(actual),"analytic_derivative"=>string(exact),"signed_error"=>string(error),
                "passed"=>error<0 && abs(error)/max(abs(exact),BigFloat(policy["relative_floor"]))<=BigFloat(policy["derivative_max_scaled"])))
        end
        Dict("scope"=>"explicitly_synthetic_inward_formula_control","records"=>records,
            "exponential_halving_ratio"=>string(errors[1]/errors[2]),
            "fourth_order_halving_verified"=>abs(errors[1]/errors[2]-16)<big"0.16",
            "passed"=>all(r["passed"] for r in records) && abs(errors[1]/errors[2]-16)<big"0.16")
    end
end

function derivative_reference(candidate,reference,theta,policy;progress=()->nothing)
    direction=Float64.(policy["direction"]);floor=policy["relative_floor"]
    hvps(f,t,d)=reduce(vcat,[permutedims(ForwardDiff.hessian(x->f(x)[i],t)*d) for i in eachindex(f(t))])
    J=ForwardDiff.jacobian(candidate,theta);H=hvps(candidate,theta,direction)
    RJ=ForwardDiff.jacobian(reference,theta);RH=hvps(reference,theta,direction)
    records=Dict{String,Any}[]
    for step in policy["steps"]
        record=setprecision(policy["reference_bits"]) do
            p=BigFloat.(theta);d=BigFloat.(direction);h=parse(BigFloat,string(step))
            columns=Any[]
            for j in eachindex(p)
                e=zeros(BigFloat,length(p));e[j]=1
                push!(columns,[reference(p+k*h*e) for k in 0:4]);progress()
            end
            jacobian=t->ForwardDiff.jacobian(reference,t)
            hvp_samples=[jacobian(p+k*h*d) for k in 0:4];progress()
            fdJ=reduce(hcat,[stencil_samples(c,h) for c in columns]);fdH=stencil_samples(hvp_samples,h)
            Dict("step"=>step,"step_decimal"=>string(h),
                "gradient_samples_decimal"=>[[string.(v) for v in c] for c in columns],
                "HVP_jacobian_samples_decimal"=>[rows(string.(m)) for m in hvp_samples],
                "gradient_decimal"=>rows(string.(fdJ)),"HVP_decimal"=>rows(string.(fdH)),
                "gradient"=>rows(Float64.(fdJ)),"HVP"=>rows(Float64.(fdH)))
        end
        record["gradient_error"]=scaled(J,matrix_rows(record["gradient"]),floor)
        record["HVP_error"]=scaled(H,matrix_rows(record["HVP"]),floor)
        push!(records,record)
    end
    Dict("scope"=>"three_source_selected_geometries_inward_all13_reference_not_production_certificate",
        "theta"=>theta,"direction"=>direction,"gradient"=>rows(J),"HVP"=>rows(H),
        "reference_gradient"=>rows(RJ),"reference_HVP"=>rows(RH),
        "independent_gradient_error"=>scaled(J,RJ,floor),"independent_HVP_error"=>scaled(H,RH,floor),
        "first_six_gradient_columns_zero"=>all(iszero,J[:,1:6]) && all(iszero,RJ[:,1:6]),
        "first_six_HVP_columns_zero"=>all(iszero,H[:,1:6]) && all(iszero,RH[:,1:6]),
        "finite_differences"=>records,
        "adjacent_gradient_error"=>scaled(matrix_rows(records[2]["gradient"]),matrix_rows(records[1]["gradient"]),floor),
        "adjacent_HVP_error"=>scaled(matrix_rows(records[2]["HVP"]),matrix_rows(records[1]["HVP"]),floor))
end
end
