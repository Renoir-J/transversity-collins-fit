using CollinsFit, JSON

"Explicitly synthetic analytic function; neither source observations nor TMD predictions."
function synthetic_sidis_job_reference_records(rate)
    parent=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT);policy=parent["derivative_policy"]
    setprecision(256) do
        theta=BigFloat.(Float64[.1+.01i for i in 1:13]);origin=theta[13]
        direction=BigFloat.(policy["direction"])
        value(t)=sum(t[i]^3 for i in 1:12)+exp(BigFloat(rate)*(t[13]-origin))+t[1]*t[2]
        function gradient(t)
            g=BigFloat[3t[i]^2 for i in 1:12]
            push!(g,BigFloat(rate)*exp(BigFloat(rate)*(t[13]-origin)))
            g[1]+=t[2];g[2]+=t[1];g
        end
        J=gradient(theta);H=BigFloat[6theta[i]*direction[i] for i in 1:12]
        push!(H,BigFloat(rate)^2*direction[13]);H[1]+=direction[2];H[2]+=direction[1]
        rows=Dict{String,Any}[]
        for level in 1:4
            row=sidis_job_reference_step(value,gradient,theta,direction,level,policy;
                base_value=value(theta),base_jacobian=J,reference_HVP=H,
                candidate_gradient=Float64.(J),candidate_HVP=Float64.(H),
                previous=isempty(rows) ? nothing : last(rows))
            push!(rows,row)
        end
        errors=Dict("gradient"=>string(CollinsFit.sidis_job_scaled(Float64.(J),J,big"1e-8")),
            "HVP"=>string(CollinsFit.sidis_job_scaled(Float64.(H),H,big"1e-8")))
        (;rows,errors,J,H,theta,direction,centre_value=value(theta))
    end
end
