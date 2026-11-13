# Independent scalar initial/charge/CSS assembly and adaptive P integrals.
# The fixed native rows/source samples are shared evidence, not a new evolution solution.
function prepare_bin_reference(nodes,weights,intervals,factors)
    refs=[prepare_generic_reference(node.context;momentum_interval=interval)
          for (node,interval) in zip(nodes,intervals)]
    sidis=first(nodes).context.leptonic.process==:SIDIS
    if sidis
        denominator=sum(weights[i]*refs[i].denominators[1] for i in eachindex(refs))
        return t->sum(weights[i]*factors[i]*refs[i].polarized_integrals(t)[1] for i in eachindex(refs))/denominator
    end
    DU=sum(weights[i,1]*(refs[i].denominators[2]+refs[i].denominators[3]) for i in eachindex(refs))
    DL=sum(weights[i,2]*(refs[i].denominators[1]+refs[i].denominators[4]) for i in eachindex(refs))
    function prediction(t)
        rows=[ref.polarized_integrals(t) for ref in refs]
        U=sum(weights[i,1]*factors[i]*(rows[i][2]+rows[i][3]) for i in eachindex(refs))
        L=sum(weights[i,2]*factors[i]*(rows[i][1]+rows[i][4]) for i in eachindex(refs))
        U/DU-(first(nodes).kind==:UL ? L/DL : (U+L)/(DU+DL))
    end
    prediction
end

"Full 13-column Jacobian and directional HVP compared with inward 256-bit references."
function bin_derivative_checks(bin,reference,theta,policy,checks,label)
    actual(t)=[predict_bin(bin,t)]
    expected(t)=[reference(t)]
    direction=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,.001]
    J=ForwardDiff.jacobian(actual,theta)
    HVP=reshape(ForwardDiff.jacobian(s->vec(ForwardDiff.jacobian(actual,theta+s[1]*direction)),[0.])[:,1],size(J))
    generic_J=setprecision(policy["finite_difference_reference_bits"]) do
        Float64.(ForwardDiff.jacobian(expected,BigFloat.(theta)))
    end
    observable_check!(checks,label*" independent full 13-column gradient",observable_scaled(J,generic_J,policy),policy["derivative_max_scaled"])
    results=[]
    for step in policy["finite_difference_steps"]
        fd,fdh=setprecision(policy["finite_difference_reference_bits"]) do
            t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
            gradient(x)=ForwardDiff.jacobian(expected,x)
            (Float64.((-3expected(t)+4expected(t+h*d)-expected(t+2h*d))/(2h)),
             Float64.((-3gradient(t)+4gradient(t+h*d)-gradient(t+2h*d))/(2h)))
        end
        e1=observable_scaled(J*direction,fd,policy);e2=observable_scaled(HVP,fdh,policy)
        observable_check!(checks,label*" inward full-13 JVP step $step",e1,policy["derivative_max_scaled"])
        observable_check!(checks,label*" inward full-13 HVP step $step",e2,policy["derivative_max_scaled"])
        push!(results,Dict("label"=>label,"step"=>step,"parameters"=>theta,"direction"=>direction,
            "actual_gradient"=>J,"reference_gradient"=>generic_J,"actual_JVP"=>J*direction,"reference_JVP"=>fd,
            "actual_HVP"=>HVP,"reference_HVP"=>fdh,"JVP_error"=>e1,"HVP_error"=>e2,
            "reference_bits"=>policy["finite_difference_reference_bits"],
            "reference_scope"=>"independent scalar initial/charge/CSS assembly, generic fixed-row products and adaptive P integration; shared native rows"))
    end
    results
end
