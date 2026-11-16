module HistoricalSIDISSourceMath
using CollinsFit, JSON, QuadGK, ForwardDiff

# Explicitly synthetic count/stencil controls under the unchanged scientific
# contract. These functions do not prepare native rows or physical predictions.
function cases(ctx)
    [Dict("id"=>p["id"],"publication"=>p["publication"],"reported_y"=>p["kinematics"]["y"],
        "view_id"=>p["dataset_id"]) for p in ctx.request.data["points"]]
end

function populate!(report,ctx)
    contract=ctx.contract;fixture=Dict("cases"=>cases(ctx))
    report["record_kind"]="explicitly_synthetic_count_and_closed_stencil_controls"
    report["physical_predictions_executed"]=false
    report["native_or_provider_calls_executed"]=false
        setprecision(256) do
            policy=contract["numerical_policy"];dp=contract["derivative_policy"]
            floor=BigFloat(string(policy["relative_floor"]));tolerance=BigFloat(string(policy["source_projection_max_scaled"]))
            err(a,b)=abs(a-b)/max(abs(b),floor)
            pi0=BigFloat(pi);limit=2pi0;normal=2pi0^2
            twice(f)=first(quadgk(s->first(quadgk(h->f(h,s),zero(pi0),limit;rtol=big"1e-45",atol=big"1e-60")),zero(pi0),limit;rtol=big"1e-45",atol=big"1e-60"))
            harmonic(phase,h,s)=sin(h+s-(phase=="sin_phi_sum_minus_pi" ? pi0 : zero(pi0)))
            bases=Dict{String,Any}()
            for (phase,sign) in (("sin_phi_sum",1),("sin_phi_sum_minus_pi",-1))
                b=(constant=twice((h,s)->harmonic(phase,h,s)),
                    collins=twice((h,s)->sin(h+s)*harmonic(phase,h,s)),
                    sivers=twice((h,s)->sin(h-s)*harmonic(phase,h,s)))
                @assert err(b.constant,zero(pi0))<=tolerance && err(b.collins,sign*normal)<=tolerance && err(b.sivers,zero(pi0))<=tolerance
                bases[phase]=b
            end
            totals=(constant=twice((h,s)->one(pi0)),collins=twice((h,s)->sin(h+s)),sivers=twice((h,s)->sin(h-s)))
            @assert err(totals.constant,4pi0^2)<=tolerance && err(totals.collins,zero(pi0))<=tolerance && err(totals.sivers,zero(pi0))<=tolerance
            primary=Dict("1006.4221v1"=>("sin_phi_sum","fixed_y_DNN_included"),
                "0802.2160"=>("sin_phi_sum_minus_pi","DNN_removed"),
                "1408.4405"=>("sin_phi_sum_minus_pi","DNN_removed"))
            cases=fixture["cases"];@assert length(cases)==length(Set(c["id"] for c in cases))==132
            rows=Any[];direct=Any[];seen=Set{String}();P=big"0.5";S=big"0.0625"
            for point in cases
                pub=point["publication"];phase,norm=primary[pub];mapping=contract["model_conventions"][pub]
                @assert mapping["phase"]==phase && mapping["normalization"]==norm
                y=BigFloat(string(point["reported_y"]));@assert 0<y<1
                D=(1-y)/(1-y+y^2/2)
                convention=norm=="DNN_removed" ? SIDISAsymmetryConvention(phase=Symbol(mapping["phase"]),normalization=:DNN_removed,source="frozen synthetic control for "*point["view_id"]) :
                    SIDISAsymmetryConvention(phase=Symbol(mapping["phase"]),normalization=:fixed_y_DNN_included,y=point["reported_y"],source="frozen synthetic control for "*point["view_id"])
                metadata=CollinsFit.sidis_convention_metadata(convention);b=bases[phase]
                for R in (big"-0.125",big"0.125")
                    @assert P*(D*abs(R)+abs(S))<1
                    numerator=b.constant+P*D*R*b.collins+P*S*b.sivers
                    total=totals.constant+P*D*R*totals.collins+P*S*totals.sivers
                    reference=norm=="DNN_removed" ? numerator/(normal*P*D) : 2numerator/(total*P)
                    implementation=BigFloat(metadata["phase_factor"]*metadata["DNN_factor"]*Float64(R))
                    error=err(implementation,reference);@assert error<=tolerance
                    wrong_DNN=norm=="DNN_removed" ? implementation*D : implementation/D
                    @assert all(abs(v-reference)>floor for v in (-implementation,wrong_DNN,2implementation))
                    push!(rows,Dict("id"=>point["id"],"publication"=>pub,"reported_y"=>point["reported_y"],
                        "R"=>Float64(R),"phase"=>phase,"normalization"=>norm,"implementation"=>Float64(implementation),
                        "reference_decimal"=>string(reference),"error"=>Float64(error),"counterfactuals_distinct"=>true,
                        "scope"=>"synthetic_count_coefficient_not_physical_prediction"))
                    if !(pub in seen)
                        N=(h,s)->1+P*D*R*sin(h+s)+P*S*sin(h-s)
                        q=twice((h,s)->N(h,s)*harmonic(phase,h,s));qtotal=twice(N)
                        @assert err(q,numerator)<=tolerance && err(qtotal,total)<=tolerance
                        push!(direct,Dict("publication"=>pub,"R"=>Float64(R),"direct_numerator_decimal"=>string(q),"basis_numerator_decimal"=>string(numerator)))
                    end
                end
                push!(seen,pub)
            end
            @assert length(rows)==264 && length(direct)==6
            report["source_factor_records"]=rows;report["direct_full_count_projections"]=direct
            report["independent_two_angle_basis_integrals"]=6;report["independent_density_normalization_integrals"]=3
            report["basis_linearity_reused_for_all_source_y"]=true
            target=Any[];Up=big"2";Un=big"1";Rp=big"0.125";Rn=big"-0.0625";y=big"0.4";D=(1-y)/(1-y+y^2/2)
            for (name,w) in (("proton",(big"1",big"0")),("deuteron",(big"0.5",big"0.5")))
                @assert contract["prediction"]["target_weights"][name]["unpolarized"]==Float64.([w...])
                @assert contract["prediction"]["target_weights"][name]["polarized"]==Float64.([w...])
                U=w[1]*Up+w[2]*Un;F=w[1]*Up*Rp+w[2]*Un*Rn;ratio=F/U
                N=(h,s)->w[1]*Up*(1+P*D*Rp*sin(h+s))+w[2]*Un*(1+P*D*Rn*sin(h+s))
                reference=twice((h,s)->N(h,s)*sin(h+s))/(normal*P*D*U)
                @assert err(ratio,reference)<=tolerance
                name=="deuteron" && @assert abs(ratio-(Rp+Rn)/2)>floor
                push!(target,Dict("target"=>name,"UU"=>Float64(U),"FUT"=>Float64(F),"ratio"=>Float64(ratio),"reference_decimal"=>string(reference)))
            end
            report["synthetic_target_records"]=target
            inward(f,x,h)=(-25f(x)+48f(x+h)-36f(x+2h)+16f(x+3h)-3f(x+4h))/(12h)
            control=dp["closed_formula_control"];x=BigFloat(string(control["x"]));steps=BigFloat.(string.(dp["inherited_initial_steps"]))
            records=Any[];exp_errors=BigFloat[]
            for degree in control["degrees"],h in steps
                f=t->t^degree;estimate=inward(f,x,h);exact=degree==0 ? zero(x) : degree*x^(degree-1);bias=degree==5 ? -24h^4 : zero(x)
                @assert abs(estimate-exact-bias)<BigFloat(string(control["monomial_residual_abs_max"]))
                push!(records,Dict("function"=>"monomial","degree"=>degree,"step_decimal"=>string(h),"estimate_decimal"=>string(estimate),"analytic_decimal"=>string(exact),"bias_decimal"=>string(bias),"passed"=>true))
            end
            for h in steps
                f=t->exp(7t);estimate=inward(f,x,h);exact=7exp(7x);error=estimate-exact
                @assert error<0 && err(estimate,exact)<=BigFloat(string(control["gradient_max_scaled"]))
                push!(exp_errors,error);push!(records,Dict("function"=>"exp7","step_decimal"=>string(h),"estimate_decimal"=>string(estimate),"analytic_decimal"=>string(exact),"signed_error_decimal"=>string(error),"passed"=>true))
            end
            expected=(steps[1]/steps[2])^control["expected_convergence_order"]
            @assert abs(abs(exp_errors[1]/exp_errors[2])/expected-1)<BigFloat(string(control["step_ratio_relative_tolerance"]))
            @assert length(records)==14
            report["closed_formula_records"]=records;report["exponential_step_error_ratio_decimal"]=string(exp_errors[1]/exp_errors[2])
            w=control["width_fixture"];t=BigFloat(string(w["point"]));f=z->(big"0.042"-z)^(-2);jac=z->2(big"0.042"-z)^(-3)
            G=jac(t);H=6(big"0.042"-t)^(-4)
            @assert err(ForwardDiff.derivative(f,t),G)<=tolerance && err(ForwardDiff.derivative(jac,t),H)<=tolerance
            levels=Any[];accepted=false
            for (i,step) in enumerate(w["steps"])
                h=BigFloat(string(step));g=inward(f,t,h);v=inward(jac,t,h);ge=err(g,G);he=err(v,H)
                ag=isempty(levels) ? nothing : err(g,parse(BigFloat,last(levels)["gradient_decimal"]))
                ah=isempty(levels) ? nothing : err(v,parse(BigFloat,last(levels)["Hessian_decimal"]))
                accepted=i>=w["minimum_refinement_level"] && max(ge,he,ag,ah)<=BigFloat(string(w["gradient_max_scaled"]))
                push!(levels,Dict("step"=>step,"gradient_decimal"=>string(g),"Hessian_decimal"=>string(v),
                    "gradient_error"=>Float64(ge),"Hessian_error"=>Float64(he),"adjacent_gradient_error"=>ag===nothing ? nothing : Float64(ag),
                    "adjacent_Hessian_error"=>ah===nothing ? nothing : Float64(ah),"accepted"=>accepted))
                accepted && break
            end
            @assert accepted
            report["synthetic_width_refinements"]=levels;report["synthetic_width_fixture_scope"]=w["scope"]
            report["status"]="SYNTHETIC_SOURCE_MATH_RECORDED_NOT_DIAGNOSTIC_GATE"
        end
    report
end

function audit(report,ctx)
    c=ctx.contract;input=Dict("cases"=>cases(ctx));r=report
    r["status"]=="SYNTHETIC_SOURCE_MATH_RECORDED_NOT_DIAGNOSTIC_GATE" &&
        r["record_kind"]=="explicitly_synthetic_count_and_closed_stencil_controls" &&
        r["physical_predictions_executed"]===r["native_or_provider_calls_executed"]===false || error("synthetic source controls cannot admit physics")
    require(b,message)=b || error(message)
    setprecision(256) do
        p=c["numerical_policy"];dp=c["derivative_policy"];floor=BigFloat(string(p["relative_floor"]));tolerance=BigFloat(string(p["source_projection_max_scaled"]))
        err(a,b)=abs(a-b)/max(abs(b),floor)
        rows=r["source_factor_records"];require(length(rows)==264,"incomplete source-factor record")
        primary=Dict("1006.4221v1"=>(1,false),"0802.2160"=>(-1,true),"1408.4405"=>(-1,true))
        for (i,(point,R)) in enumerate([(point,R) for point in input["cases"] for R in (big"-0.125",big"0.125")])
            row=rows[i];sign,removed=primary[point["publication"]];y=BigFloat(string(point["reported_y"]));D=(1-y)/(1-y+y^2/2)
            expected=sign*R*(removed ? 1 : D);ref=parse(BigFloat,row["reference_decimal"]);actual=BigFloat(row["implementation"])
            require(row["id"]==point["id"] && row["publication"]==point["publication"] && row["reported_y"]==point["reported_y"] && row["R"]==Float64(R),"changed source case")
            require(err(ref,expected)<=tolerance && row["error"]==Float64(err(actual,ref)) && err(actual,ref)<=tolerance,"changed source-factor arithmetic")
            wrong=removed ? actual*D : actual/D
            require(row["counterfactuals_distinct"]===true && all(abs(v-ref)>floor for v in (-actual,wrong,2actual)),"changed source counterfactual")
        end
        direct=r["direct_full_count_projections"];require(length(direct)==6,"incomplete direct projection list")
        for row in direct
            point=first(filter(p->p["publication"]==row["publication"],input["cases"]));sign,_=primary[row["publication"]]
            y=BigFloat(string(point["reported_y"]));D=(1-y)/(1-y+y^2/2);R=BigFloat(row["R"])
            exact=sign*big"0.5"*D*R*2BigFloat(pi)^2
            require(err(parse(BigFloat,row["direct_numerator_decimal"]),exact)<=tolerance && err(parse(BigFloat,row["basis_numerator_decimal"]),exact)<=tolerance,"changed direct projection operands")
        end
        target=r["synthetic_target_records"];require(length(target)==2,"incomplete target fixture")
        for row in target
            U,F=row["target"]=="proton" ? (big"2",big"0.25") : (big"1.5",big"0.09375")
            require(row["UU"]==Float64(U) && row["FUT"]==Float64(F) && row["ratio"]==Float64(F/U) && err(parse(BigFloat,row["reference_decimal"]),F/U)<=tolerance,"changed target sum")
        end
        inward(f,x,h)=(-25f(x)+48f(x+h)-36f(x+2h)+16f(x+3h)-3f(x+4h))/(12h)
        control=dp["closed_formula_control"];x=BigFloat(string(control["x"]));hs=BigFloat.(string.(dp["inherited_initial_steps"]));records=r["closed_formula_records"]
        require(length(records)==14,"incomplete closed formula list")
        for (i,(degree,h)) in enumerate([(degree,h) for degree in control["degrees"] for h in hs])
            row=records[i];f=t->t^degree;estimate=inward(f,x,h);exact=degree==0 ? zero(x) : degree*x^(degree-1);bias=degree==5 ? -24h^4 : zero(x)
            require(row["degree"]==degree && row["step_decimal"]==string(h) && parse(BigFloat,row["estimate_decimal"])==estimate && parse(BigFloat,row["analytic_decimal"])==exact && parse(BigFloat,row["bias_decimal"])==bias && row["passed"]===true,"changed monomial control")
            require(abs(estimate-exact-bias)<BigFloat(string(control["monomial_residual_abs_max"])),"failed monomial formula")
        end
        es=BigFloat[]
        for (i,h) in enumerate(hs)
            row=records[12+i];f=t->exp(7t);estimate=inward(f,x,h);exact=7exp(7x);e=estimate-exact
            require(row["function"]=="exp7" && row["step_decimal"]==string(h) && parse(BigFloat,row["estimate_decimal"])==estimate && parse(BigFloat,row["analytic_decimal"])==exact && parse(BigFloat,row["signed_error_decimal"])==e && row["passed"]===true,"changed exponential control")
            require(e<0 && err(estimate,exact)<=BigFloat(string(control["gradient_max_scaled"])),"failed exponential formula");push!(es,e)
        end
        require(parse(BigFloat,r["exponential_step_error_ratio_decimal"])==es[1]/es[2] && abs(abs(es[1]/es[2])/(hs[1]/hs[2])^4-1)<BigFloat(string(control["step_ratio_relative_tolerance"])),"changed formula convergence ratio")
        w=control["width_fixture"];t=BigFloat(string(w["point"]));f=z->(big"0.042"-z)^(-2);jac=z->2(big"0.042"-z)^(-3);G=jac(t);H=6(big"0.042"-t)^(-4)
        levels=r["synthetic_width_refinements"];require(w["minimum_refinement_level"]<=length(levels)<=length(w["steps"]),"incomplete width refinement")
        previous=nothing
        for (i,row) in enumerate(levels)
            h=BigFloat(string(w["steps"][i]));g=inward(f,t,h);v=inward(jac,t,h);ge=err(g,G);he=err(v,H)
            ag=previous===nothing ? nothing : err(g,previous[1]);ah=previous===nothing ? nothing : err(v,previous[2])
            accepted=i>=w["minimum_refinement_level"] && max(ge,he,ag,ah)<=BigFloat(string(w["gradient_max_scaled"]))
            require(row["step"]==w["steps"][i] && parse(BigFloat,row["gradient_decimal"])==g && parse(BigFloat,row["Hessian_decimal"])==v && row["gradient_error"]==Float64(ge) && row["Hessian_error"]==Float64(he),"changed width derivative arithmetic")
            require(row["adjacent_gradient_error"]==(ag===nothing ? nothing : Float64(ag)) && row["adjacent_Hessian_error"]==(ah===nothing ? nothing : Float64(ah)) && row["accepted"]===accepted,"changed width convergence decision")
            require(!accepted || i==length(levels),"continued after width acceptance");previous=(g,v)
        end
        require(last(levels)["accepted"]===true,"unaccepted final width reference")
    end
    Dict("source_factor_errors_recomputed"=>264,"closed_functions_reevaluated"=>14,
        "target_sum_records_recomputed"=>2,"width_refinement_decisions_recomputed"=>length(r["synthetic_width_refinements"]),
        "adaptive_Fourier_quadrature_reexecuted"=>false,"native_or_physical_prediction_reexecuted"=>false)
end
end
