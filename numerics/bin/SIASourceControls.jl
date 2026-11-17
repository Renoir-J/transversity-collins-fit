"Synthetic angle/stencil controls; no physical SIA integral or fitting admission."
module SIASourceControls
using JSON, SHA
include("SIAGeometrySelection.jl")
const S=SIAGeometrySelection
const SCOPE="explicitly_synthetic_sia_source_and_stencil_controls"
const COMPUTED="SYNTHETIC_SIA_SOURCE_CONTROLS_COMPUTED_NOT_DIAGNOSTIC_GATE"
const REVIEWED="SYNTHETIC_SIA_SOURCE_CONTROL_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE"
const FLAGS=(S.FLAGS...,"accepted_fitting_gate","accepted_prediction_diagnostic","numerical_dispatch_ready")
flags()=Dict{String,Any}(key=>false for key in FLAGS)
require(value,message)=value || error(message)
finite(x)=x isa Real && !(x isa Bool) && isfinite(x)
relative(a,b,floor)=abs(a-b)/max(abs(b),floor)
check(name,error,tolerance)=Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
    "passed"=>isfinite(error) && error<=tolerance)
boolean(name,value)=check(name,value ? 0. : 1.,0.)

function produce(data,parent,binding,projection,formula)
    S.validate_parent(parent);S.geometry_groups(data)
    policy=parent["numerical_policy"];floor=policy["relative_floor"]
    projections=projection(data["points"],policy)
    math=formula(parent["representative_derivative_policy"])
    checks=Any[boolean("inward stencil closed-function controls",math["passed"])]
    require(length(projections)==144,"all72 points and both signs are required")
    for (i,row) in enumerate(projections)
        push!(checks,check("source RF0 independent linearized Fourier projection $i",
            relative(row["projected"],row["reference"],floor),policy["source_projection_max_scaled"]))
        push!(checks,check("source angular convention application $i",row["error"],policy["source_projection_max_scaled"]))
        push!(checks,boolean("source fixture positive/distinct $i",row["positive_counts_proved_by_bound"] && row["all_counterfactuals_differ_from_correct"]))
    end
    passed=all(row["passed"] for row in checks)
    merge(flags(),Dict{String,Any}("schema_version"=>1,"scope"=>SCOPE,"theory_scope"=>"synthetic_closure",
        "status"=>passed ? COMPUTED : "SYNTHETIC_SIA_SOURCE_CONTROLS_FAILED",
        "source_controls_success"=>passed,"source_context"=>deepcopy(binding),
        "request_points"=>deepcopy(data["points"]),"source_projection"=>projections,
        "stencil_formula_control"=>math,"checks"=>checks,
        "physical_SIA_evolution_or_integrals_or_AD_executed"=>false,
        "raw_primary_tables_independently_reparsed_in_Julia"=>false,
        "scope_limit"=>"Synthetic Fourier and closed-function controls using source angular scalars; not experimental validation, independent evolution, physical integrals or a production prerequisite receipt by itself."))
end

function audit_stencil(math,dp)
    require(math["scope"]=="explicitly_synthetic_inward_formula_control" && length(math["records"])==14,"changed mathematical control scope")
    setprecision(dp["reference_bits"]) do
        errors=BigFloat[];passes=Bool[]
        grid=vcat([(degree,h) for degree in 0:5 for h in dp["steps"]],[(nothing,h) for h in dp["steps"]])
        for (record,(degree,step)) in zip(math["records"],grid)
            x=big"0.01";h=parse(BigFloat,string(step))
            require(record["x"]==string(x) && record["step"]==string(h),"changed mathematical control grid")
            samples=degree===nothing ? [exp(7*(x+k*h)) for k in 0:4] : [(x+k*h)^degree for k in 0:4]
            estimate=(-25samples[1]+48samples[2]-36samples[3]+16samples[4]-3samples[5])/(12h)
            require(record["samples"]==string.(samples) && record["estimate"]==string(estimate),"changed closed-function operands")
            if degree===nothing
                exact=7exp(7x);error=estimate-exact;push!(errors,error)
                passed=error<0 && abs(error)/max(abs(exact),BigFloat(dp["relative_floor"]))<=BigFloat(dp["derivative_max_scaled"])
                require(record["function"]=="exp7" && record["signed_error"]==string(error),"changed exponential control")
            else
                exact=degree==0 ? zero(x) : degree*x^(degree-1);bias=degree==5 ? -24h^4 : zero(h)
                passed=abs(estimate-exact-bias)<big"1e-60"
                require(record["function"]=="monomial" && record["degree"]===degree && record["known_stencil_bias"]==string(bias),"changed polynomial control")
            end
            require(record["analytic_derivative"]==string(exact) && record["passed"]===passed,"changed closed-function decision")
            push!(passes,passed)
        end
        halving=abs(errors[1]/errors[2]-16)<big"0.16"
        require(math["exponential_halving_ratio"]==string(errors[1]/errors[2]) && math["fourth_order_halving_verified"]===halving &&
            math["passed"]===(all(passes) && halving),"changed mathematical-control aggregate")
        all(passes) && halving
    end
end

function audit(data,parent,binding,report)
    S.validate_parent(parent);S.geometry_groups(data)
    require(report["schema_version"]===1 && report["status"]==COMPUTED && report["scope"]==SCOPE &&
        report["theory_scope"]=="synthetic_closure" && report["source_controls_success"]===true,"unaccepted synthetic source control")
    require(all(get(report,key,nothing)===false for key in (FLAGS...,"physical_SIA_evolution_or_integrals_or_AD_executed",
        "raw_primary_tables_independently_reparsed_in_Julia")),"synthetic control cannot authorize physical work")
    require(S.same_json(report["source_context"],binding) && S.same_json(report["request_points"],data["points"]),"changed source context or point fields")
    policy=parent["numerical_policy"];floor=policy["relative_floor"]
    checks=Any[boolean("inward stencil closed-function controls",audit_stencil(report["stencil_formula_control"],parent["representative_derivative_policy"]))]
    projections=report["source_projection"];grid=[(point,sign) for point in data["points"] for sign in (-1,1)]
    require(length(projections)==length(grid)==144,"incomplete source angular grid")
    for (i,(row,(point,sign))) in enumerate(zip(projections,grid))
        f=point["source_angular"]["value"]
        require(row["scope"]=="explicitly_synthetic_source_definition_fixture" && row["point_id"]==point["id"] &&
            row["ratio"]==point["ratio"] && row["angular_factor"]===f && row["sign"]===sign &&
            all(finite,[row[key] for key in ("D_U","D_L","C_U","C_L","linearization_at_epsilon")]) &&
            row["D_U"]==2. && row["D_L"]==1. && row["C_U"]==sign/8 && row["C_L"]==-sign/16 &&
            all(finite,row["epsilon_interval"]) && row["epsilon_interval"]==[-1.,1.] && row["linearization_at_epsilon"]==0.,"changed synthetic linearized source fixture")
        closed=setprecision(policy["source_projection_reference_bits"]) do
            a=BigFloat(f);du=BigFloat(2);dl=BigFloat(1);cu=BigFloat(sign)/8;cl=-BigFloat(sign)/16
            dc,cc=du+dl,cu+cl;dr,cr=point["ratio"]=="UL" ? (dl,cl) : (dc,cc)
            difference=cu/du-cr/dr;value=a*difference
            other=point["ratio"]=="UL" ? a*(cu/du-cc/dc) : a*(cu/du-cl/dl)
            (reference=Float64(value),positive=all(d>a*abs(c) for (d,c) in ((du,cu),(dl,cl),(dc,cc))),
                counterfactuals=Dict("omitted_Ftheta"=>Float64(difference),"repeated_Ftheta"=>Float64(a*value),
                    "factor2"=>Float64(2value),"wrong_category"=>Float64(other),"reversed_sign"=>Float64(-value),
                    "cos_phi0"=>0.,"sin_2phi0"=>0.))
        end
        require(finite(row["reference"]) && row["reference"]==closed.reference && finite(row["projected"]),"changed analytic angular coefficient")
        du,dl,cu,cl=2.,1.,sign/8,-sign/16
        dr,cr=point["ratio"]=="UL" ? (dl,cl) : (du+dl,cu+cl)
        implementation=f*(cu/du-cr/dr)
        require(row["implementation"]==implementation && Set(keys(row["counterfactuals"]))==Set(keys(closed.counterfactuals)),"changed angular implementation or counterfactuals")
        for (name,value) in closed.counterfactuals
            require(finite(row["counterfactuals"][name]) && relative(row["counterfactuals"][name],value,floor)<=policy["observable_algebra_max_scaled"],"changed angular counterfactual")
        end
        distinct=all(abs(value-row["projected"])>floor for value in values(row["counterfactuals"]))
        require(row["positive_counts_proved_by_bound"]===closed.positive && row["all_counterfactuals_differ_from_correct"]===distinct,"changed positivity or distinctness bit")
        error=relative(implementation,row["projected"],floor)
        require(row["error"]==error,"changed source application error")
        push!(checks,check("source RF0 independent linearized Fourier projection $i",relative(row["projected"],row["reference"],floor),policy["source_projection_max_scaled"]))
        push!(checks,check("source angular convention application $i",error,policy["source_projection_max_scaled"]))
        push!(checks,boolean("source fixture positive/distinct $i",closed.positive && distinct))
    end
    require(length(checks)==433 && S.same_json(report["checks"],checks) && all(row["passed"] for row in checks),"missing, changed or failed synthetic check ledger")
    merge(flags(),Dict{String,Any}("schema_version"=>1,"status"=>REVIEWED,"scope"=>SCOPE,"theory_scope"=>"synthetic_closure",
        "source_context"=>deepcopy(binding),"source_controls_success"=>true,"stored_arithmetic_checks_recomputed"=>433,
        "closed_mathematical_functions_reevaluated"=>true,"physical_SIA_evolution_or_integrals_or_AD_executed"=>false,
        "Fourier_quadrature_reexecuted_by_reader"=>false,"raw_primary_tables_independently_reparsed_in_Julia"=>false))
end
end
