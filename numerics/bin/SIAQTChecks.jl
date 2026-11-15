# Fixed-row integration checks only. Julia owns all numerical references.
using CollinsFit, JSON, SHA, ForwardDiff, LinearAlgebra, QuadGK
isdefined(@__MODULE__,:ObservableReference) || include("ObservableChecks.jl")
using .ObservableReference

const QT_CONTRACT=joinpath(@__DIR__,"../../config/sia-qt-window-contract.json")
const QT_FD_CONTRACT=joinpath(@__DIR__,"../../config/sia-qt-fd-reference-contract.json")

function sia_qt_sources(producer)
    merge(observable_source_hashes(producer),Dict(relpath(p,OBS_ROOT)=>bytes2hex(sha256(read(p)))
        for p in (QT_CONTRACT,QT_FD_CONTRACT,@__FILE__)))
end

qt_vector(x::Number)=[x]
qt_vector(x)=vec(collect(x))
function qt_check!(checks,name,actual,reference,tolerance,policy)
    a=Float64.(qt_vector(actual));b=Float64.(qt_vector(reference))
    length(a)==length(b) || error("reference size mismatch: $name")
    scaled=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
    push!(checks,Dict("name"=>name,"actual"=>a,"reference"=>b,"error"=>scaled,
        "relative_floor"=>policy["relative_floor"],"tolerance"=>tolerance,
        "passed"=>isfinite(scaled)&&scaled<=tolerance))
end
qt_assert!(checks,name,condition,policy)=qt_check!(checks,name,condition ? 0. : 1.,0.,0.,policy)
qt_rows(matrix)=[collect(matrix[i,:]) for i in axes(matrix,1)]

"Explicit synthetic response rows; no native evolution/provider is invoked."
function sia_qt_fixture(z1,settings)
    model=KPSYCSS();radial=RadialRule(Float64.(settings["radial_edges"]);order=settings["radial_order"])
    z2=Float64(settings["z2"])
    plan=leptonic_plan(:SIA,model,radial,settings["Q2"],.2,z1;z2=z2)
    grid=sort!(unique(vcat(.05,.1,z1,z2,.75,1.)))
    scales=css_mu_b2.(Ref(model),radial.b)
    times=evolution_time.(Ref(model.coupling),scales,model.Q02)
    alpha=alpha_s.(Ref(model.coupling),scales)
    blocks=CollinsFit.HoppetPointBlock[]
    for query in (z1,z2)
        rows=zeros(length(scales),length(grid))
        rows[:,findfirst(==(query),grid)].=1 .+ .01 .*radial.b
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),rows))
    end
    data=(grid,[z1,z2],scales,model.Q02,times,:KPSY_Hhat3_homogeneous,5,32,:synthetic_fixture,
        blocks,sum(sizeof(block.rows) for block in blocks),0.)
    operator=HoppetMatchingPlan(data...,:SIA,alpha,
        CollinsFit.hoppet_point_snapshot(data...;process=:SIA,alpha=alpha))
    samples=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,grid,
        [(1+.3i)*x*(1-x)^3 for i in 1:4,x in grid];input_kind=:analytic)
    first=[(1+.07i)*(1+.03b) for i in 1:10,b in radial.b]
    second=[(1+.11i)*(1+.05b) for i in 1:10,b in radial.b]
    polarized_leptonic_plan(plan,first,second,operator,samples;charges=(1,-1),
        provenance=Dict{String,Any}("fixture_scope"=>"explicitly_synthetic_fixed_rows_no_native_evolution"))
end

function qt_analytic_moments(z,interval,policy,checks)
    rule=MomentumRule(z*interval[1],z*interval[2];order=64)
    records=[]
    for k in 0:4
        reference=setprecision(policy["finite_difference_reference_bits"]) do
            a,b=BigFloat.(interval);Z=BigFloat(z)
            Z^2*(b^(k+2)-a^(k+2))/(k+2)
        end
        actual=sum(rule.weights.*(rule.nodes./z).^k)
        qt_check!(checks,"z=$z qT=$interval absolute polynomial moment $k",actual,reference,
            policy["analytic_integral_max_scaled"],policy)
        push!(records,Dict("kind"=>"polynomial","power"=>k,"actual"=>actual,"reference_decimal"=>string(reference)))
    end
    for (kind,power,rate) in (("Gaussian_UU",0,.7),("Gaussian_Collins",2,.9))
        reference=setprecision(policy["finite_difference_reference_bits"]) do
            a,b=BigFloat.(interval);Z=BigFloat(z);r=BigFloat(rate)
            power==0 ? Z^2*(exp(-r*a^2)-exp(-r*b^2))/(2r) :
                Z^2*((1+r*a^2)*exp(-r*a^2)-(1+r*b^2)*exp(-r*b^2))/(2r^2)
        end
        q=rule.nodes./z
        actual=sum(rule.weights.*q.^power.*exp.(-rate.*q.^2))
        qt_check!(checks,"z=$z qT=$interval absolute $kind integral",actual,reference,
            policy["analytic_integral_max_scaled"],policy)
        push!(records,Dict("kind"=>kind,"power"=>power,"rate"=>rate,"actual"=>actual,"reference_decimal"=>string(reference)))
    end
    records
end

qt_richardson(coarse,fine)=(4 .*fine.-coarse)./3

function qt_derivatives(prediction,reference,theta,direction,policy,fd_policy,checks,label)
    J=ForwardDiff.jacobian(prediction,theta)
    H=reduce(vcat,[permutedims(ForwardDiff.hessian(t->prediction(t)[i],theta)*direction) for i in 1:2])
    refJ=ForwardDiff.jacobian(reference,theta)
    refH=reduce(vcat,[permutedims(ForwardDiff.hessian(t->reference(t)[i],theta)*direction) for i in 1:2])
    qt_check!(checks,label*" all13 independent qT gradient",J,refJ,policy["derivative_max_scaled"],policy)
    qt_check!(checks,label*" all13 independent qT HVP",H,refH,policy["HVP_max_scaled"],policy)
    qt_check!(checks,label*" SIA transversity gradient columns are zero",J[:,1:6],zeros(2,6),0.,policy)
    records=[]
    for step in fd_policy["central_steps"]
        record=setprecision(fd_policy["reference_bits"]) do
            t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
            stencils=[];fd=zeros(BigFloat,2,13)
            for j in 1:13
                delta=zeros(BigFloat,13);delta[j]=h
                plus=reference(t+delta);minus=reference(t-delta)
                fd[:,j]=(plus-minus)/(2h)
                push!(stencils,Dict("parameter"=>j,"plus_decimal"=>string.(plus),"minus_decimal"=>string.(minus)))
            end
            plusJ=ForwardDiff.jacobian(reference,t+h*d)
            minusJ=ForwardDiff.jacobian(reference,t-h*d)
            fdh=(plusJ-minusJ)/(2h)
            Dict("step_binary64"=>step,"step_decimal"=>string(h),"reference_bits"=>precision(h),
                "gradient"=>qt_rows(Float64.(fd)),"HVP"=>qt_rows(Float64.(fdh)),
                "gradient_decimal"=>qt_rows(string.(fd)),"HVP_decimal"=>qt_rows(string.(fdh)),
                "value_stencils"=>stencils,"directional_plus_gradient_decimal"=>qt_rows(string.(plusJ)),
                "directional_minus_gradient_decimal"=>qt_rows(string.(minusJ)))
        end
        fd=reduce(vcat,permutedims.(record["gradient"]));fdh=reduce(vcat,permutedims.(record["HVP"]))
        # Retain raw misses explicitly; only the prospectively frozen fourth-order method is a new gate.
        record["raw_gradient_scaled_error"]=observable_scaled(vec(J),vec(fd),policy)
        record["raw_HVP_scaled_error"]=observable_scaled(vec(H),vec(fdh),policy)
        record["raw_gradient_target_met"]=record["raw_gradient_scaled_error"]<=policy["derivative_max_scaled"]
        record["raw_HVP_target_met"]=record["raw_HVP_scaled_error"]<=policy["HVP_max_scaled"]
        push!(records,record)
    end
    fourth=[]
    for level in 1:2
        record=setprecision(fd_policy["reference_bits"]) do
            matrices(key,i)=parse.(BigFloat,reduce(vcat,permutedims.(records[i][key])))
            gradient=qt_richardson(matrices("gradient_decimal",level),matrices("gradient_decimal",level+1))
            hvp=qt_richardson(matrices("HVP_decimal",level),matrices("HVP_decimal",level+1))
            Dict("level"=>level,"method_order"=>4,"coarse_step"=>records[level]["step_binary64"],
                "fine_step"=>records[level+1]["step_binary64"],"gradient"=>qt_rows(Float64.(gradient)),
                "HVP"=>qt_rows(Float64.(hvp)),"gradient_decimal"=>qt_rows(string.(gradient)),"HVP_decimal"=>qt_rows(string.(hvp)))
        end
        gradient=reduce(vcat,permutedims.(record["gradient"]));hvp=reduce(vcat,permutedims.(record["HVP"]))
        qt_check!(checks,label*" all13 fourth-order FD gradient level=$level",J,gradient,policy["derivative_max_scaled"],policy)
        qt_check!(checks,label*" fourth-order FD HVP level=$level",H,hvp,policy["HVP_max_scaled"],policy)
        push!(fourth,record)
    end
    qt_check!(checks,label*" adjacent fourth-order gradient stability",
        reduce(vcat,permutedims.(fourth[1]["gradient"])),reduce(vcat,permutedims.(fourth[2]["gradient"])),policy["derivative_max_scaled"],policy)
    qt_check!(checks,label*" adjacent fourth-order HVP stability",
        reduce(vcat,permutedims.(fourth[1]["HVP"])),reduce(vcat,permutedims.(fourth[2]["HVP"])),policy["HVP_max_scaled"],policy)
    Dict("gradient"=>qt_rows(J),"HVP"=>qt_rows(H),"reference_gradient"=>qt_rows(refJ),
        "reference_HVP"=>qt_rows(refH),"direction"=>direction,"finite_differences"=>records,"fourth_order_references"=>fourth,
        "acceptance_reference_method"=>"two fourth-order Richardson levels and adjacent-level stability",
        "raw_central_stencils_are_diagnostics"=>true)
end

function validate_sia_qt_windows(producer)
    sources=sia_qt_sources(producer);contract=JSON.parsefile(QT_CONTRACT)
    policy=contract["validation_policy"];settings=contract["synthetic_fixture"]
    fd_contract=JSON.parsefile(QT_FD_CONTRACT);fd_policy=fd_contract["policy"]
    fd_policy["central_steps"][1:2]==policy["finite_difference_steps"] &&
        fd_policy["reference_bits"]==policy["finite_difference_reference_bits"] &&
        all(fd_policy[k]==policy[k] for k in ("relative_floor","derivative_max_scaled","HVP_max_scaled")) ||
        error("FD reference addendum must preserve original steps, precision and acceptance targets")
    theta=Float64.(settings["parameters"]);direction=Float64.(settings["direction"])
    checks=Dict{String,Any}[];cases=[]
    angular=SIAAngularFactor(settings["angular_factor"];source="synthetic acceptance angular ratio")
    mean=SIAMeanOfRatio(settings["angular_factor"];source="synthetic event angular average",
        averaging_measure="explicit synthetic fixed event measure; no experimental admission")
    for z in Float64.(settings["z1_values"]),bounds in settings["qT_windows_GeV"]
        interval=Tuple(Float64.(bounds));context=sia_qt_fixture(z,settings)
        window=SIAQTWindow(interval...;source="KPSY-inspired explicit synthetic qT window",
            lower_limit_source=interval[1]==0 ? "caller-declared synthetic radial origin" : "caller-declared synthetic positive lower bound")
        reference=prepare_generic_reference(context;momentum_interval=interval,momentum_coordinate=:qT,angular=angular.value)
        levels=[];pair=nothing
        for order in policy["momentum_orders"]
            pair=[sia_qt_observable(context,window;ratio=r,angular=angular,order=order) for r in (:UL,:UC)]
            pol=collect(CollinsFit.observable_polarized_integrals(first(pair),KPSYParameters(theta)))
            values=observable_values(pair,theta);operands=vcat(first(pair).denominators,pol,values)
            error=isempty(levels) ? nothing : observable_scaled(operands,last(levels)["absolute_yields_and_ratios"],policy)
            accepted=order>=policy["minimum_momentum_order"] && error!==nothing && error<=policy["momentum_refinement_max_scaled"]
            push!(levels,Dict("order"=>order,"absolute_yields_and_ratios"=>operands,"adjacent_error"=>error,"accepted"=>accepted))
            accepted && break
        end
        label="z=$z qT=$interval"
        qt_assert!(checks,label*" momentum refinement accepted",last(levels)["accepted"],policy)
        obs=first(pair)
        qt_check!(checks,label*" mapped P endpoints",[obs.momentum.lower,obs.momentum.upper],z.*collect(interval),0.,policy)
        qt_check!(checks,label*" absolute projected UU radial weights",obs.uu_weights,reference.uu_weights,policy["adaptive_integral_max_scaled"],policy)
        qt_check!(checks,label*" absolute projected Collins radial weights",obs.collins_weights,reference.collins_weights,policy["adaptive_integral_max_scaled"],policy)
        refden=[reference.denominators[2]+reference.denominators[3],reference.denominators[1]+reference.denominators[4]]
        refpol=reference.polarized_integrals(theta)
        refnum=[refpol[2]+refpol[3],refpol[1]+refpol[4]]
        actualnum=collect(CollinsFit.observable_polarized_integrals(obs,KPSYParameters(theta)))
        qt_check!(checks,label*" absolute U/L unpolarized yields",obs.denominators,refden,policy["adaptive_integral_max_scaled"],policy)
        qt_check!(checks,label*" signed U/L polarized yields",actualnum,refnum,policy["adaptive_integral_max_scaled"],policy)
        prediction(t)=observable_values(pair,t)
        qt_check!(checks,label*" separate-yield UL/UC ratios",prediction(theta),reference.predict(theta),policy["adaptive_integral_max_scaled"],policy)
        alternate=[sia_qt_observable(context,window;ratio=r,angular=mean,order=last(levels)["order"]) for r in (:UL,:UC)]
        qt_check!(checks,label*" equal scalar preserves both typed angular paths",observable_values(alternate,theta),prediction(theta),0.,policy)
        qt_assert!(checks,label*" angular source kinds remain distinct",
            obs.provenance["angular_definition"]=="ratio_of_angular_means" && first(alternate).provenance["angular_definition"]=="mean_of_angular_ratio" &&
            first(alternate).provenance["angular_averaging_measure"]==mean.averaging_measure,policy)
        derivatives=qt_derivatives(prediction,reference.predict,theta,direction,policy,fd_policy,checks,label)
        qt_check!(checks,label*" equal scalar angular paths all13 gradient",ForwardDiff.jacobian(t->observable_values(alternate,t),theta),
            reduce(vcat,permutedims.(derivatives["gradient"])),0.,policy)
        # A common omitted z^2 cancels from ratios, so absolute yield checks are essential.
        wrongden=refden./z^2;wrongnum=refnum./z^2
        wrongratios=angular.value.*[wrongnum[1]/wrongden[1]-wrongnum[2]/wrongden[2],
            wrongnum[1]/wrongden[1]-sum(wrongnum)/sum(wrongden)]
        omission_error=observable_scaled(wrongden,refden,policy)
        qt_assert!(checks,label*" missing-Jacobian absolute-yield counterfactual rejects",omission_error>policy["adaptive_integral_max_scaled"],policy)
        qt_check!(checks,label*" missing common Jacobian can leave ratios unchanged",wrongratios,reference.predict(theta),policy["algebra_max_scaled"],policy)
        push!(cases,Dict("z1"=>z,"z2"=>settings["z2"],"qT_bounds_GeV"=>collect(interval),"levels"=>levels,
            "parameters"=>theta,"provenance"=>deepcopy(obs.provenance),"alternate_angular_provenance"=>deepcopy(first(alternate).provenance),
            "prediction"=>prediction(theta),"reference_prediction"=>reference.predict(theta),
            "derivatives"=>derivatives,"analytic_moments"=>qt_analytic_moments(z,interval,policy,checks),
            "missing_Jacobian_counterfactual"=>Dict("unpolarized_yields"=>wrongden,"polarized_yields"=>wrongnum,
                "ratios"=>wrongratios,"absolute_yield_error"=>omission_error)))
    end
    correct=1/3;wrong=quadgk(P->P*P^2/(1+P^2),0.,1.;rtol=1e-13)[1]/.5
    qt_check!(checks,"analytic average-of-ratios reference",wrong,1-log(2),policy["analytic_integral_max_scaled"],policy)
    qt_assert!(checks,"average-of-pointwise-ratios counterfactual rejects",abs(correct-wrong)>policy["algebra_max_scaled"],policy)
    qt_assert!(checks,"source bytes unchanged throughout complete Julia check",sia_qt_sources(producer)==sources,policy)
    Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
        "input_scope"=>"explicitly_synthetic_fixed_rows_and_uniform_momentum_acceptance",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,"experimental_data_admitted"=>false,
        "historical_fit_ready"=>false,"full_runtime_revalidated"=>false,"native_evolution_accuracy_verified"=>false,
        "contract_id"=>contract["id"],"contract_sha256"=>bytes2hex(sha256(read(QT_CONTRACT))),
        "FD_reference_contract_id"=>fd_contract["id"],"FD_reference_contract_sha256"=>bytes2hex(sha256(read(QT_FD_CONTRACT))),
        "original_raw_reference_report"=>fd_contract["reason"],
        "sources"=>sources,"checks"=>checks,"cases"=>cases,"limitations"=>vcat(contract["limitations"],fd_contract["limitations"]),
        "ratio_counterfactual"=>Dict("ratio_of_integrals"=>correct,"average_of_ratios"=>wrong),"julia_version"=>string(VERSION))
end
