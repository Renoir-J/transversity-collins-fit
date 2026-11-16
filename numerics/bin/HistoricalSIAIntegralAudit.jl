module HistoricalSIAIntegralAudit
using CollinsFit, JSON, SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const PAIRS=("pp","pm","mp","mm")
require(value,message)=value || error(message)
digest(p)=bytes2hex(sha256(read(p)))
finite(v)=v isa Real && !(v isa Bool) && isfinite(v)
scaled(a,b,floor)=maximum(abs.(Float64.(a).-Float64.(b))./max.(abs.(Float64.(b)),floor))
packed_point(p)=vcat([p["pairs"][pair][k] for pair in PAIRS for k in ("UU","Collins")],
    [p["categories"][c][k] for c in ("U","L","C") for k in ("UU","Collins")],[p["amplitude"]])
packed(rows)=reduce(vcat,packed_point.(rows))

function yields_valid(row)
    require(Set(keys(row["pairs"]))==Set(PAIRS) && Set(keys(row["categories"]))==Set(("U","L","C")),"incomplete charge/category operands")
    for item in vcat(collect(values(row["pairs"])),collect(values(row["categories"])))
        require(Set(keys(item))==Set(("UU","Collins")) && all(finite,values(item)) && item["UU"]>0,"invalid signed yields or denominator")
    end
    for k in ("UU","Collins")
        u=row["pairs"]["pm"][k]+row["pairs"]["mp"][k]
        l=row["pairs"]["pp"][k]+row["pairs"]["mm"][k]
        require(row["categories"]["U"][k]==u && row["categories"]["L"][k]==l && row["categories"]["C"][k]==u+l,
            "changed independent category sum")
    end
    require(finite(row["amplitude"]),"nonfinite SIA amplitude")
end

function closed_amplitude(row,point)
    u=row["categories"]["U"];right=row["categories"][point["ratio"]=="UL" ? "L" : "C"]
    point["source_angular"]["value"]*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"])
end

function inward_checks(derivative,dp)
    checks=Tuple{String,Float64,Float64}[];id=""
    add(name,error,tolerance)=push!(checks,(name,Float64(error),Float64(tolerance)))
    boolean(name,value)=add(name,value ? 0. : 1.,0.)
    numeric_matrix(data)=begin
        require(length(data)==2 && all(length(row)==13 && all(finite,row) for row in data),"invalid all13 derivative matrix")
        Float64.(reduce(vcat,permutedims.(data)))
    end
    decimal_matrix(data)=begin
        require(length(data)==2 && all(length(row)==13 && all(v->v isa AbstractString,row) for row in data),"invalid decimal derivative matrix")
        matrix=parse.(BigFloat,reduce(vcat,permutedims.(data)))
        require(all(isfinite,matrix),"nonfinite decimal derivative matrix");matrix
    end
    J=numeric_matrix(derivative["gradient"]);H=numeric_matrix(derivative["HVP"])
    RJ=numeric_matrix(derivative["reference_gradient"]);RH=numeric_matrix(derivative["reference_HVP"])
    gerror=scaled(vec(J),vec(RJ),dp["relative_floor"]);herror=scaled(vec(H),vec(RH),dp["relative_floor"])
    require(derivative["independent_gradient_error"]==gerror && derivative["independent_HVP_error"]==herror,"changed derivative comparison arithmetic")
    add(id*" independent all13 gradient",gerror,dp["derivative_max_scaled"])
    add(id*" independent all13 HVP",herror,dp["HVP_max_scaled"])
    gzero=all(iszero,J[:,1:6]) && all(iszero,RJ[:,1:6]);hzero=all(iszero,H[:,1:6]) && all(iszero,RH[:,1:6])
    require(derivative["first_six_gradient_columns_zero"]===gzero && derivative["first_six_HVP_columns_zero"]===hzero,"changed zero-column decision")
    boolean(id*" first six derivative columns zero",gzero && hzero)
    fd=derivative["finite_differences"];require(length(fd)==length(dp["steps"])==2,"incomplete inward stencil levels")
    for (index,record) in enumerate(fd)
        FJ,FH=setprecision(dp["reference_bits"]) do
            h=parse(BigFloat,string(dp["steps"][index]))
            require(record["step"]==dp["steps"][index] && record["step_decimal"]==string(h),"changed inward step")
            columns=record["gradient_samples_decimal"]
            require(length(columns)==13 && all(length(c)==5 && all(length(v)==2 && all(x->x isa AbstractString,v) for v in c) for c in columns),"incomplete gradient samples")
            gradients=zeros(BigFloat,2,13)
            for j in 1:13
                c=[parse.(BigFloat,v) for v in columns[j]];require(all(all(isfinite,v) for v in c),"nonfinite gradient sample")
                for i in 1:2
                    gradients[i,j]=(-25c[1][i]+48c[2][i]-36c[3][i]+16c[4][i]-3c[5][i])/(12h)
                end
            end
            samples=record["HVP_jacobian_samples_decimal"];require(length(samples)==5,"incomplete HVP stencil")
            m=decimal_matrix.(samples)
            hvps=[(-25m[1][i,j]+48m[2][i,j]-36m[3][i,j]+16m[4][i,j]-3m[5][i,j])/(12h) for i in 1:2,j in 1:13]
            require(gradients==decimal_matrix(record["gradient_decimal"]) && hvps==decimal_matrix(record["HVP_decimal"]),"changed inward derivative arithmetic")
            Float64.(gradients),Float64.(hvps)
        end
        require(FJ==numeric_matrix(record["gradient"]) && FH==numeric_matrix(record["HVP"]),"changed derivative float projection")
        ge=scaled(vec(J),vec(FJ),dp["relative_floor"]);he=scaled(vec(H),vec(FH),dp["relative_floor"])
        require(record["gradient_error"]==ge && record["HVP_error"]==he,"changed inward comparison error")
        add(id*" inward gradient level=$index",ge,dp["derivative_max_scaled"])
        add(id*" inward HVP level=$index",he,dp["HVP_max_scaled"])
    end
    ge=scaled(vec(numeric_matrix(fd[2]["gradient"])),vec(numeric_matrix(fd[1]["gradient"])),dp["relative_floor"])
    he=scaled(vec(numeric_matrix(fd[2]["HVP"])),vec(numeric_matrix(fd[1]["HVP"])),dp["relative_floor"])
    require(derivative["adjacent_gradient_error"]==ge && derivative["adjacent_HVP_error"]==he,"changed inward adjacent errors")
    add(id*" adjacent inward gradient",ge,dp["derivative_max_scaled"])
    add(id*" adjacent inward HVP",he,dp["HVP_max_scaled"])
    checks
end

function audit(requestfile,reportfile)
    request=load_historical_sia_integral_request(requestfile);r=JSON.parsefile(reportfile)
    contract=JSON.parsefile(joinpath(ROOT,"config/historical-sia-integral-contract.json"));policy=contract["numerical_policy"]
    require(r["status"]=="PASS" && r["scope"]==request.data["scope"] && r["theory_scope"]=="reconstructed_theory" &&
        r["observation_scope"]=="historical_benchmark_candidate","unaccepted SIA diagnostic scope/status")
    for key in (CollinsFit.HISTORICAL_SIA_INTEGRAL_FALSE_FLAGS...,"likelihood_constructed")
        require(get(r,key,nothing)===false,"SIA diagnostic cannot approve fitting or production")
    end
    require(r["request_sha256"]==request.sha256 && r["contract"]==contract &&
        r["contract_sha256"]==digest(joinpath(ROOT,"config/historical-sia-integral-contract.json")) &&
        r["source_fingerprint"]==request.data["source_fingerprint"] && r["request_points"]==request.data["points"] &&
        all(r[k]==request.data[k] for k in ("parameters","model_convention","momentum_window","source_preparation","limitations")) &&
        r["source_preparation_independently_replayed_in_Julia"]===false &&
        r["raw_primary_upper_cut_reparsed_in_Julia"]===false &&
        r["input_kind"]=="source_preparation_attested_rows_and_pinned_parton_inputs","changed SIA source/request identity")
    bindings=r["sources"]
    require(bindings==CollinsFit.historical_sia_integral_source_bindings(),"incomplete or stale exact SIA source bindings")
    dependencies=r["dependency_sha256"]
    require(length(dependencies)==3 && get(dependencies,abspath(requestfile),nothing)==request.sha256,"incomplete request/build dependencies")
    for (path,sha) in dependencies
        require(isfile(path) && digest(path)==sha,"changed SIA diagnostic dependency")
    end
    nativefiles=[p for p in keys(dependencies) if p!=abspath(requestfile) && get(JSON.parsefile(p),"ABI_version",nothing)==1]
    inputfiles=[p for p in keys(dependencies) if p!=abspath(requestfile) && haskey(JSON.parsefile(p),"LHAPDF_version")]
    require(length(nativefiles)==length(inputfiles)==1,"ambiguous SIA build identities")
    native=CollinsFit.verified_hoppet_build(only(nativefiles));inputs=CollinsFit.verify_parton_build(only(inputfiles))
    require(native==r["native_build"] && inputs==r["input_build"],"changed build metadata")
    require(any(basename(p)=="hoppet_bridge.f90" && h==digest(joinpath(ROOT,"fortran/hoppet_bridge.f90")) for (p,h) in native["source_hashes"]),"stale native bridge")
    require(any(basename(p)=="lhapdf.cpp" && h==digest(joinpath(ROOT,"bridges/lhapdf.cpp")) for (p,h) in inputs["source_hashes"]),"stale parton bridge")
    require(r["active_parton_input"]==Dict("provider"=>"DSS14","CT10_and_DSSV"=>"build verified; not used by SIA predictions"),"changed active provider declaration")
    require(r["preparation_reuse"]==Dict("kind"=>"unchanged_per_geometry_preparation","new_query_cache_used"=>false,
        "parameter_responses_retained_across_calls"=>false),"changed preparation method")
    points=request.data["points"];ids=[p["id"] for p in points];floor=policy["relative_floor"]
    expected=Dict{String,Tuple{Float64,Float64,Bool}}()
    add(name,error,tolerance;attestation=false)=(expected[name]=(Float64(error),Float64(tolerance),attestation))
    boolean(name,value;attestation=false)=add(name,value ? 0. : 1.,0.;attestation=attestation)
    projections=r["source_projection"];grid=[(point,sign) for point in points for sign in (-1,1)]
    require(length(projections)==length(grid),"incomplete source angular grid")
    for (i,(row,(point,sign))) in enumerate(zip(projections,grid))
        f=point["source_angular"]["value"]
        require(row["scope"]=="explicitly_synthetic_source_definition_fixture" && row["point_id"]==point["id"] &&
            row["ratio"]==point["ratio"] && row["angular_factor"]==f && row["sign"]==sign &&
            row["D_U"]==2. && row["D_L"]==1. && row["C_U"]==sign/8 && row["C_L"]==-sign/16 &&
            row["epsilon_interval"]==[-1.,1.] && row["linearization_at_epsilon"]==0.,"changed linearized source fixture")
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
        require(row["reference"]==closed.reference && finite(row["projected"]),"changed analytic angular coefficient")
        du,dl,cu,cl=2.,1.,sign/8,-sign/16
        dr,cr=point["ratio"]=="UL" ? (dl,cl) : (du+dl,cu+cl)
        implementation=f*(cu/du-cr/dr)
        require(row["implementation"]==implementation && Set(keys(row["counterfactuals"]))==Set(keys(closed.counterfactuals)),"changed angular implementation/counterfactuals")
        for (name,value) in closed.counterfactuals
            require(finite(row["counterfactuals"][name]) && scaled([row["counterfactuals"][name]],[value],floor)<=policy["observable_algebra_max_scaled"],"changed angular counterfactual")
        end
        distinct=all(abs(x-row["projected"])>floor for x in values(row["counterfactuals"]))
        require(row["positive_counts_proved_by_bound"]===closed.positive && row["all_counterfactuals_differ_from_correct"]===distinct,"changed positivity/distinctness bit")
        error=scaled([implementation],[row["projected"]],floor);require(row["error"]==error,"changed source application error")
        add("source RF0 independent linearized Fourier projection $i",scaled([row["projected"]],[row["reference"]],floor),policy["source_projection_max_scaled"])
        add("source angular convention application $i",error,policy["source_projection_max_scaled"])
        boolean("source fixture positive/distinct $i",closed.positive && distinct)
    end
    dp=contract["representative_derivative_policy"]
    math=r["stencil_formula_control"]
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
                require(record["function"]=="monomial" && record["degree"]==degree && record["known_stencil_bias"]==string(bias),"changed polynomial control")
            end
            require(record["analytic_derivative"]==string(exact) && record["passed"]===passed,"changed closed-function decision")
            push!(passes,passed)
        end
        halving=abs(errors[1]/errors[2]-16)<big"0.16"
        require(math["exponential_halving_ratio"]==string(errors[1]/errors[2]) && math["fourth_order_halving_verified"]===halving &&
            math["passed"]===(all(passes) && halving),"changed mathematical-control aggregate")
        boolean("inward stencil closed-function controls",all(passes) && halving)
    end
    function rows_valid(rows,order)
        require([p["id"] for p in rows]==ids,"changed numerical point order")
        for (row,point) in zip(rows,points)
            yields_valid(row)
            require(row["angular_factor"]==point["source_angular"]["value"] &&
                row["algebra_reference"]==closed_amplitude(row,point) && row["momentum_order"]==order &&
                row["mapped_P_bounds_GeV"]==point["mapped_P_bounds_GeV"] && all(finite,row["mapped_P_bounds_GeV"]) &&
                row["qT_bounds_GeV"]==request.data["momentum_window"]["qT_bounds_GeV"] && all(finite,row["qT_bounds_GeV"]),
                "changed integrated kinematics or algebra")
        end
    end
    function resources_valid(entry,setting,B)
        resources=entry["resources"];geometries=Set(historical_sia_integral_geometry_key.(points))
        require(length(resources)==length(geometries)==36,"incomplete resource geometry ledger")
        require(all(all(finite,a["geometry"]) for a in resources) && Set(Tuple(a["geometry"]) for a in resources)==geometries,"changed resource geometries")
        for item in resources
            require(finite(item["cold_seconds_informational"]) && item["cold_seconds_informational"]>=0 &&
                item["radial_nodes"] isa Integer && !(item["radial_nodes"] isa Bool) && item["radial_nodes"]>0 &&
                item["row_coefficient_payload_bytes"] isa Integer && !(item["row_coefficient_payload_bytes"] isa Bool) && item["row_coefficient_payload_bytes"]>0,
                "invalid informational resource metadata")
            p=item["context_provenance"]
            require(p["scope"]=="reconstructed_theory" && p["process"]=="SIA" && p["initial_source_kind"]=="pinned_inputs" &&
                p["range_policy"]=="reference" && p["production_ready"]===false && p["physical_extraction"]===false &&
                p["radial_bound_inverse_GeV"]==B && p["radial_order_per_panel"]==setting["radial_order"] &&
                p["matching_order_per_input_cell"]==setting["matching_order"],"changed context provenance")
        end
    end
    function momentum_valid(entry,label)
        levels=entry["momentum_refinements"]
        require(1<=length(levels)<=length(policy["momentum_orders"]),"incomplete momentum ledger")
        for (i,level) in enumerate(levels)
            order=level["order"]
            require(order isa Integer && !(order isa Bool) && order==policy["momentum_orders"][i],"changed momentum refinement order")
            rows_valid(level["points"],order)
            error=i==1 ? nothing : scaled(packed(level["points"]),packed(levels[i-1]["points"]),floor)
            target=error!==nothing && error<=policy["momentum_max_scaled"]
            accepted=order>=policy["minimum_momentum_order"] && target
            require(level["previous_order_max_scaled"]==error && level["target_met"]===target && level["accepted"]===accepted,
                "changed momentum arithmetic or acceptance")
            require(!accepted || i==length(levels),"continued after accepted momentum level")
        end
        require(entry["points"]==last(levels)["points"] && entry["momentum_order"]==last(levels)["order"],"changed retained momentum result")
        boolean("momentum accepted "*label,last(levels)["accepted"])
    end
    levels=r["refinements"]
    require(policy["minimum_refinement_level"]<=length(levels)<=length(policy["refinements"]),"incomplete native refinement ledger")
    for (i,level) in enumerate(levels)
        setting=policy["refinements"][i]
        require(level["level"]==i && level["setting"]==setting && level["bound_inverse_GeV"]==24,"changed native refinement settings")
        resources_valid(level,setting,24);momentum_valid(level,"native=$i")
        error=i==1 ? nothing : scaled(packed(level["points"]),packed(levels[i-1]["points"]),floor)
        target=error!==nothing && error<=policy["physical_value_max_scaled"]
        accepted=i>=policy["minimum_refinement_level"] && target
        require(level["previous_level_max_scaled"]==error && level["target_met"]===target && level["accepted_refinement"]===accepted,
            "changed native refinement arithmetic/decision")
        require(!accepted || i==length(levels),"continued past accepted native level")
    end
    boolean("all seventy-two point joint refinements accepted",last(levels)["accepted_refinement"])
    tails=r["tails"]
    require([t["bound_inverse_GeV"] for t in tails]==policy["radial_bounds_inverse_GeV"],"changed finite radial bounds")
    for entry in tails
        B=entry["bound_inverse_GeV"];resources_valid(entry,last(levels)["setting"],B)
        if B==24
            require(all(entry[k]==last(levels)[k] for k in ("resources","momentum_refinements","points","momentum_order")),"changed reused bound24 record")
        else
            momentum_valid(entry,"tail B=$B")
        end
    end
    add("all seventy-two point radial extension24to32",scaled(packed(tails[3]["points"]),packed(tails[2]["points"]),floor),policy["physical_value_max_scaled"])
    references=r["independent_references"];reuse=r["category_row_reuse_producer_attestations"]
    require(length(references)==length(reuse)==length(r["points"])==length(points)==72,"incomplete final integral/reference ledger")
    for (i,(point,final,ref,reused)) in enumerate(zip(points,r["points"],references,reuse))
        require(all(final[k]==v for (k,v) in point) && all(final[k]==v for (k,v) in tails[3]["points"][i]),"changed final source or integral values")
        require(final["unweighted_observed_minus_predicted"]==point["measured"]["central_value"]-final["amplitude"],"changed unweighted difference")
        require(ref["id"]==reused["id"]==point["id"] && ref["scope"]=="generic_fixed_row_scalar_assembly_and_adaptive_qT_reference" &&
            ref["qT_bounds_GeV"]==[0.,3.5] && ref["angular_value"]==point["source_angular"]["value"] &&
            ref["native_rows_shared"]===true && ref["initial_functions_shared"]===true && ref["independent_evolution"]===false,
            "changed independent reference scope")
        actual=Dict(k=>final[k] for k in ("pairs","categories","amplitude"));reference=ref["reference"]
        require(ref["actual"]==actual,"changed independent comparison operands");yields_valid(reference)
        require(scaled([reference["amplitude"]],[closed_amplitude(reference,point)],floor)<=policy["observable_algebra_max_scaled"],"changed independent integrated ratio")
        add("independent adaptive qT composition "*point["id"],scaled(packed_point(actual),packed_point(reference),floor),policy["physical_value_max_scaled"])
        add("source-model observable algebra "*point["id"],scaled([final["amplitude"]],[closed_amplitude(final,point)],floor),policy["observable_algebra_max_scaled"])
        boolean("positive integrated denominators "*point["id"],true)
        require(reused["same_native_row_objects"] isa Bool,"missing row-reuse attestation")
        boolean("exact native category-row reuse "*point["id"],reused["same_native_row_objects"];attestation=true)
    end
    amplitudes=[p["amplitude"] for p in r["points"]];theta=Float64.(request.data["parameters"]["values"])
    flipped=copy(theta);flipped[[7,8]].*=-1;changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
    sign=r["collins_sign_control"];probe=r["transversity_value_control"]
    require(sign["theta"]==probe["theta"]==theta && sign["original"]==probe["original"]==amplitudes &&
        sign["changed_theta"]==flipped && probe["changed_theta"]==changed &&
        length(sign["flipped"])==length(probe["changed"])==72 && all(finite,sign["flipped"]) && all(finite,probe["changed"]) &&
        probe["scope"]=="fixed_coordinate_value_probe_not_global_derivative_certificate","changed fixed-value controls")
    add("global Collins sign invariance",scaled(sign["flipped"],amplitudes,floor),policy["observable_algebra_max_scaled"])
    add("first six coordinate value independence",scaled(probe["changed"],amplitudes,floor),policy["observable_algebra_max_scaled"])
    derivatives=r["representative_derivatives"]
    require([d["geometry_source_id"] for d in derivatives]==dp["source_geometry_ids"],"changed representative derivative selection")
    for derivative in derivatives
        id=derivative["geometry_source_id"];point=only(filter(p->p["id"]==id,points));key=historical_sia_integral_geometry_key(point)
        selected=[only(filter(p->historical_sia_integral_geometry_key(p)==key && p["ratio"]==ratio,points)) for ratio in ("UL","UC")]
        require(derivative["point_ids"]==[p["id"] for p in selected] && derivative["angular_values"]==[p["source_angular"]["value"] for p in selected] &&
            derivative["scope"]=="three_source_selected_geometries_inward_all13_reference_not_production_certificate" &&
            derivative["theta"]==theta && derivative["direction"]==dp["direction"] && derivative["native_rows_shared"]===true &&
            derivative["independent_evolution"]===false,"changed derivative source/scope")
        for (suffix,error,tolerance) in inward_checks(derivative,dp)
            add(id*suffix,error,tolerance)
        end
    end
    boolean("source builds request stable during complete job",true;attestation=true)
    checks=r["checks"]
    require(length(checks)==length(expected) && Set(c["name"] for c in checks)==Set(keys(expected)),"incomplete or duplicate integral check ledger")
    for check in checks
        error,tolerance,_=expected[check["name"]]
        require(check["error"]==error && check["tolerance"]==tolerance && check["passed"]===(isfinite(error) && error<=tolerance) && check["passed"]===true,
            "changed, failed or missing numerical check")
    end
    require(r["final_bound_inverse_GeV"]==32 && r["unique_prepared_geometries"]==36,"changed final geometry/bound count")
    historical_sia_integral_request_valid(request)
    Dict("status"=>"PASS","review_status"=>"STORED_INTEGRAL_DIAGNOSTIC_ARITHMETIC_VERIFIED","scope"=>request.data["scope"],
        "accepted_prediction_diagnostic"=>true,"accepted_fitting_gate"=>false,
        "production_ready"=>false,"physical_extraction"=>false,"historical_fit_ready"=>false,"experimental_data_admitted"=>false,
        "fit_performed"=>false,"G2PERF_pass"=>false,"interval_coverage_verified"=>false,
        "request_sha256"=>request.sha256,"report_sha256"=>digest(reportfile),"source_fingerprint"=>request.data["source_fingerprint"],
        "source_bindings_verified"=>length(bindings),"stored_comparison_errors_recomputed"=>count(!last(v) for v in values(expected)),
        "producer_object_reuse_and_temporal_attestations"=>count(last(v) for v in values(expected)),
        "native_evolution_or_prediction_calls_replayed"=>false,"Fourier_quadrature_reexecuted"=>false,
        "AD_or_physical_reference_integrals_reexecuted"=>false,"closed_mathematical_functions_reevaluated"=>true,
        "producer_resource_measurements_recomputed"=>false,"raw_primary_tables_independently_reparsed_in_Julia"=>false,
        "current_build_source_and_binary_hashes_reverified"=>true,
        "source_preparation_kind"=>"Python source-replay attestation with current input hashes","auditor_sha256"=>digest(@__FILE__))
end
end
