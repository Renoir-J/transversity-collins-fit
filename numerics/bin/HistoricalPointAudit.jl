module HistoricalPointAudit
using CollinsFit, JSON, SHA

const ROOT=normpath(joinpath(@__DIR__,"../.."))
require(value,message)=value || error(message)
digest(path)=bytes2hex(sha256(read(path)))
scaled(a,b,floor)=maximum(abs.(Float64.(a).-Float64.(b))./max.(abs.(Float64.(b)),floor))
packed(rows)=reduce(vcat,[[row["UU"],row["Collins"],row["amplitude"]] for row in rows])

function audit(requestfile,reportfile)
    request=load_historical_point_request(requestfile);r=JSON.parsefile(reportfile)
    contract=JSON.parsefile(joinpath(ROOT,"config/historical-point-contract.json"));policy=contract["numerical_policy"]
    require(r["status"]=="PASS" && r["scope"]==request.data["scope"] && r["theory_scope"]=="reconstructed_theory" &&
        r["observation_scope"]=="historical_benchmark_candidate","unaccepted diagnostic scope/status")
    for key in ("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted",
                "fit_performed","likelihood_constructed","G2PERF_pass","interval_coverage_verified")
        require(get(r,key,nothing)===false,"diagnostic cannot approve a fit/production gate")
    end
    require(r["request_sha256"]==request.sha256 && r["contract"]==contract &&
        r["contract_sha256"]==digest(joinpath(ROOT,"config/historical-point-contract.json")) &&
        r["source_fingerprint"]==request.data["source_fingerprint"] &&
        r["request_points"]==request.data["points"] && r["parameters"]==request.data["parameters"] &&
        r["model_convention"]==request.data["model_convention"] &&
        r["source_preparation"]==request.data["source_preparation"] &&
        r["source_preparation_independently_replayed_in_Julia"]===false &&
        r["input_kind"]=="source_preparation_attested_rows_and_pinned_parton_inputs" &&
        r["limitations"]==contract["limitations"],"changed source/request/model identity")
    bindings=r["sources"]
    for relative in ("numerics/bin/predict_historical_points.jl","numerics/bin/HistoricalPointAudit.jl",
                     "numerics/src/HistoricalPoints.jl","numerics/bin/ObservableReference.jl",
                     "config/historical-point-contract.json","src/collinsflow/historical_points.py")
        require(haskey(bindings,relative),"missing essential source binding")
    end
    for (relative,sha) in bindings
        path=normpath(joinpath(ROOT,relative))
        require(!isabspath(relative) && !startswith(relpath(path,ROOT),"..") && isfile(path) && digest(path)==sha,"stale diagnostic source")
    end
    for (path,sha) in r["dependency_sha256"]
        require(isfile(path) && digest(path)==sha,"changed diagnostic dependency")
    end
    # Hash the current native/source/input files again without rerunning evolution.
    dependencies=r["dependency_sha256"]
    require(length(dependencies)==3 && get(dependencies,abspath(requestfile),nothing)==request.sha256,
        "incomplete request/build dependency bindings")
    nativefiles=[path for path in keys(dependencies) if path!=abspath(requestfile) &&
        get(JSON.parsefile(path),"ABI_version",nothing)==1]
    inputfiles=[path for path in keys(dependencies) if path!=abspath(requestfile) &&
        haskey(JSON.parsefile(path),"LHAPDF_version")]
    require(length(nativefiles)==length(inputfiles)==1,"ambiguous build dependency identity")
    native=CollinsFit.verified_hoppet_build(only(nativefiles));inputs=CollinsFit.verify_parton_build(only(inputfiles))
    require(native==r["native_build"] && inputs==r["input_build"],"changed recorded build metadata")
    require(any(basename(path)=="hoppet_bridge.f90" && sha==digest(joinpath(ROOT,"fortran/hoppet_bridge.f90"))
        for (path,sha) in native["source_hashes"]),"stale native bridge")
    require(any(basename(path)=="lhapdf.cpp" && sha==digest(joinpath(ROOT,"bridges/lhapdf.cpp"))
        for (path,sha) in inputs["source_hashes"]),"stale parton bridge")
    points=request.data["points"];ids=[p["id"] for p in points];floor=policy["relative_floor"]
    expected=Dict{String,Tuple{Float64,Float64,Bool}}()
    add(name,error,tolerance;attestation=false)=(expected[name]=(Float64(error),Float64(tolerance),attestation))
    boolean(name,value;attestation=false)=add(name,value ? 0. : 1.,0.;attestation=attestation)
    ys=unique(Float64[p["kinematics"]["y"] for p in points]);pairs=[(y,sign) for y in ys for sign in (-1,1)]
    projections=r["source_projection"];require(length(projections)==length(pairs),"incomplete source projection grid")
    for (i,(row,(y,sign))) in enumerate(zip(projections,pairs))
        require(row["scope"]=="explicitly_synthetic_source_definition_fixture" && row["y"]==y &&
            row["R"]==sign/8 && row["spin"]==.5 && row["Sivers"]==1/16,"changed synthetic source-definition fixture")
        closed=setprecision(policy["source_projection_reference_bits"]) do
            yy=BigFloat(y);factor=2(1-yy)/(1+(1-yy)^2);R=BigFloat(sign)/8
            (expected=Float64(factor*R),positive=BigFloat(1)/2*(abs(factor*R)+BigFloat(1)/16)<1,
             counterfactuals=Dict("wrong_phase"=>Float64(-factor*R),"omitted_DNN"=>Float64(R),
                 "repeated_DNN"=>Float64(factor^2*R),"factor2"=>Float64(2factor*R)))
        end
        require(row["reference"]==closed.expected,"changed analytic Fourier coefficient")
        implementation=2(1-y)/(1+(1-y)^2)*(sign/8)
        require(row["implementation_factor_times_R"]==implementation,"changed source-factor operand")
        require(Set(keys(row["counterfactuals"]))==Set(keys(closed.counterfactuals)),"incomplete source counterfactuals")
        for (name,value) in closed.counterfactuals
            require(scaled([row["counterfactuals"][name]],[value],floor)<=policy["observable_algebra_max_scaled"],"changed counterfactual operand")
        end
        distinct=all(abs(x-row["projected"])>floor for x in values(row["counterfactuals"]))
        require(row["positive_counts_proved_by_bound"]===closed.positive && row["all_counterfactuals_distinct"]===distinct,"changed count/phase control")
        error=scaled([implementation],[row["projected"]],floor)
        require(row["error"]==error,"changed projection comparison error")
        add("source Eq1 independent two-angle projection $i",scaled([row["projected"]],[row["reference"]],floor),policy["source_projection_max_scaled"])
        add("source convention application $i",error,policy["source_projection_max_scaled"])
        boolean("source fixture positive/distinct $i",closed.positive && distinct)
    end
    function rows_valid(rows)
        require([row["id"] for row in rows]==ids,"changed numerical point order")
        for row in rows
            require(all(v->v isa Real && !(v isa Bool) && isfinite(v),[row["UU"],row["Collins"],row["amplitude"]]) && row["UU"]>0,"invalid signed numerical operand")
        end
    end
    levels=r["refinements"]
    require(policy["minimum_refinement_level"]<=length(levels)<=length(policy["refinements"]),"invalid refinement budget")
    for (i,level) in enumerate(levels)
        require(level["level"]==i && level["setting"]==policy["refinements"][i] && level["bound_inverse_GeV"]==24,"changed frozen refinement setting")
        rows_valid(level["points"])
        error=i==1 ? nothing : scaled(packed(level["points"]),packed(levels[i-1]["points"]),floor)
        require(level["previous_level_max_scaled"]==error,"changed adjacent refinement error")
        target=error!==nothing && error<=policy["physical_value_max_scaled"]
        accepted=i>=policy["minimum_refinement_level"] && target
        require(level["target_met"]===target && level["accepted_refinement"]===accepted,"changed refinement acceptance bit")
        require(!accepted || i==length(levels),"continued past accepted frozen refinement")
    end
    boolean("all eight point joint refinements accepted",last(levels)["accepted_refinement"])
    tails=r["tails"]
    require([row["bound_inverse_GeV"] for row in tails]==policy["radial_bounds_inverse_GeV"],"changed radial bounds")
    for row in tails;rows_valid(row["points"]);end
    require(tails[2]["points"]==last(levels)["points"],"changed retained bound24 values")
    add("all eight point radial extension24to32",scaled(packed(tails[end]["points"]),packed(tails[end-1]["points"]),floor),policy["physical_value_max_scaled"])
    references=r["independent_references"];reuse=r["charge_row_reuse_producer_attestations"]
    require(length(references)==length(reuse)==length(points),"incomplete reference/reuse ledger")
    for (i,(point,final,reference,reused)) in enumerate(zip(points,r["points"],references,reuse))
        for key in keys(point);require(final[key]==point[key],"changed measured/source field");end
        numerical=tails[end]["points"][i]
        for key in keys(numerical);require(final[key]==numerical[key],"changed final prediction field");end
        require(reference["id"]==reused["id"]==point["id"],"changed reference identity")
        require(reference["scope"]=="independent_fixed_momentum_scalar_assembly" &&
            reference["momentum_interval"]===nothing && reference["native_rows_shared"]===true &&
            reference["initial_functions_shared"]===true && reference["independent_evolution"]===false,
            "changed fixed-row reference scope")
        actual=[final["UU"],final["Collins"],final["amplitude"]]
        require(reference["actual"]==actual && length(reference["reference"])==3,"changed independent comparison operands")
        y=point["kinematics"]["y"];factor=2(1-y)/(1+(1-y)^2)
        reference_values=reference["reference"]
        require(all(value->value isa Real && !(value isa Bool) && isfinite(value),reference_values) &&
            reference_values[1]>0 &&
            scaled([reference["reference_structure_ratio"]],[reference_values[2]/reference_values[1]],floor)<=policy["observable_algebra_max_scaled"] &&
            scaled([reference_values[3]],[factor*reference_values[2]/reference_values[1]],floor)<=policy["observable_algebra_max_scaled"],
            "changed independent ratio/normalization operands")
        algebra=factor*final["Collins"]/final["UU"]
        require(final["DNN_factor"]==factor && final["algebra_reference"]==algebra,"changed observable algebra operands")
        require(final["unweighted_observed_minus_predicted"]==point["measured"]["central_value"]-final["amplitude"],"changed unweighted difference")
        add("independent fixed-row scalar composition "*point["id"],scaled(actual,reference["reference"],floor),policy["physical_value_max_scaled"])
        add("source-model observable algebra "*point["id"],scaled([final["amplitude"]],[algebra],floor),policy["observable_algebra_max_scaled"])
        require(reused["same_native_row_objects"] isa Bool,"missing object-reuse attestation")
        boolean("exact native charge-row reuse "*point["id"],reused["same_native_row_objects"];attestation=true)
    end
    sign=r["global_sign_control"];require(sign["original"]==[p["amplitude"] for p in r["points"]],"changed global-sign base")
    add("global transversity Collins sign invariance",scaled(sign["flipped"],sign["original"],floor),policy["observable_algebra_max_scaled"])
    boolean("source builds request stable during complete job",true;attestation=true)
    checks=r["checks"]
    require(length(checks)==length(expected) && Set(c["name"] for c in checks)==Set(keys(expected)),"incomplete diagnostic check ledger")
    for check in checks
        error,tolerance,_=expected[check["name"]]
        require(check["error"]==error && check["tolerance"]==tolerance &&
            check["passed"]===(isfinite(error) && error<=tolerance) && check["passed"]===true,"changed or failed check bit/target")
    end
    require(r["final_bound_inverse_GeV"]==32 && r["unique_prepared_geometries"]==length(unique(historical_geometry_key.(points))),"changed final geometry/bound count")
    historical_request_valid(request)
    Dict("status"=>r["status"],"review_status"=>"STORED_DIAGNOSTIC_ARITHMETIC_VERIFIED",
        "scope"=>request.data["scope"],"accepted_prediction_diagnostic"=>true,"accepted_fitting_gate"=>false,
        "production_ready"=>false,"physical_extraction"=>false,"historical_fit_ready"=>false,"experimental_data_admitted"=>false,
        "fit_performed"=>false,"G2PERF_pass"=>false,"interval_coverage_verified"=>false,
        "report_sha256"=>digest(reportfile),"request_sha256"=>request.sha256,
        "source_fingerprint"=>request.data["source_fingerprint"],"source_bindings_verified"=>length(bindings),
        "stored_comparison_errors_recomputed"=>count(!last(value) for value in values(expected)),
        "producer_object_reuse_and_temporal_attestations"=>count(last(value) for value in values(expected)),
        "native_evolution_or_prediction_calls_replayed"=>false,"two_angle_quadrature_reexecuted"=>false,
        "current_build_source_and_binary_hashes_reverified"=>true,
        "raw_primary_tables_independently_reparsed_in_Julia"=>false,
        "source_preparation_kind"=>"Python source-replay attestation with current input hashes",
        "auditor_sha256"=>digest(@__FILE__))
end
end
