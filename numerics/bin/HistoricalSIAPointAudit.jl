module HistoricalSIAPointAudit
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

function audit(requestfile,reportfile)
    request=load_historical_sia_point_request(requestfile);r=JSON.parsefile(reportfile)
    contract=JSON.parsefile(joinpath(ROOT,"config/historical-sia-point-contract.json"));policy=contract["numerical_policy"]
    require(r["status"]=="PASS" && r["scope"]==request.data["scope"] && r["theory_scope"]=="reconstructed_theory" &&
        r["observation_scope"]=="historical_benchmark_candidate","unaccepted SIA diagnostic scope/status")
    for key in (CollinsFit.HISTORICAL_SIA_FALSE_FLAGS...,"likelihood_constructed")
        require(get(r,key,nothing)===false,"SIA diagnostic cannot approve fitting or production")
    end
    require(r["request_sha256"]==request.sha256 && r["contract"]==contract &&
        r["contract_sha256"]==digest(joinpath(ROOT,"config/historical-sia-point-contract.json")) &&
        r["source_fingerprint"]==request.data["source_fingerprint"] && r["request_points"]==request.data["points"] &&
        all(r[k]==request.data[k] for k in ("parameters","model_convention","source_preparation","limitations")) &&
        r["source_preparation_independently_replayed_in_Julia"]===false &&
        r["input_kind"]=="source_preparation_attested_rows_and_pinned_parton_inputs","changed SIA source/request identity")
    bindings=r["sources"]
    require(bindings==CollinsFit.historical_sia_source_bindings(),"incomplete or stale exact SIA source bindings")
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
    function rows_valid(rows)
        require([p["id"] for p in rows]==ids,"changed numerical point order")
        for (row,point) in zip(rows,points)
            yields_valid(row)
            algebra=closed_amplitude(row,point);q=point["kinematics"]["pT_GeV"]/point["kinematics"]["z1"]
            require(row["angular_factor"]==point["source_angular"]["value"] && row["algebra_reference"]==algebra &&
                row["fixed_model_qT_GeV"]==q && row["fixed_model_qT_below_source_upper_cut"]===(q<3.5),"changed fixed-mean angular/kinematic algebra")
        end
    end
    levels=r["refinements"]
    require(policy["minimum_refinement_level"]<=length(levels)<=length(policy["refinements"]),"invalid declared refinement budget")
    for (i,level) in enumerate(levels)
        require(level["level"]==i && level["setting"]==policy["refinements"][i] && level["bound_inverse_GeV"]==24,"changed frozen precision setting")
        rows_valid(level["points"])
        error=i==1 ? nothing : scaled(packed(level["points"]),packed(levels[i-1]["points"]),floor)
        target=error!==nothing && error<=policy["physical_value_max_scaled"];accepted=i>=policy["minimum_refinement_level"] && target
        require(level["previous_level_max_scaled"]==error && level["target_met"]===target && level["accepted_refinement"]===accepted,"changed refinement arithmetic/decision")
        require(!accepted || i==length(levels),"continued past accepted refinement")
    end
    boolean("all eighteen point joint refinements accepted",last(levels)["accepted_refinement"])
    tails=r["tails"];require([t["bound_inverse_GeV"] for t in tails]==policy["radial_bounds_inverse_GeV"],"changed radial extension bounds")
    for t in tails;rows_valid(t["points"]);end
    require(tails[2]["points"]==last(levels)["points"],"changed retained bound24 values")
    add("all eighteen point radial extension24to32",scaled(packed(tails[end]["points"]),packed(tails[end-1]["points"]),floor),policy["physical_value_max_scaled"])
    references=r["independent_references"];reuse=r["category_row_reuse_producer_attestations"]
    require(length(references)==length(reuse)==length(r["points"])==length(points),"incomplete final reference/reuse ledger")
    for (i,(point,final,reference,reused)) in enumerate(zip(points,r["points"],references,reuse))
        for key in keys(point);require(final[key]==point[key],"changed source/measured SIA field");end
        for (key,value) in tails[end]["points"][i];require(final[key]==value,"changed final SIA prediction");end
        require(reference["id"]==reused["id"]==point["id"] && reference["scope"]=="independent_fixed_momentum_scalar_assembly" &&
            reference["implementation_family"]=="generic_fixed_row_scalar_assembly_and_adaptive_P_reference" &&
            reference["momentum_interval"]===nothing && reference["native_rows_shared"]===true &&
            reference["initial_functions_shared"]===true && reference["independent_evolution"]===false &&
            reference["angular_value"]==point["source_angular"]["value"],"changed independent reference identity/scope")
        actual=Dict(k=>final[k] for k in ("pairs","categories","amplitude"))
        require(reference["actual"]==actual,"changed independent comparison operands")
        ref=reference["reference"];yields_valid(ref)
        require(scaled([ref["amplitude"]],[closed_amplitude(ref,point)],floor)<=policy["observable_algebra_max_scaled"],"changed independent normalized ratio")
        require(final["unweighted_observed_minus_predicted"]==point["measured"]["central_value"]-final["amplitude"],"changed unweighted difference")
        add("independent fixed-row scalar composition "*point["id"],scaled(packed_point(actual),packed_point(ref),floor),policy["physical_value_max_scaled"])
        add("source-model observable algebra "*point["id"],scaled([final["amplitude"]],[closed_amplitude(final,point)],floor),policy["observable_algebra_max_scaled"])
        boolean("positive charge/category denominators "*point["id"],true)
        require(reused["same_native_row_objects"] isa Bool,"missing category-reuse attestation")
        boolean("exact native category-row reuse "*point["id"],reused["same_native_row_objects"];attestation=true)
    end
    amplitudes=[p["amplitude"] for p in r["points"]];theta=Float64.(request.data["parameters"]["values"])
    flipped=copy(theta);flipped[[7,8]].*=-1
    sign=r["collins_sign_control"]
    require(sign["theta"]==theta && sign["changed_theta"]==flipped && sign["original"]==amplitudes && length(sign["flipped"])==length(amplitudes),"changed Collins sign probe")
    add("global Collins sign invariance",scaled(sign["flipped"],amplitudes,floor),policy["observable_algebra_max_scaled"])
    changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
    probe=r["transversity_value_control"]
    require(probe["theta"]==theta && probe["changed_theta"]==changed && probe["original"]==amplitudes && length(probe["changed"])==length(amplitudes) &&
        probe["scope"]=="fixed_coordinate_value_probe_not_global_derivative_certificate","changed transversity-coordinate probe")
    add("first six coordinate value independence",scaled(probe["changed"],amplitudes,floor),policy["observable_algebra_max_scaled"])
    boolean("source builds request stable during complete job",true;attestation=true)
    checks=r["checks"]
    require(length(checks)==length(expected)==185 && Set(c["name"] for c in checks)==Set(keys(expected)),"incomplete or duplicate diagnostic check ledger")
    for check in checks
        error,tolerance,_=expected[check["name"]]
        require(check["error"]==error && check["tolerance"]==tolerance && check["passed"]===(isfinite(error) && error<=tolerance) &&
            check["passed"]===true,"changed or failed numerical check/threshold")
    end
    require(r["final_bound_inverse_GeV"]==32 && r["unique_prepared_geometries"]==length(unique(historical_sia_geometry_key.(points)))==9,
        "changed final geometry/bound count")
    historical_sia_request_valid(request)
    Dict("status"=>r["status"],"review_status"=>"STORED_SIA_DIAGNOSTIC_ARITHMETIC_VERIFIED","scope"=>request.data["scope"],
        "accepted_prediction_diagnostic"=>true,"accepted_fitting_gate"=>false,
        "production_ready"=>false,"physical_extraction"=>false,"historical_fit_ready"=>false,"experimental_data_admitted"=>false,
        "fit_performed"=>false,"G2PERF_pass"=>false,"interval_coverage_verified"=>false,
        "request_sha256"=>request.sha256,"report_sha256"=>digest(reportfile),"source_fingerprint"=>request.data["source_fingerprint"],
        "source_bindings_verified"=>length(bindings),"stored_comparison_errors_recomputed"=>count(!last(v) for v in values(expected)),
        "producer_object_reuse_and_temporal_attestations"=>count(last(v) for v in values(expected)),
        "native_evolution_or_prediction_calls_replayed"=>false,"Fourier_quadrature_reexecuted"=>false,
        "raw_primary_tables_independently_reparsed_in_Julia"=>false,"current_build_source_and_binary_hashes_reverified"=>true,
        "source_preparation_kind"=>"Python source-replay attestation with current input hashes","auditor_sha256"=>digest(@__FILE__))
end
end
