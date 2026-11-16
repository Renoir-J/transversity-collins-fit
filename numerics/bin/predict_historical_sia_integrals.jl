using CollinsFit, JSON, SHA, Dates, LinearAlgebra
module HistoricalSIAPointComponents
include("predict_historical_sia_points.jl")
end
include("HistoricalSIAIntegralReference.jl")
const HIR=HistoricalSIAIntegralReference
const HSP=HistoricalSIAPointComponents
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/historical-sia-integral-contract.json")
const SCOPE="historical_reconstructed_sia_integral_prediction_diagnostic"
const PAIRS=("pp","pm","mp","mm")
const CHARGES=((1,1),(1,-1),(-1,1),(-1,-1))
digest(p)=bytes2hex(sha256(read(p)))
scaled(a,b,floor)=maximum(abs.(a.-b)./max.(abs.(b),floor))
packed(rows)=reduce(vcat,HSP.packed_point.(rows))

function atomic_report(path,data)
    temp=path*".writing";writejson(temp,data);mv(temp,path;force=true)
end

"Absolute pair integrals using the parent's prepared P*dP Fourier weights."
function pair_integrals(observable,theta)
    context=observable.context;plan=context.leptonic;p=KPSYParameters(theta;gh=plan.model.gh)
    first,second=CollinsFit.evolved_polarized_columns(context,p)
    damping=exp.(-plan.fixed_exponent.-CollinsFit.collins_width(plan.model,p.values[13]).*plan.fragmentation_shape)
    Dict(name=>begin
        A,D=CollinsFit.polarized_flavor_legs(context,first,second;charges=charges)
        Dict("UU"=>CollinsFit.observable_unpolarized(context,observable.uu_weights,charges...),
             "Collins"=>CollinsFit.finite_transform(observable.collins_weights,CollinsFit.flavor_product(plan,A,D).*damping))
    end for (name,charges) in zip(PAIRS,CHARGES))
end

function predictions(request,contexts,theta,order)
    model=historical_sia_integral_model(request,contexts;order=order)
    amplitudes=predict_joint(model,theta);rows=Dict{String,Any}[];pairs_by_geometry=Dict{Tuple,Any}()
    for (i,(point,observation)) in enumerate(zip(request.data["points"],model.observations))
        key=historical_sia_integral_geometry_key(point)
        pairs=get!(pairs_by_geometry,key) do;pair_integrals(observation.plan,theta);end
        categories=HSP.category_values(pairs);u=categories["U"];right=categories[point["ratio"]=="UL" ? "L" : "C"]
        angular=point["source_angular"]["value"]
        push!(rows,Dict("id"=>point["id"],"pairs"=>pairs,"categories"=>categories,"amplitude"=>amplitudes[i],
            "angular_factor"=>angular,"algebra_reference"=>angular*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"]),
            "momentum_order"=>order,"qT_bounds_GeV"=>copy(request.data["momentum_window"]["qT_bounds_GeV"]),
            "mapped_P_bounds_GeV"=>[observation.plan.momentum.lower,observation.plan.momentum.upper]))
    end
    (model=model,points=rows)
end

function run_checks!(report,requestfile,nativefile,inputfile,output)
    contract=JSON.parsefile(CONTRACT);policy=contract["numerical_policy"];dp=contract["representative_derivative_policy"]
    report["contract"]=contract;report["contract_sha256"]=digest(CONTRACT)
    request=load_historical_sia_integral_request(requestfile)
    for key in ("parameters","model_convention","momentum_window","source_preparation","limitations")
        report[key]=request.data[key]
    end
    report["request_points"]=request.data["points"];report["request_sha256"]=request.sha256
    report["source_fingerprint"]=request.data["source_fingerprint"]
    report["source_preparation_independently_replayed_in_Julia"]=false
    report["raw_primary_upper_cut_reparsed_in_Julia"]=false
    sources=CollinsFit.historical_sia_integral_source_bindings();report["sources"]=sources
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    any(basename(p)=="hoppet_bridge.f90" && h==digest(joinpath(ROOT,"fortran/hoppet_bridge.f90")) for (p,h) in native["source_hashes"]) || error("stale native bridge")
    any(basename(p)=="lhapdf.cpp" && h==digest(joinpath(ROOT,"bridges/lhapdf.cpp")) for (p,h) in inputs["source_hashes"]) || error("stale parton bridge")
    dependencies=Dict(p=>digest(p) for p in (requestfile,nativefile,inputfile));report["dependency_sha256"]=dependencies
    report["native_build"]=native;report["input_build"]=inputs
    report["active_parton_input"]=Dict("provider"=>"DSS14","CT10_and_DSSV"=>"build verified; not used by SIA predictions")
    report["preparation_reuse"]=Dict("kind"=>"unchanged_per_geometry_preparation","new_query_cache_used"=>false,
        "parameter_responses_retained_across_calls"=>false)
    checks=report["checks"];floor=policy["relative_floor"]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error) && error<=tolerance))
    boolean(name,value)=check(name,value ? 0. : 1.,0.)
    function checkpoint(stage)
        report["current_stage"]=stage;atomic_report(output,report)
        atomic_report(output*".progress.json",Dict("status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>SCOPE,
            "production_ready"=>false,"physical_extraction"=>false,"current_stage"=>stage,
            "refinements_completed"=>count(r->haskey(r,"points"),report["refinements"]),
            "checks_completed"=>length(checks)))
    end
    report["stencil_formula_control"]=HIR.mathematical_controls(dp)
    boolean("inward stencil closed-function controls",report["stencil_formula_control"]["passed"])
    projections=HSP.source_projection(request.data["points"],policy);report["source_projection"]=projections
    for (i,row) in enumerate(projections)
        check("source RF0 independent linearized Fourier projection $i",scaled([row["projected"]],[row["reference"]],floor),policy["source_projection_max_scaled"])
        check("source angular convention application $i",row["error"],policy["source_projection_max_scaled"])
        boolean("source fixture positive/distinct $i",row["positive_counts_proved_by_bound"] && row["all_counterfactuals_differ_from_correct"])
    end
    ff=load_dss14(inputs["input_root"]);css=KPSYCSS();theta=Float64.(request.data["parameters"]["values"])
    function contexts_for(setting,B,resources,label)
        contexts=Dict{Tuple,PolarizedLeptonicPlan}()
        radial=RadialRule(parton_radial_edges(css,ff,inputs["input_root"],B);order=setting["radial_order"])
        for point in request.data["points"]
            key=historical_sia_integral_geometry_key(point);haskey(contexts,key) && continue
            checkpoint(Dict("stage"=>"preparing fixed geometry","label"=>label,"geometry"=>collect(key),"bound"=>B))
            k=point["kinematics"];context=nothing
            elapsed=@elapsed context=prepare_polarized_sia(ff,nativefile,css,radial,k["Q2_GeV2"],0.,k["z1"],k["z2"],1,1;
                matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference)
            context.first_samples.input_kind==:pinned_inputs && context.first_samples.kind==:KPSY_Hhat3_homogeneous &&
                context.second_samples===nothing || error("pinned fragmentation-only initial samples required")
            contexts[key]=context
            push!(resources,Dict("geometry"=>collect(key),"radial_nodes"=>length(radial.b),
                "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes,
                "cold_seconds_informational"=>elapsed,"context_provenance"=>context.provenance))
            println(JSON.json(Dict("stage"=>"historical SIA integral context prepared","label"=>label,
                "geometry"=>collect(key),"bound"=>B,"setting"=>setting,"scope"=>SCOPE)));flush(stdout)
            checkpoint(Dict("stage"=>"fixed geometry prepared","label"=>label,"geometries"=>length(contexts),"bound"=>B))
        end
        contexts
    end
    function momentum_refine!(entry,contexts,label)
        levels=entry["momentum_refinements"];result=nothing
        for order in policy["momentum_orders"]
            checkpoint(Dict("stage"=>"momentum refinement","label"=>label,"order"=>order))
            result=predictions(request,contexts,theta,order)
            error=isempty(levels) ? nothing : scaled(packed(result.points),packed(last(levels)["points"]),floor)
            accepted=order>=policy["minimum_momentum_order"] && error!==nothing && error<=policy["momentum_max_scaled"]
            push!(levels,Dict("order"=>order,"points"=>result.points,"previous_order_max_scaled"=>error,
                "target_met"=>error!==nothing && error<=policy["momentum_max_scaled"],"accepted"=>accepted))
            checkpoint(Dict("stage"=>"momentum level completed","label"=>label,"order"=>order,"accepted"=>accepted))
            accepted && break
        end
        boolean("momentum accepted "*label,last(levels)["accepted"])
        last(levels)["accepted"] || error("declared momentum refinement exhausted for "*label)
        entry["points"]=result.points;entry["momentum_order"]=last(levels)["order"]
        result
    end
    accepted=false;fine_contexts=nothing;fine_result=nothing;setting=nothing
    for (level,current) in enumerate(policy["refinements"])
        entry=Dict{String,Any}("level"=>level,"setting"=>current,"bound_inverse_GeV"=>24,"resources"=>Any[],"momentum_refinements"=>Any[])
        push!(report["refinements"],entry)
        contexts=contexts_for(current,24,entry["resources"],"native=$level")
        result=momentum_refine!(entry,contexts,"native=$level")
        error=level==1 ? nothing : scaled(packed(result.points),packed(report["refinements"][level-1]["points"]),floor)
        accepted=level>=policy["minimum_refinement_level"] && error!==nothing && error<=policy["physical_value_max_scaled"]
        merge!(entry,Dict("previous_level_max_scaled"=>error,"target_met"=>error!==nothing && error<=policy["physical_value_max_scaled"],"accepted_refinement"=>accepted))
        fine_contexts,fine_result,setting=contexts,result,current
        println(JSON.json(Dict("stage"=>"historical SIA integral native refinement","level"=>level,"error"=>error,"accepted"=>accepted)));flush(stdout)
        checkpoint(Dict("stage"=>"native level completed","level"=>level,"accepted"=>accepted))
        accepted && break
    end
    boolean("all seventy-two point joint refinements accepted",accepted)
    accepted || error("declared native refinement exhausted; later checks were not executed")
    tails=Any[];report["tails"]=tails;final_contexts=nothing;final_result=nothing
    for B in policy["radial_bounds_inverse_GeV"]
        entry=Dict{String,Any}("bound_inverse_GeV"=>B,"resources"=>Any[],"momentum_refinements"=>Any[])
        push!(tails,entry)
        if B==24
            for key in ("resources","momentum_refinements","points","momentum_order")
                entry[key]=last(report["refinements"])[key]
            end
            contexts,result=fine_contexts,fine_result
        else
            contexts=contexts_for(setting,B,entry["resources"],"tail B=$B")
            result=momentum_refine!(entry,contexts,"tail B=$B")
        end
        final_contexts,final_result=contexts,result
        checkpoint(Dict("stage"=>"finite radial bound completed","bound"=>B))
    end
    check("all seventy-two point radial extension24to32",scaled(packed(tails[3]["points"]),packed(tails[2]["points"]),floor),policy["physical_value_max_scaled"])
    references=Any[];reuse=Any[];report["independent_references"]=references;report["category_row_reuse_producer_attestations"]=reuse
    for (point,observation,row) in zip(request.data["points"],final_result.model.observations,final_result.points)
        checkpoint(Dict("stage"=>"independent adaptive qT reference","point_id"=>point["id"]))
        context=observation.plan.context;angular=point["source_angular"]["value"]
        reference=HSP.ObservableReference.prepare_generic_reference(context;momentum_interval=(0.,3.5),momentum_coordinate=:qT,angular=angular)
        pol=reference.polarized_integrals(theta);den=reference.denominators
        pairs=Dict(name=>Dict("UU"=>den[i],"Collins"=>pol[i]) for (i,name) in enumerate(PAIRS))
        u=Dict("UU"=>den[2]+den[3],"Collins"=>pol[2]+pol[3]);l=Dict("UU"=>den[1]+den[4],"Collins"=>pol[1]+pol[4])
        categories=Dict("U"=>u,"L"=>l,"C"=>Dict(k=>u[k]+l[k] for k in ("UU","Collins")))
        ref=Dict("pairs"=>pairs,"categories"=>categories,"amplitude"=>reference.predict(theta)[point["ratio"]=="UL" ? 1 : 2])
        actual=Dict(k=>row[k] for k in ("pairs","categories","amplitude"))
        push!(references,Dict("id"=>point["id"],"actual"=>actual,"reference"=>ref,"angular_value"=>angular,
            "scope"=>reference.scope,"qT_bounds_GeV"=>[0.,3.5],"native_rows_shared"=>true,
            "initial_functions_shared"=>true,"independent_evolution"=>false))
        check("independent adaptive qT composition "*point["id"],scaled(HSP.packed_point(actual),HSP.packed_point(ref),floor),policy["physical_value_max_scaled"])
        check("source-model observable algebra "*point["id"],scaled([row["amplitude"]],[row["algebra_reference"]],floor),policy["observable_algebra_max_scaled"])
        boolean("positive integrated denominators "*point["id"],all(row["pairs"][name]["UU"]>0 for name in PAIRS) && all(row["categories"][c]["UU"]>0 for c in ("U","L","C")))
        base=final_contexts[historical_sia_integral_geometry_key(point)]
        shared=context.first_operator===base.first_operator && context.second_operator===nothing
        push!(reuse,Dict("id"=>point["id"],"same_native_row_objects"=>shared))
        boolean("exact native category-row reuse "*point["id"],shared)
    end
    amplitudes=[r["amplitude"] for r in final_result.points]
    flipped=copy(theta);flipped[[7,8]].*=-1;changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
    flipped_amplitudes=predict_joint(final_result.model,flipped);changed_amplitudes=predict_joint(final_result.model,changed)
    report["collins_sign_control"]=Dict("theta"=>theta,"changed_theta"=>flipped,"original"=>amplitudes,"flipped"=>flipped_amplitudes)
    report["transversity_value_control"]=Dict("theta"=>theta,"changed_theta"=>changed,"original"=>amplitudes,"changed"=>changed_amplitudes,
        "scope"=>"fixed_coordinate_value_probe_not_global_derivative_certificate")
    check("global Collins sign invariance",scaled(flipped_amplitudes,amplitudes,floor),policy["observable_algebra_max_scaled"])
    check("first six coordinate value independence",scaled(changed_amplitudes,amplitudes,floor),policy["observable_algebra_max_scaled"])
    derivatives=Any[];report["representative_derivatives"]=derivatives
    for id in dp["source_geometry_ids"]
        point=only(filter(p->p["id"]==id,request.data["points"]));key=historical_sia_integral_geometry_key(point)
        indices=[only(findall(p->historical_sia_integral_geometry_key(p)==key && p["ratio"]==ratio,request.data["points"])) for ratio in ("UL","UC")]
        plans=[final_result.model.observations[i].plan for i in indices]
        angulars=[request.data["points"][i]["source_angular"]["value"] for i in indices]
        ref=HSP.ObservableReference.prepare_generic_reference(final_contexts[key];momentum_interval=(0.,3.5),momentum_coordinate=:qT,angular=1.)
        candidate=t->[predict_observable(plan,t) for plan in plans];reference=t->angulars.*ref.predict(t)
        progress=()->checkpoint(Dict("stage"=>"inward derivative reference","geometry_source_id"=>id))
        record=HIR.derivative_reference(candidate,reference,theta,dp;progress=progress)
        record["geometry_source_id"]=id;record["point_ids"]=[request.data["points"][i]["id"] for i in indices]
        record["angular_values"]=angulars;record["native_rows_shared"]=true;record["independent_evolution"]=false
        push!(derivatives,record)
        check(id*" independent all13 gradient",record["independent_gradient_error"],dp["derivative_max_scaled"])
        check(id*" independent all13 HVP",record["independent_HVP_error"],dp["HVP_max_scaled"])
        boolean(id*" first six derivative columns zero",record["first_six_gradient_columns_zero"] && record["first_six_HVP_columns_zero"])
        for (i,r) in enumerate(record["finite_differences"])
            check(id*" inward gradient level=$i",r["gradient_error"],dp["derivative_max_scaled"])
            check(id*" inward HVP level=$i",r["HVP_error"],dp["HVP_max_scaled"])
        end
        check(id*" adjacent inward gradient",record["adjacent_gradient_error"],dp["derivative_max_scaled"])
        check(id*" adjacent inward HVP",record["adjacent_HVP_error"],dp["HVP_max_scaled"])
    end
    report["joint_model"]=joint_model_metadata(final_result.model)
    report["points"]=[merge(p,r,Dict("unweighted_observed_minus_predicted"=>p["measured"]["central_value"]-r["amplitude"])) for (p,r) in zip(request.data["points"],final_result.points)]
    report["final_bound_inverse_GeV"]=32;report["unique_prepared_geometries"]=length(final_contexts)
    historical_sia_integral_request_valid(request)
    boolean("source builds request stable during complete job",CollinsFit.historical_sia_integral_source_bindings()==sources &&
        all(digest(p)==h for (p,h) in dependencies) && CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs)
    report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    report
end

function main(args=ARGS)
    length(args)==4 || error("usage: predict_historical_sia_integrals.jl REQUEST NATIVE_BUILD PARTON_BUILD FRESH_OUTPUT")
    request,native,inputs,output=args
    any(ispath,(output,output*".progress.json",output*".writing")) && error("use fresh output paths; preserve prior reports")
    BLAS.set_num_threads(1)
    report=Dict{String,Any}("status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>SCOPE,"theory_scope"=>"reconstructed_theory",
        "observation_scope"=>"historical_benchmark_candidate","input_kind"=>"source_preparation_attested_rows_and_pinned_parton_inputs",
        "production_ready"=>false,"physical_extraction"=>false,"historical_fit_ready"=>false,"experimental_data_admitted"=>false,
        "fit_performed"=>false,"likelihood_constructed"=>false,"G2PERF_pass"=>false,"interval_coverage_verified"=>false,
        "julia_version"=>string(VERSION),"checks"=>Dict{String,Any}[],"refinements"=>Any[])
    atomic_report(output,report)
    try
        run_checks!(report,request,native,inputs,output)
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));atomic_report(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"scope"=>SCOPE,"checks"=>length(report["checks"]),
        "error"=>get(report,"error",nothing),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report
end
if abspath(PROGRAM_FILE)==(@__FILE__);main()["status"]=="PASS" || exit(2);end
