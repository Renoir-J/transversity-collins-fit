using CollinsFit, JSON, SHA, Dates, LinearAlgebra, QuadGK
include("PolarizedReference.jl")
include("ObservableReference.jl")

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/historical-sia-point-contract.json")
const SCOPE="historical_reconstructed_sia_point_prediction_diagnostic"
const PAIR_CHARGES=((1,1),(1,-1),(-1,1),(-1,-1))
const PAIR_NAMES=("pp","pm","mp","mm")
file_digest(path)=bytes2hex(sha256(read(path)))
scaled(a,b,floor)=maximum(abs.(a.-b)./max.(abs.(b),floor))

function category_values(pairs)
    u=Dict(k=>pairs["pm"][k]+pairs["mp"][k] for k in ("UU","Collins"))
    l=Dict(k=>pairs["pp"][k]+pairs["mm"][k] for k in ("UU","Collins"))
    Dict("U"=>u,"L"=>l,"C"=>Dict(k=>u[k]+l[k] for k in ("UU","Collins")))
end

function packed_point(point)
    vcat([point["pairs"][pair][k] for pair in PAIR_NAMES for k in ("UU","Collins")],
         [point["categories"][cat][k] for cat in ("U","L","C") for k in ("UU","Collins")],
         [point["amplitude"]])
end
packed(points)=reduce(vcat,packed_point.(points))

function source_projection(points,policy)
    records=Dict{String,Any}[]
    for point in points,sign in (-1,1)
        angular=point["source_angular"]["value"];ratio=point["ratio"]
        value=setprecision(policy["source_projection_reference_bits"]) do
            f=BigFloat(angular);du=BigFloat(2);dl=BigFloat(1);cu=BigFloat(sign)/8;cl=-BigFloat(sign)/16
            dc,cc=du+dl,cu+cl;dr,cr=ratio=="UL" ? (dl,cl) : (dc,cc)
            # Derivative at epsilon=0 of (N_U/<N_U>)/(N_R/<N_R>).
            derivative(phi;harmonic=x->cos(2x))=begin
                nu,nr=one(f),one(f)
                dnu=f*cu/du*harmonic(phi);dnr=f*cr/dr*harmonic(phi)
                (dnu*nr-nu*dnr)/nr^2
            end
            twopi=2BigFloat(pi)
            projection(harmonic)=quadgk(phi->derivative(phi;harmonic=harmonic)*cos(2phi),zero(f),twopi;
                rtol=big"1e-24",atol=big"1e-30")[1]/BigFloat(pi)
            base=cu/du-cr/dr;expected=f*base
            other=ratio=="UL" ? f*(cu/du-cc/dc) : f*(cu/du-cl/dl)
            (projected=Float64(projection(x->cos(2x))),expected=Float64(expected),
             positive=all(d>f*abs(c) for (d,c) in ((du,cu),(dl,cl),(dc,cc))),
             counterfactuals=Dict("omitted_Ftheta"=>Float64(base),"repeated_Ftheta"=>Float64(f*expected),
                "factor2"=>Float64(2expected),"wrong_category"=>Float64(other),"reversed_sign"=>Float64(-expected),
                "cos_phi0"=>Float64(projection(cos)),"sin_2phi0"=>Float64(projection(x->sin(2x)))))
        end
        u=(UU=2.,Collins=sign/8);l=(UU=1.,Collins=-sign/16);c=(UU=3.,Collins=sign/16)
        a=SIAMeanOfRatio(angular;source="explicitly synthetic RF0 linearization fixture",averaging_measure="uniform mathematical phi0")
        implementation=CollinsFit.sia_double_ratio((U=u,L=l,C=c),Symbol(ratio),a)
        push!(records,Dict("scope"=>"explicitly_synthetic_source_definition_fixture","point_id"=>point["id"],
            "ratio"=>ratio,"angular_factor"=>angular,"sign"=>sign,"D_U"=>2.,"D_L"=>1.,"C_U"=>sign/8,"C_L"=>-sign/16,
            "epsilon_interval"=>[-1.,1.],"linearization_at_epsilon"=>0.,"positive_counts_proved_by_bound"=>value.positive,
            "projected"=>value.projected,"reference"=>value.expected,"implementation"=>implementation,
            "error"=>scaled([implementation],[value.projected],policy["relative_floor"]),
            "counterfactuals"=>value.counterfactuals,
            "all_counterfactuals_differ_from_correct"=>all(abs(x-value.projected)>policy["relative_floor"] for x in values(value.counterfactuals))))
    end
    records
end

function predictions(request,contexts,theta)
    model=historical_sia_point_model(request,contexts);values=predict_joint(model,theta)
    rows=Dict{String,Any}[]
    for (i,(point,observation)) in enumerate(zip(request.data["points"],model.observations))
        context=observation.plan.context
        pairs=Dict{String,Any}()
        for (name,charges) in zip(PAIR_NAMES,PAIR_CHARGES)
            v=polarized_structure(with_pion_charges(context,charges...),theta)
            pairs[name]=Dict("UU"=>v.UU,"Collins"=>v.Collins)
        end
        categories=category_values(pairs);u=categories["U"];right=categories[point["ratio"]=="UL" ? "L" : "C"]
        f=point["source_angular"]["value"];q=point["kinematics"]["pT_GeV"]/point["kinematics"]["z1"]
        push!(rows,Dict("id"=>point["id"],"pairs"=>pairs,"categories"=>categories,"amplitude"=>values[i],
            "angular_factor"=>f,"algebra_reference"=>f*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"]),
            "fixed_model_qT_GeV"=>q,"fixed_model_qT_below_source_upper_cut"=>q<3.5))
    end
    (model=model,points=rows)
end

function run_checks!(report,requestfile,nativefile,inputfile,progressfile)
    sources=CollinsFit.historical_sia_source_bindings();report["sources"]=sources
    contract=JSON.parsefile(CONTRACT);policy=contract["numerical_policy"];report["contract"]=contract
    report["contract_sha256"]=file_digest(CONTRACT)
    request=load_historical_sia_point_request(requestfile)
    for key in ("parameters","model_convention","source_preparation","limitations")
        report[key]=request.data[key]
    end
    report["request_points"]=request.data["points"];report["request_sha256"]=request.sha256
    report["source_fingerprint"]=request.data["source_fingerprint"]
    report["source_preparation_independently_replayed_in_Julia"]=false
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    any(basename(p)=="hoppet_bridge.f90" && h==file_digest(joinpath(ROOT,"fortran/hoppet_bridge.f90")) for (p,h) in native["source_hashes"]) || error("stale native bridge")
    any(basename(p)=="lhapdf.cpp" && h==file_digest(joinpath(ROOT,"bridges/lhapdf.cpp")) for (p,h) in inputs["source_hashes"]) || error("stale input bridge")
    dependencies=Dict(p=>file_digest(p) for p in (requestfile,nativefile,inputfile))
    report["dependency_sha256"]=dependencies;report["native_build"]=native;report["input_build"]=inputs
    report["active_parton_input"]=Dict("provider"=>"DSS14","CT10_and_DSSV"=>"build verified; not used by SIA predictions")
    checks=report["checks"];floor=policy["relative_floor"]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error) && error<=tolerance))
    boolean(name,value)=check(name,value ? 0. : 1.,0.)
    projections=source_projection(request.data["points"],policy);report["source_projection"]=projections
    for (i,row) in enumerate(projections)
        check("source RF0 independent linearized Fourier projection $i",scaled([row["projected"]],[row["reference"]],floor),policy["source_projection_max_scaled"])
        check("source angular convention application $i",row["error"],policy["source_projection_max_scaled"])
        boolean("source fixture positive/distinct $i",row["positive_counts_proved_by_bound"] && row["all_counterfactuals_differ_from_correct"])
    end
    ff=load_dss14(inputs["input_root"]);css=KPSYCSS();theta=Float64.(request.data["parameters"]["values"])
    function contexts_for(setting,B)
        contexts=Dict{Tuple,PolarizedLeptonicPlan}();resources=Dict{String,Any}[]
        radial=RadialRule(parton_radial_edges(css,ff,inputs["input_root"],B);order=setting["radial_order"])
        for point in request.data["points"]
            key=historical_sia_geometry_key(point);haskey(contexts,key) && continue
            k=point["kinematics"];context=nothing
            elapsed=@elapsed context=prepare_polarized_sia(ff,nativefile,css,radial,k["Q2_GeV2"],k["pT_GeV"],k["z1"],k["z2"],1,1;
                matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference)
            context.first_samples.input_kind==:pinned_inputs && context.first_samples.kind==:KPSY_Hhat3_homogeneous &&
                context.second_samples===nothing || error("source job requires pinned fragmentation-only initial samples")
            contexts[key]=context
            push!(resources,Dict("geometry"=>collect(key),"radial_nodes"=>length(radial.b),
                "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes,
                "cold_seconds_informational"=>elapsed,"context_provenance"=>context.provenance))
            println(JSON.json(Dict("stage"=>"historical SIA context prepared","geometry"=>collect(key),"bound"=>B,
                "setting"=>setting,"input_kind"=>"pinned_inputs","scope"=>SCOPE)));flush(stdout)
        end
        contexts,resources
    end
    levels=report["refinements"];accepted=false;fine_contexts=nothing;fine_prediction=nothing;setting=nothing
    for (level,current) in enumerate(policy["refinements"])
        contexts,resources=contexts_for(current,24);result=predictions(request,contexts,theta)
        error=isempty(levels) ? nothing : scaled(packed(result.points),packed(last(levels)["points"]),floor)
        accepted=level>=policy["minimum_refinement_level"] && error!==nothing && error<=policy["physical_value_max_scaled"]
        push!(levels,Dict("level"=>level,"setting"=>current,"bound_inverse_GeV"=>24,"points"=>result.points,
            "previous_level_max_scaled"=>error,"target_met"=>error!==nothing && error<=policy["physical_value_max_scaled"],
            "accepted_refinement"=>accepted,"resources"=>resources))
        fine_contexts,fine_prediction,setting=contexts,result,current
        writejson(progressfile,Dict("status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>SCOPE,"production_ready"=>false,
            "physical_extraction"=>false,"refinements"=>levels))
        println(JSON.json(Dict("stage"=>"historical SIA refinement","level"=>level,"error"=>error,"accepted"=>accepted)));flush(stdout)
        accepted && break
    end
    boolean("all eighteen point joint refinements accepted",accepted)
    tails=Dict{String,Any}[];final_contexts=nothing;final_prediction=nothing
    for B in policy["radial_bounds_inverse_GeV"]
        contexts,resources=B==24 ? (fine_contexts,last(levels)["resources"]) : contexts_for(setting,B)
        result=B==24 ? fine_prediction : predictions(request,contexts,theta)
        push!(tails,Dict("bound_inverse_GeV"=>B,"points"=>result.points,"resources"=>resources))
        final_contexts,final_prediction=contexts,result
    end
    report["tails"]=tails
    check("all eighteen point radial extension24to32",scaled(packed(tails[end]["points"]),packed(tails[end-1]["points"]),floor),policy["physical_value_max_scaled"])
    references=Dict{String,Any}[];reuse=Dict{String,Any}[];model=final_prediction.model
    for (point,observation,row) in zip(request.data["points"],model.observations,final_prediction.points)
        context=observation.plan.context;f=point["source_angular"]["value"]
        independent=ObservableReference.prepare_generic_reference(context;angular=f)
        pol=independent.polarized_integrals(theta);den=independent.denominators
        pairs=Dict(name=>Dict("UU"=>den[j],"Collins"=>pol[j]) for (j,name) in enumerate(PAIR_NAMES))
        # Assemble reference categories explicitly, without the implementation's category helper.
        u=Dict("UU"=>den[2]+den[3],"Collins"=>pol[2]+pol[3]);l=Dict("UU"=>den[1]+den[4],"Collins"=>pol[1]+pol[4])
        categories=Dict("U"=>u,"L"=>l,"C"=>Dict(k=>u[k]+l[k] for k in ("UU","Collins")))
        reference=Dict("pairs"=>pairs,"categories"=>categories,"amplitude"=>independent.predict(theta)[point["ratio"]=="UL" ? 1 : 2])
        actual=Dict(k=>row[k] for k in ("pairs","categories","amplitude"))
        push!(references,Dict("id"=>point["id"],"actual"=>actual,"reference"=>reference,"angular_value"=>f,
            "scope"=>"independent_fixed_momentum_scalar_assembly","implementation_family"=>string(independent.scope),
            "momentum_interval"=>nothing,"native_rows_shared"=>true,"initial_functions_shared"=>true,"independent_evolution"=>false))
        check("independent fixed-row scalar composition "*point["id"],scaled(packed_point(actual),packed_point(reference),floor),policy["physical_value_max_scaled"])
        check("source-model observable algebra "*point["id"],scaled([row["amplitude"]],[row["algebra_reference"]],floor),policy["observable_algebra_max_scaled"])
        positive=all(row["pairs"][name]["UU"]>0 for name in PAIR_NAMES) && all(row["categories"][c]["UU"]>0 for c in ("U","L","C"))
        boolean("positive charge/category denominators "*point["id"],positive)
        base=final_contexts[historical_sia_geometry_key(point)]
        shared=context.first_operator===base.first_operator && context.second_operator===nothing &&
            all(with_pion_charges(context,charges...).first_operator===base.first_operator for charges in PAIR_CHARGES)
        push!(reuse,Dict("id"=>point["id"],"same_native_row_objects"=>shared))
        boolean("exact native category-row reuse "*point["id"],shared)
    end
    report["independent_references"]=references;report["category_row_reuse_producer_attestations"]=reuse
    values=[p["amplitude"] for p in final_prediction.points]
    flipped=copy(theta);flipped[[7,8]].*=-1
    changed=copy(theta);changed[1:6]=[0.,0.,1.,1.,1.,1.]
    flipped_values=predict_joint(model,flipped);changed_values=predict_joint(model,changed)
    report["collins_sign_control"]=Dict("theta"=>theta,"changed_theta"=>flipped,"original"=>values,"flipped"=>flipped_values)
    report["transversity_value_control"]=Dict("theta"=>theta,"changed_theta"=>changed,"original"=>values,"changed"=>changed_values,
        "scope"=>"fixed_coordinate_value_probe_not_global_derivative_certificate")
    check("global Collins sign invariance",scaled(flipped_values,values,floor),policy["observable_algebra_max_scaled"])
    check("first six coordinate value independence",scaled(changed_values,values,floor),policy["observable_algebra_max_scaled"])
    report["joint_model"]=joint_model_metadata(model)
    report["points"]=[merge(point,row,Dict("unweighted_observed_minus_predicted"=>point["measured"]["central_value"]-row["amplitude"]))
        for (point,row) in zip(request.data["points"],final_prediction.points)]
    report["final_bound_inverse_GeV"]=32;report["unique_prepared_geometries"]=length(final_contexts)
    historical_sia_request_valid(request)
    boolean("source builds request stable during complete job",CollinsFit.historical_sia_source_bindings()==sources &&
        all(file_digest(p)==h for (p,h) in dependencies) && CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs)
    report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    report
end

function main()
    length(ARGS)==4 || error("usage: predict_historical_sia_points.jl REQUEST NATIVE_BUILD PARTON_BUILD FRESH_OUTPUT")
    request,native,inputs,output=ARGS
    (ispath(output) || ispath(output*".progress.json")) && error("use a fresh output path; preserve earlier reports")
    BLAS.set_num_threads(1)
    report=Dict{String,Any}("status"=>"FAIL","scope"=>SCOPE,"theory_scope"=>"reconstructed_theory",
        "observation_scope"=>"historical_benchmark_candidate","input_kind"=>"source_preparation_attested_rows_and_pinned_parton_inputs",
        "production_ready"=>false,"physical_extraction"=>false,"historical_fit_ready"=>false,"experimental_data_admitted"=>false,
        "fit_performed"=>false,"likelihood_constructed"=>false,"G2PERF_pass"=>false,"interval_coverage_verified"=>false,
        "julia_version"=>string(VERSION),"checks"=>Dict{String,Any}[],"refinements"=>Dict{String,Any}[])
    try
        run_checks!(report,request,native,inputs,output*".progress.json")
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"scope"=>SCOPE,"checks"=>length(report["checks"]),
        "error"=>get(report,"error",nothing),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
