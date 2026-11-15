using CollinsFit, JSON, SHA, Dates, LinearAlgebra, QuadGK
include("PolarizedReference.jl")
include("ObservableReference.jl")

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/historical-point-contract.json")
const SCOPE="historical_reconstructed_point_prediction_diagnostic"
const SOURCES=vcat([CONTRACT,joinpath(ROOT,"config/initial-contract.json"),
    joinpath(ROOT,"config/historical-sources.json"),joinpath(ROOT,"fortran/hoppet_bridge.f90"),
    joinpath(ROOT,"bridges/lhapdf.cpp"),joinpath(ROOT,"numerics/Project.toml"),
    joinpath(ROOT,"numerics/Manifest.toml"),joinpath(ROOT,"requirements.lock")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)),
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/bin");join=true)),
    filter(p->endswith(p,".py"),readdir(joinpath(ROOT,"src/collinsflow");join=true)),
    filter(p->endswith(p,"-contract.json"),readdir(joinpath(ROOT,"config");join=true)))
source_hashes()=Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in unique(SOURCES))
file_digest(path)=bytes2hex(sha256(read(path)))
scaled(a,b,floor)=maximum(abs.(a.-b)./max.(abs.(b),floor))

function source_projection(points,policy)
    records=Dict{String,Any}[]
    for y in unique(Float64[p["kinematics"]["y"] for p in points]),sign in (-1,1)
        value=setprecision(policy["source_projection_reference_bits"]) do
            yy=BigFloat(y);R=BigFloat(sign)/8;spin=BigFloat(1)/2;S=BigFloat(1)/16
            factor=2(1-yy)/(1+(1-yy)^2);twopi=2BigFloat(pi)
            asymmetry(h,s)=begin
                modulation=factor*R*sin(h+s)+S*sin(h-s)
                up,down=1+spin*modulation,1-spin*modulation
                (up-down)/(up+down)/spin
            end
            inner(s)=quadgk(h->asymmetry(h,s)*sin(h+s),zero(yy),twopi;rtol=big"1e-24",atol=big"1e-30")[1]
            coefficient=quadgk(inner,zero(yy),twopi;rtol=big"1e-24",atol=big"1e-30")[1]/(2BigFloat(pi)^2)
            expected=factor*R
            (projected=Float64(coefficient),expected=Float64(expected),
             positive=spin*(abs(factor*R)+abs(S))<1,
             counterfactuals=Dict("wrong_phase"=>Float64(-expected),"omitted_DNN"=>Float64(R),
                 "repeated_DNN"=>Float64(factor*expected),"factor2"=>Float64(2expected)))
        end
        convention=SIDISAsymmetryConvention(phase=:sin_phi_sum,normalization=:fixed_y_DNN_included,y=y,
            source="explicitly synthetic HallA Eq1 two-angle count fixture")
        implementation=CollinsFit.sidis_convention_metadata(convention)["DNN_factor"]*(sign/8)
        push!(records,Dict("scope"=>"explicitly_synthetic_source_definition_fixture","y"=>y,"R"=>sign/8,
            "spin"=>.5,"Sivers"=>1/16,"positive_counts_proved_by_bound"=>value.positive,
            "projected"=>value.projected,"reference"=>value.expected,"implementation_factor_times_R"=>implementation,
            "error"=>scaled([implementation],[value.projected],policy["relative_floor"]),
            "counterfactuals"=>value.counterfactuals,
            "all_counterfactuals_distinct"=>all(abs(x-value.projected)>policy["relative_floor"] for x in values(value.counterfactuals))))
    end
    records
end

function predictions(request,contexts,theta)
    model=historical_point_model(request,contexts);values=predict_joint(model,theta)
    points=Dict{String,Any}[]
    for (i,(point,observation)) in enumerate(zip(request.data["points"],model.observations))
        structure=polarized_structure(observation.plan.context,theta)
        y=point["kinematics"]["y"];DNN=2(1-y)/(1+(1-y)^2)
        push!(points,Dict("id"=>point["id"],"UU"=>structure.UU,"Collins"=>structure.Collins,
            "amplitude"=>values[i],"DNN_factor"=>DNN,"algebra_reference"=>DNN*structure.Collins/structure.UU))
    end
    (model=model,points=points)
end

function run_checks!(report,requestfile,nativefile,inputfile,progressfile)
    sources=source_hashes();report["sources"]=sources
    contract=JSON.parsefile(CONTRACT);policy=contract["numerical_policy"];report["contract"]=contract
    report["contract_sha256"]=file_digest(CONTRACT)
    request=load_historical_point_request(requestfile)
    report["request_sha256"]=request.sha256;report["source_fingerprint"]=request.data["source_fingerprint"]
    report["request_points"]=request.data["points"];report["parameters"]=request.data["parameters"]
    report["model_convention"]=request.data["model_convention"]
    report["source_preparation"]=request.data["source_preparation"]
    report["source_preparation_independently_replayed_in_Julia"]=false
    report["limitations"]=contract["limitations"]
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    any(basename(path)=="hoppet_bridge.f90" && digest==file_digest(joinpath(ROOT,"fortran/hoppet_bridge.f90"))
        for (path,digest) in native["source_hashes"]) || error("retained native bridge does not match current project source")
    any(basename(path)=="lhapdf.cpp" && digest==file_digest(joinpath(ROOT,"bridges/lhapdf.cpp"))
        for (path,digest) in inputs["source_hashes"]) || error("parton bridge does not match current project source")
    dependencies=Dict(path=>file_digest(path) for path in (requestfile,nativefile,inputfile))
    report["dependency_sha256"]=dependencies;report["native_build"]=native;report["input_build"]=inputs
    checks=report["checks"];floor=policy["relative_floor"]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>isfinite(error) && error<=tolerance))
    boolean(name,passed)=check(name,passed ? 0. : 1.,0.)
    projections=source_projection(request.data["points"],policy);report["source_projection"]=projections
    for (i,row) in enumerate(projections)
        check("source Eq1 independent two-angle projection $i",scaled([row["projected"]],[row["reference"]],floor),policy["source_projection_max_scaled"])
        check("source convention application $i",row["error"],policy["source_projection_max_scaled"])
        boolean("source fixture positive/distinct $i",row["positive_counts_proved_by_bound"] && row["all_counterfactuals_distinct"])
    end
    pdf=nothing
    try
        pdf=CT10Input(inputfile;policy=:reference);ff=load_dss14(inputs["input_root"]);hel=load_dssv08(inputs["input_root"])
        css=KPSYCSS();theta=Float64.(request.data["parameters"]["values"])
        function contexts_for(setting,B)
            contexts=Dict{Tuple,PolarizedLeptonicPlan}();resources=Dict{String,Any}[]
            radial=RadialRule(parton_radial_edges(css,ff,inputs["input_root"],B);order=setting["radial_order"])
            for point in request.data["points"]
                key=historical_geometry_key(point);haskey(contexts,key) && continue
                k=point["kinematics"];context=nothing
                elapsed=@elapsed context=prepare_polarized_sidis(pdf,hel,ff,nativefile,css,radial,
                    k["Q2_GeV2"],k["pT_GeV"],k["x"],k["z"],1;
                    matching_order=setting["matching_order"],intervals=setting["intervals"],
                    range_policy=:reference,target=NucleonTarget(:neutron))
                context.first_samples.input_kind==:pinned_inputs && context.second_samples.input_kind==:pinned_inputs || error("real-source job requires pinned initial samples")
                contexts[key]=context
                push!(resources,Dict("geometry"=>collect(key),"radial_nodes"=>length(radial.b),
                    "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes+context.second_operator.coefficient_payload_bytes,
                    "cold_seconds_informational"=>elapsed,"context_provenance"=>context.provenance))
                println(JSON.json(Dict("stage"=>"historical context prepared","geometry"=>collect(key),"bound"=>B,
                    "setting"=>setting,"input_kind"=>"pinned_inputs","scope"=>SCOPE)));flush(stdout)
            end
            contexts,resources
        end
        packed(points)=reduce(vcat,[[row["UU"],row["Collins"],row["amplitude"]] for row in points])
        levels=report["refinements"];accepted=false;fine_contexts=nothing;fine_prediction=nothing;setting=nothing
        for (level,current) in enumerate(policy["refinements"])
            contexts,resources=contexts_for(current,policy["radial_refinement_bound_inverse_GeV"])
            result=predictions(request,contexts,theta)
            values=packed(result.points)
            error=isempty(levels) ? nothing : scaled(values,packed(last(levels)["points"]),floor)
            accepted=level>=policy["minimum_refinement_level"] && error!==nothing && error<=policy["physical_value_max_scaled"]
            push!(levels,Dict("level"=>level,"setting"=>current,"bound_inverse_GeV"=>24,
                "points"=>result.points,"previous_level_max_scaled"=>error,
                "target_met"=>error!==nothing && error<=policy["physical_value_max_scaled"],
                "accepted_refinement"=>accepted,"resources"=>resources))
            fine_contexts,fine_prediction,setting=contexts,result,current
            writejson(progressfile,Dict("status"=>"IN_PROGRESS_NOT_A_GATE","scope"=>SCOPE,
                "production_ready"=>false,"physical_extraction"=>false,"refinements"=>levels))
            println(JSON.json(Dict("stage"=>"historical refinement","level"=>level,"error"=>error,"accepted"=>accepted)));flush(stdout)
            accepted && break
        end
        boolean("all eight point joint refinements accepted",accepted)
        tails=Dict{String,Any}[];final_contexts=nothing;final_prediction=nothing
        for B in policy["radial_bounds_inverse_GeV"]
            contexts,resources=B==24 ? (fine_contexts,last(levels)["resources"]) : contexts_for(setting,B)
            result=B==24 ? fine_prediction : predictions(request,contexts,theta)
            push!(tails,Dict("bound_inverse_GeV"=>B,"points"=>result.points,"resources"=>resources))
            final_contexts,final_prediction=contexts,result
        end
        report["tails"]=tails
        check("all eight point radial extension24to32",scaled(packed(tails[end]["points"]),packed(tails[end-1]["points"]),floor),policy["physical_value_max_scaled"])
        model=final_prediction.model;references=Dict{String,Any}[];charge_reuse=Dict{String,Any}[]
        for (i,(point,observation)) in enumerate(zip(request.data["points"],model.observations))
            row=final_prediction.points[i];context=observation.plan.context
            independent=ObservableReference.prepare_generic_reference(context;y=point["kinematics"]["y"])
            reference=independent.predict(theta)
            refvalues=[independent.denominators[1],independent.polarized_integrals(theta)[1],reference[2]]
            actual=[row["UU"],row["Collins"],row["amplitude"]]
            push!(references,Dict("id"=>point["id"],"actual"=>actual,"reference"=>refvalues,
                "reference_structure_ratio"=>reference[1],"scope"=>"independent_fixed_momentum_scalar_assembly",
                "implementation_family"=>string(independent.scope),"momentum_interval"=>nothing,
                "native_rows_shared"=>true,"initial_functions_shared"=>true,"independent_evolution"=>false))
            check("independent fixed-row scalar composition "*point["id"],scaled(actual,refvalues,floor),policy["physical_value_max_scaled"])
            check("source-model observable algebra "*point["id"],scaled([row["amplitude"]],[row["algebra_reference"]],floor),policy["observable_algebra_max_scaled"])
            base=final_contexts[historical_geometry_key(point)]
            shared=context.first_operator===base.first_operator && context.second_operator===base.second_operator
            push!(charge_reuse,Dict("id"=>point["id"],"same_native_row_objects"=>shared))
            boolean("exact native charge-row reuse "*point["id"],shared)
        end
        report["independent_references"]=references
        report["charge_row_reuse_producer_attestations"]=charge_reuse
        values=[row["amplitude"] for row in final_prediction.points]
        flipped=copy(theta);flipped[[1,2,7,8]].*=-1
        flipped_values=predict_joint(model,flipped)
        report["global_sign_control"]=Dict("original"=>values,"flipped"=>flipped_values)
        check("global transversity Collins sign invariance",scaled(flipped_values,values,floor),policy["observable_algebra_max_scaled"])
        report["joint_model"]=joint_model_metadata(model)
        report["points"]=[merge(point,row,Dict("unweighted_observed_minus_predicted"=>point["measured"]["central_value"]-row["amplitude"]))
            for (point,row) in zip(request.data["points"],final_prediction.points)]
        report["final_bound_inverse_GeV"]=last(policy["radial_bounds_inverse_GeV"])
        report["unique_prepared_geometries"]=length(final_contexts)
        historical_request_valid(request)
        boolean("source builds request stable during complete job",source_hashes()==sources &&
            all(file_digest(path)==digest for (path,digest) in dependencies) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs)
        report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    finally
        pdf===nothing || close(pdf)
    end
    report
end

function main()
    length(ARGS)==4 || error("usage: predict_historical_points.jl REQUEST NATIVE_BUILD PARTON_BUILD FRESH_OUTPUT")
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
