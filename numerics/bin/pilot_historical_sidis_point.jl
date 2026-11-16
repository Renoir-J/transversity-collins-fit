using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff
include("PolarizedReference.jl")
include("ObservableReference.jl")

const PILOT_ROOT=normpath(joinpath(@__DIR__,"../.."))
pilot_digest(path)=bytes2hex(sha256(read(path)))

function pilot_preflight!(report,requestfile,case_id,nativefile,inputfile)
    policy=historical_sidis_pilot_policy();request=load_historical_sidis_point_request(requestfile)
    point=sidis_pilot_case(request,policy,case_id);sources=sidis_pilot_current_sources(PILOT_ROOT)
    sidis_pilot_source_fingerprint(sources)==request.data["source_fingerprint"] || error("stale full SIDIS source fingerprint")
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    any(basename(path)=="hoppet_bridge.f90" && digest==pilot_digest(joinpath(PILOT_ROOT,"fortran/hoppet_bridge.f90"))
        for (path,digest) in native["source_hashes"]) || error("native bridge differs from pilot source")
    any(basename(path)=="lhapdf.cpp" && digest==pilot_digest(joinpath(PILOT_ROOT,"bridges/lhapdf.cpp"))
        for (path,digest) in inputs["source_hashes"]) || error("parton bridge differs from pilot source")
    report["sources"]=sources;report["source_fingerprint"]=request.data["source_fingerprint"]
    report["policy_sha256"]=pilot_digest(CollinsFit.HISTORICAL_SIDIS_PILOT_POLICY)
    report["parent_contract_sha256"]=request.data["contract_sha256"]
    report["request_sha256"]=request.sha256;report["request_points"]=132;report["point"]=deepcopy(point)
    report["case_id"]=case_id;report["parameters"]=deepcopy(request.data["parameters"])
    report["dependency_sha256"]=Dict(p=>pilot_digest(p) for p in (requestfile,nativefile,inputfile))
    report["native_build"]=native;report["input_build"]=inputs
    report["native_build_receipt"]=nativefile;report["parton_build_receipt"]=inputfile
    report["numerical_setting"]=deepcopy(policy["numerics"])
    report["limitations"]=vcat(policy["limitations"],request.data["limitations"])
    (;policy,request,point,native,inputs)
end

function run_sidis_pilot!(report,requestfile,case_id,nativefile,inputfile,progressfile)
    prepared=pilot_preflight!(report,requestfile,case_id,nativefile,inputfile)
    (;policy,request,point,inputs)=prepared;k=point["kinematics"];settings=policy["numerics"]["setting"]
    function progress(stage)
        report["current_stage"]=stage
        partial=deepcopy(report);partial["status"]="IN_PROGRESS_NOT_A_GATE";writejson(progressfile,partial)
        println(JSON.json(Dict("scope"=>report["scope"],"case_id"=>case_id,"stage"=>stage)));flush(stdout)
    end
    progress("source bindings verified; preparing one fixed geometry")
    pdf=nothing
    try
        pdf=CT10Input(inputfile;policy=:reference);ff=load_dss14(inputs["input_root"]);hel=load_dssv08(inputs["input_root"])
        css=KPSYCSS();bound=policy["numerics"]["radial_bound_inverse_GeV"]
        radial=RadialRule(parton_radial_edges(css,ff,inputs["input_root"],bound);order=settings["radial_order"])
        context,seconds=sidis_pilot_measure(()->prepare_polarized_sidis(pdf,hel,ff,nativefile,css,radial,
            k["Q2_GeV2"],k["pT_GeV"],k["x"],k["z"],1;matching_order=settings["matching_order"],
            intervals=settings["intervals"],range_policy=:reference,target=historical_sidis_target(point["target"]),
            construction=:batched_probes))
        context.first_samples.input_kind==context.second_samples.input_kind==:pinned_inputs || error("pilot requires pinned initial samples")
        report["resources"]=Dict("radial_nodes"=>length(radial.b),
            "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes+context.second_operator.coefficient_payload_bytes,
            "cold_preparation_seconds_informational"=>seconds,"context_provenance"=>deepcopy(context.provenance))
        progress("fixed geometry prepared; evaluating value and unvalidated all13 AD samples")
        observable=historical_sidis_point_observable(request,point["id"],context)
        theta=Float64.(request.data["parameters"]["values"])
        parent=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT);direction=Float64.(parent["derivative_policy"]["direction"])
        value(t)=predict_observable(observable,t)
        evaluations=sidis_pilot_evaluations(value,theta,direction,policy["numerics"]["derivative_evaluations"]["warm_repetitions"])
        charged=observable.context;structure=polarized_structure(charged,theta)
        report["values"]=Dict("UU"=>structure.UU,"FUT"=>structure.Collins,"amplitude"=>evaluations.amplitude)
        report["source_factor"]=sidis_pilot_source_factor(point)
        report["AD_samples"]=evaluations.samples;report["evaluation_times"]=evaluations.times
        progress("AD samples measured; assembling separate fixed-row scalar value reference")
        independent,reference_seconds=sidis_pilot_measure(()->ObservableReference.prepare_generic_reference(charged;y=k["y"]))
        ref=independent.predict(theta)
        reference=Dict("UU"=>independent.denominators[1],"FUT"=>independent.polarized_integrals(theta)[1],
            "amplitude"=>sidis_pilot_source_factor(point)*ref[1])
        report["independent_reference"]=reference
        report["reference_scope"]=Dict("family"=>independent.scope,"native_rows_shared"=>true,
            "initial_functions_shared"=>true,"independent_evolution"=>false,"derivative_reference_executed"=>false,
            "preparation_seconds_informational"=>reference_seconds)
        checks=sidis_pilot_value_checks(report["values"],reference,point,policy);report["checks"]=checks
        finite=evaluations.finite;repeatable=evaluations.repeatable
        report["finite_AD_samples"]=finite;report["repeatable_samples"]=repeatable
        historical_sidis_request_valid(request)
        sidis_pilot_current_sources(PILOT_ROOT)==report["sources"] || error("pilot sources changed during execution")
        all(pilot_digest(p)==h for (p,h) in report["dependency_sha256"]) || error("pilot dependency bytes changed")
        report["status"]=all(c["passed"] for c in checks) && finite && repeatable ?
            "PILOT_COMPLETED_NOT_DIAGNOSTIC_GATE" : "FAIL"
    finally
        pdf===nothing || close(pdf)
    end
    report
end

function main()
    length(ARGS)==5 || error("usage: pilot_historical_sidis_point.jl REQUEST CASE NATIVE_BUILD PARTON_BUILD FRESH_OUTPUT")
    request,case_id,native,inputs,output=ARGS
    (ispath(output) || ispath(output*".progress.json")) && error("use fresh SIDIS pilot paths; preserve earlier evidence")
    BLAS.set_num_threads(1)
    report=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_PILOT_SCOPE,
        "theory_scope"=>"reconstructed_theory","observation_scope"=>"historical_benchmark_candidate",
        (k=>false for k in CollinsFit.HISTORICAL_SIDIS_PILOT_FALSE_FLAGS)...,
        "julia_version"=>string(VERSION),"started_at"=>string(now(UTC)),"checks"=>Dict{String,Any}[],
        "full132_native_tail_reference_derivative_policy_completed"=>false)
    try
        run_sidis_pilot!(report,request,case_id,native,inputs,output*".progress.json")
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(output,report)
    end
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"case_id"=>case_id,
        "error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]=="PILOT_COMPLETED_NOT_DIAGNOSTIC_GATE" || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
