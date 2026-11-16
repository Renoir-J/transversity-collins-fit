module FixedContextStorageRecords
using CollinsFit, JSON, SHA, ForwardDiff, Dates

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const SCOPE="reconstructed_fixed_leptonic_storage_roundtrip"
const FLAGS=CollinsFit.FIXED_STORAGE_FALSE_FLAGS
digest(path)=open(sha256,path) |> bytes2hex

function bindings(sidisfile,siafile,nativefile,inputfile)
    contract=fixed_context_storage_contract()
    sidis=load_historical_sidis_point_request(sidisfile)
    sia=load_historical_sia_point_request(siafile)
    sources=sidis_pilot_current_sources(ROOT);fingerprint=sidis_pilot_source_fingerprint(sources)
    fingerprint==sidis.data["source_fingerprint"]==sia.data["source_fingerprint"] || error("stale storage source requests")
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    paths=Dict{String,String}()
    for request in (sidis,sia)
        for (p,h) in request.data["source_inputs"]
            haskey(paths,p) && paths[p]!=h && error("inconsistent shared source input")
            paths[p]=h
        end
        paths[request.path]=request.sha256
    end
    paths[abspath(nativefile)]=digest(nativefile);paths[abspath(inputfile)]=digest(inputfile)
    theta=Float64.(sidis.data["parameters"]["values"])
    theta==Float64.(sia.data["parameters"]["values"])==collect(kpsy_published_parameters().values) || error("changed published center")
    policy=JSON.parsefile(joinpath(ROOT,"config/historical-sidis-point-contract.json"))
    direction=Float64.(policy["derivative_policy"]["direction"])
    physics=Dict(p=>h for (p,h) in sources if startswith(p,"config/") && endswith(p,".json"))
    common=Dict("source_fingerprint"=>fingerprint,
        "physics_contract_bundle_sha256"=>sidis_pilot_source_fingerprint(physics),
        "native_build_sha256"=>digest(nativefile),"parton_build_sha256"=>digest(inputfile),
        "raw_source_binding_sha256"=>sidis_pilot_source_fingerprint(paths),
        "storage_contract_sha256"=>digest(CollinsFit.FIXED_STORAGE_CONTRACT))
    (;contract,sidis,sia,nativefile=abspath(nativefile),inputfile=abspath(inputfile),native,inputs,
      sources,fingerprint,paths,theta,direction,common)
end

function selected(c,case_id)
    case_id in c.contract["validation"]["physical_cases"] || error("case is outside prospective storage selection")
    process=startswith(case_id,"KPSY15_SIDIS_") ? :SIDIS : :SIA
    request=process==:SIDIS ? c.sidis : c.sia
    found=filter(p->p["id"]==case_id,request.data["points"])
    length(found)==1 || error("selected source row missing")
    point=only(found)
    points=if process==:SIDIS
        [point]
    else
        key=historical_sia_geometry_key(point)
        matched=filter(p->historical_sia_geometry_key(p)==key,request.data["points"])
        length(matched)==2 && Set(p["ratio"] for p in matched)==Set(("UL","UC")) || error("both source SIA ratios required")
        matched
    end
    setting=c.contract["validation"]["native_setting"]
    geometry=CollinsFit.numerical_snapshot((case_id,string(process),points,setting))
    bound=copy(c.common);bound["geometry_binding_sha256"]=geometry
    (;case_id,process,request,point,points,setting,bound)
end

function header(c,s)
    report=Dict{String,Any}("scope"=>SCOPE,"theory_scope"=>"reconstructed_theory",
        "input_kind"=>"pinned_inputs","case_id"=>s.case_id,"process"=>string(s.process),
        "source_fingerprint"=>c.fingerprint,"sources"=>c.sources,"input_hashes"=>c.paths,
        "storage_contract_sha256"=>digest(CollinsFit.FIXED_STORAGE_CONTRACT),"bindings"=>s.bound,
        "source_points"=>deepcopy(s.points),"native_setting"=>deepcopy(s.setting),
        "parameters"=>c.theta,"direction"=>c.direction,"julia_version"=>string(VERSION),
        "limitations"=>deepcopy(c.contract["limitations"]),
        "reconstruction_limitations"=>[
            "Original author fitter,modified HOPPET and unrounded states remain unlocated by bounded searches",
            "Experimental error composition,covariance and projection overlap are not admitted",
            "Homogeneous Collins evolution and reference continuation outside strict input domains remain reconstructed",
            "Stored identity is not an independent evolution,parameter-domain or interval-coverage benchmark"],
        "source_preparation_semantics"=>"Current Python source attestation and Julia source/header byte verification; no new raw-table replay",
        "raw_primary_tables_independently_replayed_in_Julia"=>false,
        "independent_evolution_or_derivative_accuracy_accepted"=>false,
        "full132_or_full72_diagnostic_accepted"=>false,"interval_coverage_verified"=>false)
    for f in FLAGS;report[f]=false;end
    report
end

function valid_header(report,c,s)
    expected=header(c,s)
    all(get(report,k,nothing)==v for (k,v) in expected) || error("stale or changed storage case header")
    all(report[f]===false for f in FLAGS) || error("storage case cannot self-admit")
end

function stable(c)
    historical_sidis_request_valid(c.sidis);historical_sia_request_valid(c.sia)
    sidis_pilot_current_sources(ROOT)==c.sources || error("source changed during storage case")
    all(isfile(p) && digest(p)==h for (p,h) in c.paths) || error("input changed during storage case")
    CollinsFit.verified_hoppet_build(c.nativefile)==c.native && CollinsFit.verify_parton_build(c.inputfile)==c.inputs || error("changed storage native/provider build")
end

function prepare(c,s)
    k=s.point["kinematics"];setting=s.setting;pdf=nothing
    try
        ff=load_dss14(c.inputs["input_root"]);css=KPSYCSS()
        radial=RadialRule(parton_radial_edges(css,ff,c.inputs["input_root"],setting["bound_inverse_GeV"]);order=setting["radial_order"])
        if s.process==:SIDIS
            pdf=CT10Input(c.inputfile;policy=:reference);hel=load_dssv08(c.inputs["input_root"])
            return prepare_polarized_sidis(pdf,hel,ff,c.nativefile,css,radial,
                k["Q2_GeV2"],k["pT_GeV"],k["x"],k["z"],1;
                matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference,
                target=historical_sidis_target(s.point["target"]),construction=:batched_probes)
        end
        prepare_polarized_sia(ff,c.nativefile,css,radial,k["Q2_GeV2"],k["pT_GeV"],k["z1"],k["z2"],1,1;
            matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference,construction=:batched_probes)
    finally
        pdf===nothing || close(pdf)
    end
end

function observable(c,s,point,context)
    if s.process==:SIDIS
        return historical_sidis_point_observable(c.sidis,point["id"],context)
    end
    k=point["kinematics"];p=context.leptonic
    p.process==:SIA && p.Q2==k["Q2_GeV2"] && p.P==k["pT_GeV"] && p.z1==k["z1"] && p.z2==k["z2"] || error("changed SIA storage geometry")
    a=point["source_angular"]
    angular=SIAMeanOfRatio(a["value"];source="source-bound BaBar RF0 storage check; request SHA256="*c.sia.sha256,
        averaging_measure=a["averaging_measure"])
    sia_observable(context;ratio=Symbol(point["ratio"]),angular=angular,
        source="conditional fixed-center storage identity only; no fitting admission")
end

bits(x::Float64)=bitstring(x)
function metrics(c,s,context)
    result=Dict{String,Any}[]
    for point in s.points
        obs=observable(c,s,point,context);value=t->predict_observable(obs,t)
        v=value(c.theta);derivatives=sidis_job_derivatives(value,c.theta,c.direction)
        v isa Float64 && isfinite(v) || error("finite Float64 storage value required")
        gradient=derivatives.gradient;hvp=derivatives.HVP
        s.process==:SIA && (!all(iszero,gradient[1:6]) || !all(iszero,hvp[1:6])) && error("SIA first-six columns changed")
        push!(result,Dict("id"=>point["id"],"value"=>v,"gradient"=>gradient,"HVP"=>hvp,
            "bits"=>Dict("value"=>bits(v),"gradient"=>bits.(gradient),"HVP"=>bits.(hvp))))
    end
    result
end

function metrics_valid(actual,expected,s)
    length(actual)==length(expected)==length(s.points) || error("missing storage observable")
    for (a,b,p) in zip(actual,expected,s.points)
        a["id"]==b["id"]==p["id"] || error("changed storage output order")
        for record in (a,b)
            length(record["gradient"])==length(record["HVP"])==13 || error("incomplete all13 storage derivatives")
            record["value"] isa Float64 && isfinite(record["value"]) &&
                all(v->v isa Float64 && isfinite(v),record["gradient"]) &&
                all(v->v isa Float64 && isfinite(v),record["HVP"]) || error("finite Float64 storage metrics required")
            record["bits"]["value"]==bits(Float64(record["value"])) &&
                record["bits"]["gradient"]==bits.(Float64.(record["gradient"])) &&
                record["bits"]["HVP"]==bits.(Float64.(record["HVP"])) || error("changed saved IEEE operands")
        end
        a["bits"]==b["bits"] || error("cold/reloaded value or full13 derivatives differ in bits")
    end
    true
end

end
