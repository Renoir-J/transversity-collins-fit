# Pinned-input observables with explicit reference scale extrapolation.
# No experiment, fit, acceptance correction or production admission.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, QuadGK, Statistics
include("ObservableChecks.jl")

"Integrate original point structure functions, then apply the paper's ratios."
function direct_point_reference(context,theta;interval=nothing)
    process=context.leptonic.process;calls=Ref(0)
    contexts=process==:SIDIS ? [context] : [with_pion_charges(context,c...) for c in ((1,1),(1,-1),(-1,1),(-1,-1))]
    function point(p)
        calls[]+=1
        reduce(vcat,[[r.UU,r.Collins] for c in contexts for r in (polarized_structure(with_momentum(c,p),theta),)])
    end
    rows=interval===nothing ? point(context.leptonic.P) :
        quadgk(p->p.*point(p),interval...;rtol=1e-12,atol=1e-14)[1]
    if process==:SIDIS
        ratio=rows[2]/rows[1]
        values=[ratio,(2(1-OBS_Y)/(1+(1-OBS_Y)^2))*ratio]
    else
        U,L=rows[4]+rows[6],rows[2]+rows[8]
        DU,DL=rows[3]+rows[5],rows[1]+rows[7]
        values=OBS_ANGULAR.*[U/DU-L/DL,U/DU-(U+L)/(DU+DL)]
    end
    (values=values,integrated_charge_rows=rows,adaptive_evaluations=calls[])
end

function pinned_observable_checks(nativefile,inputfile,progressfile)
    sources=observable_source_hashes(@__FILE__)
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    hashes=Dict(p=>bytes2hex(sha256(read(p))) for p in (nativefile,inputfile))
    contract=JSON.parsefile(OBS_CONTRACT);policy=contract["validation_policy"]
    parent=JSON.parsefile(joinpath(OBS_ROOT,contract["parent_contract"]))["validation_policy"]
    checks=Dict{String,Any}[]
    check(name,e,t)=observable_check!(checks,name,e,t)
    scaled(a,b)=observable_scaled(a,b,policy)
    model=KPSYCSS();theta=collect(kpsy_published_parameters().values)
    pdf=CT10Input(inputfile;policy=:reference);hel=load_dssv08(inputs["input_root"]);ff=load_dss14(inputs["input_root"])
    levels=Dict{String,Any}();fine=Dict{Symbol,Any}();coarse=Dict{Symbol,Any}()
    momentum=Dict{String,Any}();tails=Dict{String,Any}();derivatives=[];direct=[];resources=[]
    prepare(process,radial,setting)=process==:SIDIS ?
        prepare_polarized_sidis(pdf,hel,ff,nativefile,model,radial,2.4,.3,.1,.3,1;
            matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference) :
        prepare_polarized_sia(ff,nativefile,model,radial,110.,.3,.3,.5,1,-1;
            matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference)
    try
        for process in (:SIDIS,:SIA)
            label=string(process);levels[label]=Any[]
            for (level,setting) in enumerate(parent["refinements"])
                radial=RadialRule(parton_radial_edges(model,ff,inputs["input_root"],24.);order=setting["radial_order"])
                context=nothing;cold=@elapsed context=prepare(process,radial,setting)
                values=all_observable_values(context,theta)
                error=isempty(levels[label]) ? nothing : scaled(values,last(levels[label])["values"])
                accepted=level>=3 && error!==nothing && error<=policy["physical_value_max_scaled"]
                meta=context.provenance
                explicit=meta["scope"]=="reconstructed_theory" && !meta["production_ready"] && !meta["physical_extraction"] &&
                    meta["range_policy"]=="reference" && meta["initial_source_kind"]=="pinned_inputs" &&
                    meta["source_range_audits"]["DSS14_outside_strict_nodes"]>0
                process==:SIDIS && (explicit &= meta["source_range_audits"]["CT10_outside_strict_nodes"]>0 && meta["second_initial_source_kind"]=="pinned_inputs")
                check("$process level $level explicit pinned/reference scope",explicit ? 0. : 1.,0.)
                check("$process level $level complete actual radial scales",context.first_operator.Q2==css_mu_b2.(Ref(model),radial.b) ? 0. : 1.,0.)
                check("$process level $level finite signed observables",all(isfinite,values) ? 0. : 1.,0.)
                push!(levels[label],Dict("setting"=>setting,"radial_nodes"=>length(radial.b),"values"=>values,
                    "previous_level_max_error"=>error,"refinement_accepted"=>accepted,"provenance"=>meta,
                    "preparation_seconds_informational"=>cold))
                coarse[process]=get(coarse,process,context);fine[process]=context
                writejson(progressfile,Dict("status"=>"IN_PROGRESS","scope"=>"reconstructed_theory","production_ready"=>false,
                    "sources"=>sources,"refinements"=>levels,"build_sha256"=>hashes[nativefile],"input_build_sha256"=>hashes[inputfile]))
                println(JSON.json(Dict("completed_process"=>label,"level"=>level,"radial_nodes"=>length(radial.b),"error"=>error,"accepted"=>accepted)));flush(stdout)
                accepted && break
            end
            check("$process full pinned observable refinement admitted",last(levels[label])["refinement_accepted"] ? 0. : 1.,0.)
            check("$process adjacent matching/radial/native observable difference",last(levels[label])["previous_level_max_error"],policy["physical_value_max_scaled"])
            context=fine[process]
            for interval in (nothing,OBS_BINS...)
                pair=observable_pair(context;interval=interval)
                actual=observable_values(pair,theta)
                reference=direct_point_reference(context,theta;interval=interval)
                check("$process interval $interval independent point/adaptive-P composition",scaled(actual,reference.values),policy["analytic_integral_max_scaled"])
                check("$process interval $interval observable metadata preserves scope",all(o->o.provenance["scope"]=="reconstructed_theory" &&
                    o.provenance["production_ready"]===false && o.provenance["physical_extraction"]===false,pair) ? 0. : 1.,0.)
                push!(direct,Dict("process"=>label,"momentum_bounds_GeV"=>interval,"values"=>actual,"reference_values"=>reference.values,
                    "integrated_charge_rows"=>reference.integrated_charge_rows,"adaptive_evaluations"=>reference.adaptive_evaluations,
                    "reference_scope"=>"direct original point structure functions; independent adaptive P integration and manual ratios"))
            end
            momentum[label]=momentum_refinements(context,theta,policy,checks)
            tail=[];setting=last(levels[label])["setting"]
            for B in parent["radial_bounds"]
                tailcontext=B==24 ? context : prepare(process,
                    RadialRule(parton_radial_edges(model,ff,inputs["input_root"],B);order=setting["radial_order"]),setting)
                push!(tail,Dict("B_inverse_GeV"=>B,"values"=>all_observable_values(tailcontext,theta)))
            end
            tails[label]=tail
            check("$process pinned finite-bound observable extension 24 to 32",scaled(tail[end]["values"],tail[end-1]["values"]),policy["physical_value_max_scaled"])
            append!(derivatives,observable_derivatives(coarse[process],theta,policy,checks))
            # Identical point values, parameters, source rows and output pair.
            # Resource observations never certify a production performance gate.
            pair=observable_pair(context)
            observable_values(pair,theta);direct_point_reference(context,theta)
            projected=[@elapsed observable_values(pair,theta) for _ in 1:5]
            original=[@elapsed direct_point_reference(context,theta) for _ in 1:5]
            rowbytes=context.first_operator.coefficient_payload_bytes+
                (context.second_operator===nothing ? 0 : context.second_operator.coefficient_payload_bytes)
            push!(resources,Dict("process"=>label,"scope"=>"identical single-point two-observable workload",
                "production_ready"=>false,"G2PERF_pass"=>false,"peak_RSS_measured"=>false,
                "radial_nodes"=>length(context.leptonic.radial.b),"row_coefficient_payload_bytes"=>rowbytes,
                "observable_pair_seconds"=>projected,"original_point_composition_seconds"=>original,
                "median_observable_pair_seconds"=>median(projected),"median_original_point_composition_seconds"=>median(original)))
        end
        check("source/build stability",observable_source_hashes(@__FILE__)==sources && all(bytes2hex(sha256(read(p)))==h for (p,h) in hashes) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs ? 0. : 1.,0.)
        observable_result(checks,sources,contract;build=native,build_sha256=hashes[nativefile],input_build_sha256=hashes[inputfile],
            refinements=levels,tails=tails,momentum_refinements=momentum,direct_point_references=direct,derivative_checks=derivatives,
            resource_pilots=resources,parameters=theta,
            input_scope="strict pinned CT10/DSSV08/DSS14 initial samples; explicit reference mu_b scale extrapolation; no experimental admission")
    finally
        close(pdf)
    end
end

nativefile=isempty(ARGS) ? joinpath(OBS_ROOT,".runtime/hoppet/build.json") : ARGS[1]
inputfile=length(ARGS)<2 ? joinpath(OBS_ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(OBS_ROOT,"runs/observables/inputs.json") : ARGS[3]
try
    report=pinned_observable_checks(nativefile,inputfile,output*".progress.json");writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
        "physical_extraction"=>false,"sources"=>observable_source_hashes(@__FILE__),"error"=>sprint(showerror,exception)))
    rethrow()
end
