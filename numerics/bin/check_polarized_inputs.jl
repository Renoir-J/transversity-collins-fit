# Pinned-input complete structure functions with explicitly labeled reference
# scale extrapolation. No experiment, observable adapter or fitting admission.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, Statistics

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/polarized-process-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/hoppet-contract.json"),
    joinpath(ROOT,"config/hoppet-point-contract.json"),joinpath(ROOT,"config/hoppet-batch-contract.json"),
    joinpath(ROOT,"config/matched-rows-contract.json"),joinpath(ROOT,"config/css-contract.json"),
    joinpath(ROOT,"config/leptonic-contract.json"),joinpath(ROOT,"config/initial-contract.json"),
    joinpath(ROOT,"config/parton-contract.json"),joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"fortran/hoppet_bridge.f90"),joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)
rejected(f,message)=try f();false catch exception;occursin(message,sprint(showerror,exception)) end

function input_checks(nativefile,inputfile,progressfile)
    sources=source_hashes();native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    builds=Dict(p=>bytes2hex(sha256(read(p))) for p in (nativefile,inputfile))
    contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"];floor=policy["relative_floor"]
    scaled(a,b)=maximum(abs.(a.-b)./max.(abs.(b),floor))
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error)&&error<=tolerance))
    model=KPSYCSS();ff=load_dss14(inputs["input_root"]);hel=load_dssv08(inputs["input_root"])
    pdf=CT10Input(inputfile;policy=:reference);strict=CT10Input(inputfile)
    theta=collect(kpsy_published_parameters().values);momenta=[0.,.1,.3,.5,.7]
    levels=Dict{String,Any}();fine=Dict{Symbol,Any}();coarse=Dict{Symbol,Any}();resources=[]
    prepare(process,radial,setting)=process==:SIDIS ?
        prepare_polarized_sidis(pdf,hel,ff,nativefile,model,radial,2.4,0.,.1,.3,1;
            matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference) :
        prepare_polarized_sia(ff,nativefile,model,radial,110.,0.,.3,.5,1,-1;
            matching_order=setting["matching_order"],intervals=setting["intervals"],range_policy=:reference)
    function predictions(context,t)
        reduce(vcat,[[r.UU,r.Collins] for P in momenta for r in (polarized_structure(with_momentum(context,P),t),)])
    end
    try
        small=RadialRule([0.,.1,1.,4.];order=4)
        check("strict full SIDIS rejects unsupported scales",rejected(()->prepare_polarized_sidis(strict,hel,ff,nativefile,model,small,2.4,.3,.1,.3,1),"outside") ? 0. : 1.,0.)
        check("strict full SIA rejects unsupported scales",rejected(()->prepare_polarized_sia(ff,nativefile,model,small,110.,.3,.3,.5,1,-1),"outside strict range") ? 0. : 1.,0.)
        check("thrust frame is rejected",rejected(()->prepare_polarized_sia(ff,nativefile,model,small,110.,.3,.3,.5,1,-1;frame=:thrust),"hadron-2 frame") ? 0. : 1.,0.)
        for process in (:SIDIS,:SIA)
            label=string(process);levels[label]=Any[]
            for (level,setting) in enumerate(policy["refinements"])
                radial=RadialRule(parton_radial_edges(model,ff,inputs["input_root"],24.);order=setting["radial_order"])
                context=nothing;cold=@elapsed context=prepare(process,radial,setting)
                values=predictions(context,theta)
                previous=isempty(levels[label]) ? nothing : last(levels[label])["values"]
                error=previous===nothing ? nothing : scaled(values,previous)
                accepted=level>=3 && error!==nothing && error<=policy["physical_value_max_scaled"]
                meta=context.provenance
                metadata=meta["scope"]=="reconstructed_theory" && !meta["production_ready"] && !meta["physical_extraction"] &&
                    meta["range_policy"]=="reference" && meta["initial_source_kind"]=="pinned_inputs" &&
                    meta["source_range_audits"]["DSS14_outside_strict_nodes"]>0
                process==:SIDIS && (metadata &= meta["source_range_audits"]["CT10_outside_strict_nodes"]>0 && meta["second_initial_source_kind"]=="pinned_inputs")
                check("$label level $level explicit reference source scope",metadata ? 0. : 1.,0.)
                check("$label level $level actual radial scales",context.first_operator.Q2==css_mu_b2.(Ref(model),radial.b) ? 0. : 1.,0.)
                check("$label level $level finite signed structure functions",all(isfinite,values) ? 0. : 1.,0.)
                push!(levels[label],Dict("setting"=>setting,"radial_nodes"=>length(radial.b),"momenta_GeV"=>momenta,"values"=>values,
                    "previous_level_max_error"=>error,"refinement_accepted"=>accepted,"provenance"=>meta,
                    "preparation_seconds_informational"=>cold))
                writejson(progressfile,Dict("status"=>"IN_PROGRESS","scope"=>"reconstructed_theory","production_ready"=>false,"sources"=>sources,"refinements"=>levels))
                coarse[process]=get(coarse,process,context);fine[process]=context
                println(JSON.json(Dict("completed_process"=>label,"level"=>level,"radial_nodes"=>length(radial.b),"error"=>error,"accepted"=>accepted)));flush(stdout)
                accepted && break
            end
            check("$label full physical refinement",last(levels[label])["previous_level_max_error"],policy["physical_value_max_scaled"])
            context=fine[process]
            row_bytes=context.first_operator.coefficient_payload_bytes+
                (context.second_operator===nothing ? 0 : context.second_operator.coefficient_payload_bytes)
            fixed_context=with_momentum(context,.3)
            prediction(t)=polarized_structure(fixed_context,t).Collins
            prediction(theta);warm=[@elapsed prediction(theta) for _ in 1:5]
            push!(resources,Dict("process"=>label,"scope"=>"single_kinematic_point_reconstructed_structure_function",
                "production_ready"=>false,"G2PERF_pass"=>false,"row_coefficient_payload_bytes"=>row_bytes,
                "radial_nodes"=>length(context.leptonic.radial.b),"warm_prediction_seconds"=>warm,
                "median_warm_prediction_seconds"=>median(warm)))
            check("$label P=0 exact Collins zero",abs(polarized_structure(context,theta).Collins),0.)
            # Check charge reuse against fresh provider matching in the declared
            # reference policy; this is distinct from the analytic Mellin oracle.
            setting=last(levels[label])["setting"]
            if process==:SIDIS
                minus=with_pion_charges(context,-1)
                direct=prepare_unpolarized_sidis(pdf,ff,model,context.leptonic.radial,2.4,0.,.1,.3,-1;
                    matching_order=setting["matching_order"],range_policy=:reference)
                check("SIDIS charge reuse versus fresh unpolarized matching",scaled(minus.unpolarized,leptonic_unpolarized(direct.plan,direct.first,direct.second)),policy["algebra_max_scaled"])
                neutron=polarized_leptonic_plan(context.leptonic,context.unpolarized_first,context.unpolarized_second,
                    context.first_operator,context.first_samples;second_operator=context.second_operator,second_samples=context.second_samples,
                    target=NucleonTarget(:neutron),charges=context.charges,provenance=context.provenance)
                mixture=polarized_leptonic_plan(context.leptonic,context.unpolarized_first,context.unpolarized_second,
                    context.first_operator,context.first_samples;second_operator=context.second_operator,second_samples=context.second_samples,
                    target=NucleonTarget("algebra test mixture",(.6,.4),(.5,-.2);source="validation fixture; not a nuclear prescription"),charges=context.charges,provenance=context.provenance)
                p,n,m=[polarized_structure(with_momentum(c,.3),theta) for c in (context,neutron,mixture)]
                check("pinned target denominator impulse algebra",scaled(m.UU,.6p.UU+.4n.UU),policy["algebra_max_scaled"])
                check("pinned target numerator separate spin weights",scaled(m.Collins,.5p.Collins-.2n.Collins),policy["algebra_max_scaled"])
            else
                same=with_pion_charges(context,-1,1)
                direct=prepare_unpolarized_sia(ff,model,context.leptonic.radial,110.,0.,.3,.5,-1,1;
                    matching_order=setting["matching_order"],range_policy=:reference)
                check("SIA charge reuse versus fresh unpolarized matching",scaled(same.unpolarized,leptonic_unpolarized(direct.plan,direct.first,direct.second)),policy["algebra_max_scaled"])
                neutral=polarized_structure(with_momentum(with_pion_charges(context,0,0),.3),theta).Collins
                charged=sum(polarized_structure(with_momentum(with_pion_charges(context,c1,c2),.3),theta).Collins for c1 in (-1,1),c2 in (-1,1))/4
                check("pinned SIA neutral bilinear charge average",scaled(neutral,charged),policy["algebra_max_scaled"])
            end
        end
        tails=Dict{String,Any}()
        for process in (:SIDIS,:SIA)
            label=string(process);setting=last(levels[label])["setting"];tail=[]
            for B in policy["radial_bounds"]
                context=B==24 ? fine[process] : prepare(process,
                    RadialRule(parton_radial_edges(model,ff,inputs["input_root"],B);order=setting["radial_order"]),setting)
                push!(tail,Dict("B_inverse_GeV"=>B,"values"=>predictions(context,theta)))
            end
            tails[label]=tail
            check("$label pinned finite-bound extension 24 to 32",scaled(tail[end]["values"],tail[end-1]["values"]),policy["physical_value_max_scaled"])
        end
        derivatives=[]
        for process in (:SIDIS,:SIA)
            context=with_momentum(coarse[process],.3)
            prediction(t)=polarized_structure(context,t).Collins
            gradient(t)=ForwardDiff.gradient(prediction,t)
            d=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,.001]
            J=gradient(theta);hv=ForwardDiff.jacobian(s->gradient(theta+s[1]*d),[0.])[:,1]
            check("$process full published-boundary gradient/HVP finite",all(isfinite,J)&&all(isfinite,hv) ? 0. : 1.,0.)
            for step in policy["finite_difference_steps"]
                fd,fdh=setprecision(policy["finite_difference_reference_bits"]) do
                    t=BigFloat.(theta);h=BigFloat(step);direction=BigFloat.(d)
                    (Float64((-3prediction(t)+4prediction(t+h*direction)-prediction(t+2h*direction))/(2h)),
                     Float64.((-3gradient(t)+4gradient(t+h*direction)-gradient(t+2h*direction))/(2h)))
                end
                e1,e2=scaled(fd,dot(J,d)),scaled(fdh,hv)
                check("$process pinned full 13-parameter inward JVP step $step",e1,policy["derivative_max_scaled"])
                check("$process pinned full 13-parameter inward HVP step $step",e2,policy["derivative_max_scaled"])
                push!(derivatives,Dict("process"=>string(process),"step"=>step,"JVP_error"=>e1,"HVP_error"=>e2,
                    "radial_order"=>context.leptonic.radial.order,"radial_nodes"=>length(context.leptonic.radial.b),
                    "intervals"=>context.first_operator.intervals,"reference_bits"=>policy["finite_difference_reference_bits"]))
            end
        end
        check("source/build stability",source_hashes()==sources && all(bytes2hex(sha256(read(p)))==h for (p,h) in builds) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs ? 0. : 1.,0.)
        Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
            "physical_extraction"=>false,"contract_id"=>contract["id"],"checks"=>checks,"sources"=>sources,"build_sha256"=>builds[nativefile],
            "build"=>native,"input_build_sha256"=>builds[inputfile],"printed_parameters"=>theta,"refinements"=>levels,
            "tails"=>tails,"derivative_checks"=>derivatives,"resource_pilots"=>resources,"julia_version"=>string(VERSION),
            "completed_at"=>string(now(UTC)),"input_scope"=>"strict pinned initial samples; explicit reference mu_b scale extrapolation; no experimental fit",
            "limitations"=>contract["limitations"])
    finally
        close(pdf);close(strict)
    end
end

nativefile=isempty(ARGS) ? joinpath(ROOT,".runtime/hoppet/build.json") : ARGS[1]
inputfile=length(ARGS)<2 ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(ROOT,"runs/polarized-inputs/validation.json") : ARGS[3]
try
    report=input_checks(nativefile,inputfile,output*".progress.json");writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"error"=>sprint(showerror,exception)))
    rethrow()
end
