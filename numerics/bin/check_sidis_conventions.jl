# Source phase/normalization checks, distinct from data or fit admission.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, QuadGK
include("ObservableChecks.jl")

const SIDIS_CONTRACT=joinpath(OBS_ROOT,"config/sidis-convention-contract.json")
const SIDIS_SOURCE="explicit same-axis SIDIS source convention diagnostic; no experimental admission"
const SIDIS_PHASES=(:sin_phi_sum,:sin_phi_sum,:sin_phi_sum_minus_pi,:sin_phi_sum_minus_pi)
const SIDIS_NORMS=(:DNN_removed,:fixed_y_DNN_included,:DNN_removed,:fixed_y_DNN_included)

function convention_sources()
    merge(observable_source_hashes(@__FILE__),Dict(relpath(SIDIS_CONTRACT,OBS_ROOT)=>bytes2hex(sha256(read(SIDIS_CONTRACT)))))
end

function convention_plans(context;momentum=nothing)
    [sidis_observable(context,SIDISAsymmetryConvention(phase=phase,normalization=normalization,
        y=normalization==:DNN_removed ? nothing : 2/5,source=SIDIS_SOURCE);momentum)
        for (phase,normalization) in zip(SIDIS_PHASES,SIDIS_NORMS)]
end

"Synthetic positive counts and harmonic projection, independent of the implementation factor."
function fourier_projection_checks(factors,policy,checks)
    rows=setprecision(256) do
        result=[];twopi=2BigFloat(pi);dnn=BigFloat(15)/17
        for ratio in (BigFloat(1)/8,-BigFloat(1)/8),spin in (BigFloat(0),BigFloat(pi)/7,-BigFloat(pi)/3)
            counts(phi)=1+dnn*ratio*sin(phi+spin)
            for i in 1:4
                shifted=i>=3;removed=isodd(i)
                harmonic(phi)=sin(phi+spin-(shifted ? BigFloat(pi) : BigFloat(0)))
                numerator=quadgk(phi->counts(phi)*harmonic(phi),BigFloat(0),twopi;rtol=BigFloat("1e-60"))[1]
                denominator=quadgk(phi->harmonic(phi)^2,BigFloat(0),twopi;rtol=BigFloat("1e-60"))[1]
                raw=numerator/denominator
                projection=removed ? raw/dnn : raw
                exact=(shifted ? -1 : 1)*(removed ? ratio : dnn*ratio)
                label="synthetic ratio $ratio spin $spin convention $i"
                observable_check!(checks,label*" positive counts",1-abs(dnn*ratio)>0 ? 0. : 1.,0.)
                observable_check!(checks,label*" exact Fourier coefficient",
                    observable_scaled([Float64(projection)],[Float64(exact)],policy),policy["algebra_max_scaled"])
                observable_check!(checks,label*" implementation factor from Fourier projection",
                    observable_scaled([factors[i]*Float64(ratio)],[Float64(projection)],policy),policy["algebra_max_scaled"])
                push!(result,Dict("scope"=>"synthetic_mathematical_projection","ratio"=>Float64(ratio),
                    "target_spin_azimuth"=>Float64(spin),"y"=>2/5,"DNN_reference"=>Float64(dnn),
                    "phase"=>string(SIDIS_PHASES[i]),"normalization"=>string(SIDIS_NORMS[i]),
                    "raw_projected_amplitude"=>Float64(raw),"normalized_projection"=>Float64(projection),
                    "exact_reference"=>Float64(exact),"implementation_factor"=>factors[i],
                    "reference_bits"=>256,"positive_count_lower_bound"=>Float64(1-abs(dnn*ratio))))
            end
        end
        result
    end
    rows
end

function sidis_convention_checks(nativefile,inputfile,progressfile)
    sources=convention_sources()
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    hashes=Dict(p=>bytes2hex(sha256(read(p))) for p in (nativefile,inputfile))
    contract=JSON.parsefile(SIDIS_CONTRACT);policy=JSON.parsefile(OBS_CONTRACT)["validation_policy"]
    checks=Dict{String,Any}[]
    check(name,error,tolerance)=observable_check!(checks,name,error,tolerance)
    scaled(a,b)=observable_scaled(a,b,policy)
    model=KPSYCSS();theta=collect(kpsy_published_parameters().values)
    direction=[.02,.03,.01,.02,.01,.01,.02,.02,.01,.02,.03,.03,.001]
    pdf=CT10Input(inputfile;policy=:reference);hel=load_dssv08(inputs["input_root"]);ff=load_dss14(inputs["input_root"])
    comparisons=[];derivatives=[];fourier=[]
    try
        radial=RadialRule(parton_radial_edges(model,ff,inputs["input_root"],24.);order=8)
        context=prepare_polarized_sidis(pdf,hel,ff,nativefile,model,radial,2.4,.3,.1,.3,1;
            matching_order=4,intervals=512,range_policy=:reference)
        metadata=context.provenance
        explicit=metadata["scope"]=="reconstructed_theory" && !metadata["production_ready"] && !metadata["physical_extraction"] &&
            metadata["range_policy"]=="reference" && metadata["initial_source_kind"]==metadata["second_initial_source_kind"]=="pinned_inputs" &&
            metadata["source_range_audits"]["CT10_outside_strict_nodes"]>0 && metadata["source_range_audits"]["DSS14_outside_strict_nodes"]>0
        check("explicit pinned inputs and reference scale extrapolation",explicit ? 0. : 1.,0.)
        check("complete actual radial evolution scales",context.first_operator.Q2==css_mu_b2.(Ref(model),radial.b) ? 0. : 1.,0.)
        fourier=fourier_projection_checks([o.angular for o in convention_plans(context)],policy,checks)
        for interval in (nothing,(.1,.3))
            momentum=interval===nothing ? nothing : MomentumRule(interval...;order=128)
            plans=convention_plans(context;momentum)
            prediction(t)=[predict_observable(o,t) for o in plans]
            generic=prepare_generic_reference(context;momentum_interval=interval,y=2/5)
            function reference(t)
                r=generic.predict(t)
                [r[1],r[2],-r[1],-r[2]]
            end
            values=prediction(theta);references=reference(theta)
            check("interval $interval independent scalar/adaptive-P values",scaled(values,references),policy["algebra_max_scaled"])
            check("interval $interval finite signed nonzero values",all(isfinite,values)&&all(!iszero,values) ? 0. : 1.,0.)
            J=ForwardDiff.jacobian(prediction,theta)
            hv=reshape(ForwardDiff.jacobian(s->vec(ForwardDiff.jacobian(prediction,theta+s[1]*direction)),[0.])[:,1],size(J))
            Jr=setprecision(policy["finite_difference_reference_bits"]) do
                Float64.(ForwardDiff.jacobian(reference,BigFloat.(theta)))
            end
            check("interval $interval independent full13 gradients",scaled(J,Jr),policy["derivative_max_scaled"])
            check("interval $interval finite full13 gradients and HVP",all(isfinite,J)&&all(isfinite,hv) ? 0. : 1.,0.)
            canonical=sidis_observable(context;definition=:structure_function_ratio,momentum,source=SIDIS_SOURCE)
            oldraw=sidis_observable(context;definition=:sin_phi_sum_amplitude,y=2/5,momentum,source=SIDIS_SOURCE)
            check("interval $interval legacy canonical and raw values unchanged",values[1]==predict_observable(canonical,theta)&&values[2]==predict_observable(oldraw,theta) ? 0. : 1.,0.)
            provenance=[deepcopy(o.provenance) for o in plans]
            valid=all(m->m["sidis_convention_contract_id"]==contract["id"] && m["scope"]=="reconstructed_theory" &&
                m["production_ready"]===false && m["physical_extraction"]===false && m["source"]==SIDIS_SOURCE,provenance)
            check("interval $interval source-labeled reconstructed metadata",valid ? 0. : 1.,0.)
            counterfactuals=Dict("correct_shifted_removed"=>values[3],"incorrect_unshifted_removed"=>values[1],
                "incorrect_second_DNN"=>(15/17)*values[3])
            check("interval $interval unshifted phase counterfactual differs",!isapprox(values[3],values[1];rtol=policy["algebra_max_scaled"],atol=0.) ? 0. : 1.,0.)
            check("interval $interval repeated DNN counterfactual differs",!isapprox(values[3],counterfactuals["incorrect_second_DNN"];rtol=policy["algebra_max_scaled"],atol=0.) ? 0. : 1.,0.)
            push!(comparisons,Dict("momentum_bounds_GeV"=>interval,"values"=>values,"reference_values"=>references,
                "actual_gradient"=>J,"reference_gradient"=>Jr,"metadata"=>provenance,"counterfactuals"=>counterfactuals,
                "reference_scope"=>generic.scope,"momentum_order"=>(interval===nothing ? nothing : 128)))
            for step in policy["finite_difference_steps"]
                fd,fdh=setprecision(policy["finite_difference_reference_bits"]) do
                    t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
                    gradient(x)=ForwardDiff.jacobian(reference,x)
                    (Float64.((-3reference(t)+4reference(t+h*d)-reference(t+2h*d))/(2h)),
                     Float64.((-3gradient(t)+4gradient(t+h*d)-gradient(t+2h*d))/(2h)))
                end
                jvp=J*direction;e1=scaled(fd,jvp);e2=scaled(fdh,hv)
                check("interval $interval inward full13 JVP step $step",e1,policy["derivative_max_scaled"])
                check("interval $interval inward full13 HVP step $step",e2,policy["derivative_max_scaled"])
                push!(derivatives,Dict("momentum_bounds_GeV"=>interval,"step"=>step,
                    "actual_JVP"=>jvp,"reference_JVP"=>fd,"actual_HVP"=>hv,"reference_HVP"=>fdh,
                    "JVP_error"=>e1,"HVP_error"=>e2,"direction"=>direction,"parameters"=>theta,
                    "reference_bits"=>policy["finite_difference_reference_bits"],"reference_scope"=>generic.scope,
                    "radial_order"=>radial.order,"radial_nodes"=>length(radial.b),"intervals"=>512,"matching_order"=>4,
                    "momentum_order"=>(interval===nothing ? nothing : 128)))
                println(JSON.json(Dict("completed_interval"=>interval,"derivative_step"=>step,"JVP_error"=>e1,"HVP_error"=>e2)));flush(stdout)
            end
            writejson(progressfile,Dict("status"=>"IN_PROGRESS","scope"=>"reconstructed_theory","production_ready"=>false,
                "sources"=>sources,"point_and_bin_comparisons"=>comparisons,"derivative_checks"=>derivatives))
        end
        check("source/build stability",convention_sources()==sources && all(bytes2hex(sha256(read(p)))==h for (p,h) in hashes) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs ? 0. : 1.,0.)
        Dict("status"=>(all(c["passed"] for c in checks) ? "PASS" : "FAIL"),"scope"=>"reconstructed_theory",
            "contract_id"=>contract["id"],"parent_contract_id"=>"KPSY-LEPTONIC-OBS-002", "production_ready"=>false,
            "physical_extraction"=>false,"G2PERF_pass"=>false,"checks"=>checks,"sources"=>sources,
            "build"=>native,"build_sha256"=>hashes[nativefile],"input_build_sha256"=>hashes[inputfile],
            "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)),"parameters"=>theta,"direction"=>direction,
            "synthetic_fourier_projections"=>fourier,"point_and_bin_comparisons"=>comparisons,"derivative_checks"=>derivatives,
            "context_provenance"=>context.provenance,"limitations"=>contract["limitations"],
            "input_scope"=>"Pinned CT10/DSSV08/DSS14; explicit reference scale extrapolation; shared fixed native rows. Coarse adapter-composition check; existing independent observable refinements remain required. No experimental data.")
    finally
        close(pdf)
    end
end

nativefile=isempty(ARGS) ? joinpath(OBS_ROOT,".runtime/hoppet/build.json") : ARGS[1]
inputfile=length(ARGS)<2 ? joinpath(OBS_ROOT,".runtime/parton-inputs/build.json") : ARGS[2]
output=length(ARGS)<3 ? joinpath(OBS_ROOT,"runs/observables/sidis-source.json") : ARGS[3]
try
    report=sidis_convention_checks(nativefile,inputfile,output*".progress.json");writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
        "physical_extraction"=>false,"sources"=>convention_sources(),"error"=>sprint(showerror,exception)))
    rethrow()
end
