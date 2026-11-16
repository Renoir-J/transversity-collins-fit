using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff
include("PolarizedReference.jl")
include("ObservableReference.jl")
include("HistoricalSIDISJobRecords.jl")
const SJ=HistoricalSIDISJobRecords

function run_sidis_point!(report,ctx,phase,index,run,progressfile)
    point=ctx.request.data["points"][index];k=point["kinematics"]
    theta=Float64.(ctx.request.data["parameters"]["values"])
    direction=Float64.(ctx.contract["derivative_policy"]["direction"])
    report["dependency_sha256"]=Dict(p=>SJ.digest(p) for p in SJ.point_dependencies(ctx,run,phase,index))
    report["source_factor"]=sidis_pilot_source_factor(point)
    report["native_construction"]="batched_probes"
    function progress(stage)
        report["current_stage"]=stage
        partial=deepcopy(report);partial["status"]="IN_PROGRESS_NOT_A_GATE";writejson(progressfile,partial)
        println(JSON.json(Dict("phase"=>phase,"point_id"=>point["id"],"stage"=>stage)));flush(stdout)
    end
    pdf=nothing
    try
        pdf=CT10Input(ctx.inputfile;policy=:reference)
        ff=load_dss14(ctx.inputs["input_root"]);hel=load_dssv08(ctx.inputs["input_root"]);css=KPSYCSS()
        function prepare(level,bound)
            settings=ctx.contract["numerical_policy"]["refinements"][level]
            progress(Dict("operation"=>"prepare fixed geometry","native_level"=>level,"bound_inverse_GeV"=>bound))
            radial=RadialRule(parton_radial_edges(css,ff,ctx.inputs["input_root"],bound);order=settings["radial_order"])
            context,seconds=sidis_pilot_measure(()->prepare_polarized_sidis(pdf,hel,ff,ctx.nativefile,css,radial,
                k["Q2_GeV2"],k["pT_GeV"],k["x"],k["z"],1;matching_order=settings["matching_order"],
                intervals=settings["intervals"],range_policy=:reference,target=historical_sidis_target(point["target"]),construction=:batched_probes))
            context.first_samples.input_kind==context.second_samples.input_kind==:pinned_inputs || error("physical point requires pinned initial samples")
            resource=Dict{String,Any}("radial_nodes"=>length(radial.b),
                "row_coefficient_payload_bytes"=>context.first_operator.coefficient_payload_bytes+context.second_operator.coefficient_payload_bytes,
                "cold_preparation_seconds_informational"=>seconds,"context_provenance"=>deepcopy(context.provenance))
            context,resource
        end
        function values(context)
            observable=historical_sidis_point_observable(ctx.request,point["id"],context)
            structure=polarized_structure(observable.context,theta)
            amplitude=predict_observable(observable,theta)
            value=Dict{String,Any}("UU"=>structure.UU,"FUT"=>structure.Collins,"amplitude"=>amplitude)
            reference=report["source_factor"]*structure.Collins/structure.UU
            all(isfinite,Base.values(value)) && structure.UU>0 &&
                abs(amplitude-reference)/max(abs(reference),1e-8)<=ctx.contract["numerical_policy"]["observable_algebra_max_scaled"] || error("invalid point/source observable algebra")
            value
        end
        if phase in ("native123","native4")
            if phase=="native4"
                selected=SJ.selection(ctx,run,"native-select3")
                selected["status"]=="NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" && selected["selected_level"]==0 || error("fourth level was not requested globally")
                SJ.audit_point(ctx,run,"native123",index;complete_invocation=true)
            end
            report["levels"]=Dict{String,Any}[]
            for level in (phase=="native123" ? (1:3) : (4:4))
                context,resource=prepare(level,24)
                push!(report["levels"],Dict("level"=>level,"bound_inverse_GeV"=>24,
                    "setting"=>deepcopy(ctx.contract["numerical_policy"]["refinements"][level]),"values"=>values(context),"resources"=>resource))
                context=nothing;GC.gc()
                progress(Dict("operation"=>"native observation retained without diagnostic admission","level"=>level))
            end
        else
            selected=SJ.selection(ctx,run,"native-selection")
            level=selected["selected_level"];level in (3,4) || error("missing common native level")
            report["native_level"]=level
            original=nothing
            if phase=="reference4"
                decision=SJ.selection(ctx,run,"derivative-select3")
                decision["status"]=="DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED" && decision["selected_reference_level"]==0 || error("fourth inward step was not requested globally")
                original=SJ.audit_point(ctx,run,"tailref123",index;complete_invocation=true)
            else
                native_phase=level==3 ? "native123" : "native4"
                native=SJ.audit_point(ctx,run,native_phase,index;complete_invocation=true)
                coarse=only(r for r in native["levels"] if r["level"]==level)
                context12,resource12=prepare(level,12)
                report["tails"]=Dict{String,Any}[
                    Dict("bound_inverse_GeV"=>12,"values"=>values(context12),"resources"=>resource12),
                    Dict("bound_inverse_GeV"=>24,"values"=>deepcopy(coarse["values"]),"resources"=>deepcopy(coarse["resources"]))]
                context12=nothing;GC.gc()
            end
            context,resource=prepare(level,32)
            observable=historical_sidis_point_observable(ctx.request,point["id"],context)
            charged=observable.context;report["values"]=values(context)
            report["kernel_signature"]=bytes2hex(sha256(join((point["id"],context.first_operator.signature,context.second_operator.signature,
                CollinsFit.numerical_snapshot(context.first_samples.bases),CollinsFit.numerical_snapshot(context.second_samples.bases),
                CollinsFit.numerical_snapshot(context.unpolarized_first),CollinsFit.numerical_snapshot(context.unpolarized_second)),";")))
            report["final_preparation_resources"]=resource
            if phase=="tailref123"
                push!(report["tails"],Dict("bound_inverse_GeV"=>32,"values"=>deepcopy(report["values"]),"resources"=>resource))
                a,b=report["tails"][3]["values"],report["tails"][2]["values"]
                report["tail_error"]=maximum(abs(a[key]-b[key])/max(abs(b[key]),1e-8) for key in ("UU","FUT","amplitude"))
                report["tail_error"]<=ctx.contract["numerical_policy"]["physical_value_max_scaled"] || error("unchanged bound24-to32 physical tail target failed")
            end
            progress("preparing separate fixed-row scalar/BigFloat reference")
            independent=ObservableReference.prepare_generic_reference(charged;y=k["y"])
            report["reference_scope"]=Dict("family"=>independent.scope,"native_rows_shared"=>true,
                "initial_functions_shared"=>true,"independent_evolution"=>false,
                "parameter_arithmetic_reference_bits"=>256,"fit_domain_certified"=>false)
            scalar=Dict("UU"=>independent.denominators[1],"FUT"=>independent.polarized_integrals(theta)[1],
                "amplitude"=>report["source_factor"]*independent.predict(theta)[1])
            report["independent_scalar_reference"]=scalar
            report["scalar_checks"]=SJ.scalar_checks(report["values"],scalar,point,ctx)
            all(c["passed"] for c in report["scalar_checks"]) || error("unchanged independent physical scalar target failed")
            value(t)=predict_observable(observable,t)
            candidate=sidis_job_derivatives(value,theta,direction)
            candidate_gradient,candidate_HVP=candidate.gradient,candidate.HVP
            length(candidate_gradient)==length(candidate_HVP)==13 && all(isfinite,candidate_gradient) && all(isfinite,candidate_HVP) || error("nonfinite full13 candidate derivatives")
            report["candidate_derivatives"]=Dict("value"=>value(theta),"gradient"=>candidate_gradient,"HVP"=>candidate_HVP,"direction"=>direction)
            report["reference_steps"]=Dict{String,Any}[]
            setprecision(ctx.contract["derivative_policy"]["reference_bits"]) do
                t=BigFloat.(theta);d=BigFloat.(direction)
                factor=point["publication"]=="1006.4221v1" ? begin
                    y=BigFloat(k["y"]);2*(1-y)/(1+(1-y)^2)
                end : BigFloat(-1)
                reference(x)=factor*independent.predict(x)[1]
                jacobian(x)=ForwardDiff.gradient(reference,x)
                progress("independent256-bit all13 center Jacobian and directional HVP")
                refvalue=reference(t);derivatives=sidis_job_derivatives(reference,t,d)
                refgradient,refHVP=derivatives.gradient,derivatives.HVP
                report["reference_centre"]=Dict("value"=>CollinsFit.sidis_job_decimal(refvalue),
                    "gradient"=>CollinsFit.sidis_job_decimal_vector(refgradient),"HVP"=>CollinsFit.sidis_job_decimal_vector(refHVP),
                    "reference_bits"=>256,"parameters_decimal"=>CollinsFit.sidis_job_decimal_vector(t),
                    "direction_decimal"=>CollinsFit.sidis_job_decimal_vector(d))
                floor=BigFloat(ctx.contract["derivative_policy"]["relative_floor"])
                report["centre_errors_decimal"]=Dict("gradient"=>CollinsFit.sidis_job_decimal(CollinsFit.sidis_job_scaled(candidate_gradient,refgradient,floor)),
                    "HVP"=>CollinsFit.sidis_job_decimal(CollinsFit.sidis_job_scaled(candidate_HVP,refHVP,floor)))
                if original!==nothing
                    original["kernel_signature"]==report["kernel_signature"] &&
                        original["reference_centre"]==report["reference_centre"] &&
                        original["candidate_derivatives"]==report["candidate_derivatives"] && original["values"]==report["values"] &&
                        original["independent_scalar_reference"]==report["independent_scalar_reference"] || error("fourth reference re-preparation changed fixed rows or center operands")
                end
                previous=original===nothing ? nothing : original["reference_steps"][3]
                for refinement in (phase=="tailref123" ? (1:3) : (4:4))
                    progress(Dict("operation"=>"independent all13 inward reference","reference_level"=>refinement))
                    row=sidis_job_reference_step(reference,jacobian,t,d,refinement,ctx.contract["derivative_policy"];
                        base_value=refvalue,base_jacobian=refgradient,reference_HVP=refHVP,
                        candidate_gradient,candidate_HVP,previous)
                    push!(report["reference_steps"],row);previous=row
                    progress(Dict("operation"=>"retained complete decimal reference operands","reference_level"=>refinement))
                end
            end
        end
        SJ.stable(ctx)
        all(SJ.digest(p)==h for (p,h) in report["dependency_sha256"]) || error("point numerical inputs changed during execution")
        report["status"]=SJ.RECORDED
    finally
        pdf===nothing || close(pdf)
    end
    report
end

function main()
    length(ARGS)==7 || error("REQUEST PHASE INDEX NATIVE_BUILD PARTON_BUILD RUN FRESH_REPORT")
    request,phase,index_text,native,inputs,run,output=ARGS;index=parse(Int,index_text)
    phase in SJ.POINT_PHASES || error("unknown point phase")
    abspath(output)==SJ.point_report(abspath(run),phase,index) || error("use canonical fresh point report path")
    (ispath(output)||ispath(output*".progress.json")) && error("preserve earlier point reports")
    BLAS.set_num_threads(1)
    report=Dict{String,Any}("status"=>"FAIL","scope"=>CollinsFit.HISTORICAL_SIDIS_JOB_SCOPE,
        "phase"=>phase,"point_index"=>index,SJ.false_flags()...)
    try
        ctx=SJ.bindings(request,native,inputs);merge!(report,SJ.header(ctx,phase;index))
        run_sidis_point!(report,ctx,phase,index,abspath(run),output*".progress.json")
    catch exception
        report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));writejson(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"phase"=>phase,"point_index"=>index,"error"=>get(report,"error",nothing))));flush(stdout)
    report["status"]==SJ.RECORDED || exit(2)
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
