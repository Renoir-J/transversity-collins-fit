using CollinsFit, JSON, SHA, Dates
include("SIAJobBindings.jl")
include("SIASourceControls.jl")
module FrozenSIAControlsProducer
end

function load_legacy(frozen)
    Base.include(FrozenSIAControlsProducer,joinpath(frozen,"numerics/bin/predict_historical_sia_integrals.jl"))
    FrozenSIAControlsProducer
end

function main(args=ARGS)
    length(args)==4 || error("usage: run_sia_source_controls.jl FROZEN SOURCE_PLAN EXTERNAL_SHA FRESH_REPORT")
    frozen,plan,sha,output=args
    context=SIAJobBindings.load(frozen,plan,sha)
    output=SIAJobBindings.fresh_output(output,plan)
    report=merge(SIASourceControls.flags(),Dict{String,Any}("status"=>"FAIL","scope"=>SIASourceControls.SCOPE,
        "theory_scope"=>"synthetic_closure","source_controls_success"=>false,"started_at"=>string(now(UTC))))
    SIAJobBindings.atomic_report(output,report)
    try
        legacy=load_legacy(context.frozen)
        report=SIASourceControls.produce(context.request.data,context.parent,context.binding,
            legacy.HSP.source_projection,legacy.HIR.mathematical_controls)
        SIAJobBindings.recheck(context)
    catch exception
        report["status"]="FAIL";report["source_controls_success"]=false
        report["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        report["completed_at"]=string(now(UTC));SIAJobBindings.atomic_report(output,report)
    end
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"fit_performed"=>false)))
    report["status"]==SIASourceControls.COMPUTED ? 0 : 2
end
if abspath(PROGRAM_FILE)==(@__FILE__);exit(main());end
