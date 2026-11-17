using CollinsFit, JSON, SHA, Dates
include("SIAJobBindings.jl")
include("SIASourceControls.jl")
module FrozenSIAControlsProducer
end
const LEGACY_ROOT=Ref{Union{Nothing,String}}(nothing)

function load_legacy(frozen)
    if LEGACY_ROOT[]===nothing
        Base.include(FrozenSIAControlsProducer,joinpath(frozen,"numerics/bin/predict_historical_sia_integrals.jl"))
        LEGACY_ROOT[]=frozen
    else
        LEGACY_ROOT[]==frozen || error("do not mix frozen numerical libraries in one worker")
    end
    FrozenSIAControlsProducer
end

function produce_loaded(legacy,data,parent,binding)
    SIASourceControls.produce(data,parent,binding,legacy.HSP.source_projection,legacy.HIR.mathematical_controls)
end

function produce_controls(context)
    legacy=load_legacy(context.frozen)
    # Enter the new method world once, outside the mathematical loop. The
    # included files were verified before loading; no data selects executable code.
    Base.invokelatest(produce_loaded,legacy,context.request.data,context.parent,context.binding)
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
        report=produce_controls(context)
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
