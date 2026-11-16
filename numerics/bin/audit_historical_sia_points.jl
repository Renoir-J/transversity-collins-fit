include("HistoricalSIAPointAudit.jl")
using CollinsFit, JSON, SHA, Dates
function main(args=ARGS)
    length(args)==3 || error("usage: audit_historical_sia_points.jl REQUEST REPORT FRESH_OUTPUT")
    request,report,output=args
    ispath(output) && error("use a fresh audit output")
    result=Dict{String,Any}("status"=>"FAIL","review_status"=>"FAIL",
        "scope"=>"historical_reconstructed_sia_point_prediction_diagnostic","accepted_prediction_diagnostic"=>false,
        "accepted_fitting_gate"=>false,"production_ready"=>false,"physical_extraction"=>false,
        "historical_fit_ready"=>false,"experimental_data_admitted"=>false,"fit_performed"=>false,
        "G2PERF_pass"=>false,"interval_coverage_verified"=>false)
    try
        result=HistoricalSIAPointAudit.audit(request,report)
    catch exception
        result["error"]=sprint(showerror,exception,catch_backtrace())
    finally
        result["completed_at"]=string(now(UTC));writejson(output,result)
    end
    println(JSON.json(result))
    result
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    main()["status"]=="PASS" || exit(2)
end
