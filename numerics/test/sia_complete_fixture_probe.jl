# Fresh-process synthetic saved-ledger reader; never a physical worker.
include("SIACompleteJobFixture.jl")
using JSON
function main(args=ARGS)
    length(args)==2 || error("usage: SYNTHETIC_LEDGER FRESH_OUTPUT")
    directory,output=args;ispath(output) && error("retain occupied probe output")
    result=Dict{String,Any}("status"=>"FAIL","fixture_scope"=>SIACompleteJobFixture.FIXTURE,"theory_scope"=>"synthetic_closure",
        "physical_SIA_jobs_started"=>false,"fit_performed"=>false,"production_ready"=>false,"accepted_fitting_gate"=>false)
    try
        merge!(result,SIACompleteJobFixture.inspect(directory));result["status"]="SYNTHETIC_FULL72_SAVED_LEDGER_REVIEWED_NOT_PHYSICAL_ADMISSION"
    catch error
        result["error"]=sprint(showerror,error,catch_backtrace())
    end
    write(output,JSON.json(result)*"\n")
    println(JSON.json(Dict("status"=>result["status"],"fixture_scope"=>result["fixture_scope"],"fit_performed"=>false)))
    result["status"]=="FAIL" ? 2 : 0
end
if abspath(PROGRAM_FILE)==(@__FILE__);exit(main());end
