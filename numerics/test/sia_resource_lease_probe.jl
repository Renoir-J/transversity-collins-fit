# Explicit private-lock integration probe. It cannot dispatch numerical jobs.
include("../bin/SIACompleteJobs.jl")
using JSON
const C=SIACompleteJobs

function main(args=ARGS)
    length(args)==2 || error("usage: PRIVATE_FIXTURE FRESH_REPORT")
    input,output=args
    ispath(output) && error("retain occupied probe output")
    fixture=JSON.parsefile(input)
    fixture["scope"]=="explicitly_synthetic_private_resource_cases_not_dispatch" || error("private fixture scope required")
    results=Any[]
    for entry in fixture["cases"]
        path=entry["expected_path"]
        startswith(path,"/tmp/collinsflow-sia-") && error("probe cannot use a production resource lock")
        result=Dict{String,Any}("name"=>entry["name"],"accepted"=>false)
        try
            proof=C.resource_lease_identity(entry["proof"],path)
            result["accepted"]=true
            result["owner_pid"]=proof.pid
        catch error
            result["error"]=sprint(showerror,error)
        end
        push!(results,result)
    end
    report=Dict("scope"=>fixture["scope"],"cases"=>results,"production_locks_acquired"=>false,
        "physical_Julia_jobs_started"=>false,"fit_performed"=>false,"production_ready"=>false,"accepted_fitting_gate"=>false)
    write(output,JSON.json(report)*"\n")
    println(JSON.json(Dict("scope"=>report["scope"],"cases"=>length(results),"physical_Julia_jobs_started"=>false)))
    0
end
if abspath(PROGRAM_FILE)==(@__FILE__);exit(main());end
