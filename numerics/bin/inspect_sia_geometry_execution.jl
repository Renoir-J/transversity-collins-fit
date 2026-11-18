# Read-only source/runtime/prerequisite inspection; no physical job or fitting admission.
include("SIACompleteJobs.jl")
using JSON
function main(args=ARGS)
    length(args)==2 || error("usage: inspect_sia_geometry_execution.jl MANIFEST EXTERNAL_SHA")
    c=SIACompleteJobs.load_manifest(args...)
    SIACompleteJobs.hashes(SIACompleteJobs.base_inputs(c))
    SIACompleteJobs.source_controls(c)
    SIACompleteJobs.recheck(c)
    println(JSON.json(Dict("status"=>"SOURCE_AND_SYNTHETIC_PREREQUISITES_REVIEWED_NOT_PHYSICAL_DISPATCH",
        "adapter_head"=>c.m["adapter_head"],"source_plan_sha256"=>c.source.plan_sha,"manifest_sha256"=>c.sha,
        "physical_jobs_started"=>false,"fit_performed"=>false,"accepted_fitting_gate"=>false)))
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
