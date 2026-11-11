using CollinsFit, JSON
length(ARGS) == 1 || error("usage: julia --project=numerics numerics/bin/run_job.jl job.json")
CollinsFit.run_job(JSON.parsefile(ARGS[1]))

