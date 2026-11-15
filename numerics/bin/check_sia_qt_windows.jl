using CollinsFit, JSON
include("SIAQTChecks.jl")

length(ARGS)==1 || error("usage: check_sia_qt_windows.jl output.json")
output=abspath(only(ARGS))
ispath(output) && error("refusing to overwrite an earlier qT-window report")
try
    report=validate_sia_qt_windows(@__FILE__)
    writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],
        "input_scope"=>report["input_scope"],"checks"=>length(report["checks"]),
        "failed"=>filter(c->!c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory",
        "input_scope"=>"explicitly_synthetic_fixed_rows","production_ready"=>false,
        "physical_extraction"=>false,"G2PERF_pass"=>false,"sources"=>sia_qt_sources(@__FILE__),
        "error"=>sprint(showerror,exception)))
    rethrow()
end
