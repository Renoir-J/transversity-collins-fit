# Independent source-program comparisons for reconstructed input components.
using CollinsFit, JSON, SHA, Dates, ForwardDiff

const ROOT = normpath(joinpath(@__DIR__,"../.."))
const CONTRACT = joinpath(ROOT,"config/parton-contract.json")
const SOURCES = vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml"),
    joinpath(ROOT,"scripts/prepare_parton_inputs.py"),joinpath(ROOT,"src/collinsflow/parton.py")],
    [joinpath(ROOT,p) for p in ("bridges/lhapdf.cpp","bridges/check_lhapdf.cpp",
                               "fortran/check_dssv08.f90","fortran/check_dss14.f90")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes() = Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function reference_rows(executable,cwd,queries,nf;args=String[])
    payload = string(length(queries))*"\n"*join([join(string.(q)," ") for q in queries],"\n")*"\n"
    command = Cmd(Cmd(vcat([executable],args));dir=cwd)
    output = read(pipeline(command;stdin=IOBuffer(payload)),String)
    lines = split(output,'\n';keepempty=false)
    rows = [parse.(Float64,split(line)[2:end]) for line in lines if startswith(line,"RESULT ")]
    length(rows) == length(queries) && all(r -> length(r)==nf && all(isfinite,r),rows) ||
        error("incomplete/nonfinite source-program reference output")
    hcat(rows...),filter(line -> !startswith(line,"RESULT "),lines)
end

function parton_checks(buildfile)
    sources = source_hashes()
    build = CollinsFit.verify_parton_build(buildfile)
    build_hash = bytes2hex(sha256(read(buildfile)))
    contract = JSON.parsefile(CONTRACT)
    policy = contract["validation_policy"]
    floor = policy["relative_floor"]
    scaled(a,b) = maximum(abs.(a.-b)./max.(abs.(b),floor))
    checks = Dict{String,Any}[]
    check(name,value,tolerance) = push!(checks,Dict("name"=>name,"error"=>value,
        "tolerance"=>tolerance,"passed"=>isfinite(value) && value<=tolerance))
    rejected(f) = try f(); false catch; true end
    root = build["input_root"]
    dssv,dss = load_dssv08(root),load_dss14(root)
    # Include exact knots, the true source boundaries, and off-grid coordinates.
    dssv_queries = [(x,q) for x in (1e-5,5e-5,.0013,.013,.08,.137,.4,.651,.9,.999,1.)
        for q in (.8,1.21,2.4,5.,31.,100.,1234.,7e5,1e6)]
    ref = build["references"]["DSSV08"]
    dssv_reference,dssv_log = reference_rows(ref["executable"],ref["cwd"],dssv_queries,6)
    native = hcat([weighted_components(dssv,x,q) for (x,q) in dssv_queries]...)
    check("DSSV08 native interpolation versus unchanged Fortran",scaled(native,dssv_reference),policy["legacy_max_scaled"])
    conversion = Float64[]
    for (i,(x,q)) in enumerate(dssv_queries)
        uv,dv,ub,db,_,g = dssv_reference[:,i]
        expected = [uv+ub,dv+db,ub,db,g]./x
        actual = [parton_helicity(dssv,p,x,q) for p in (2,1,-2,-1,21)]
        push!(conversion,scaled(actual,expected))
    end
    check("DSSV08 valence-to-quark and x removal",maximum(conversion),policy["legacy_max_scaled"])
    check("DSSV08 source range and unsupported flavor rejection",all((
        rejected(() -> weighted_components(dssv,1e-6,2.4)),
        rejected(() -> weighted_components(dssv,.1,.79;policy=:reference)),
        rejected(() -> weighted_components(dssv,.1,1.01e6)),
        rejected(() -> parton_helicity(dssv,3,.1,2.4)))) ? 0. : 1.,0.)
    dss_queries = [(c,z,q) for c in (-1,0,1) for z in (.05,.075,.15,.237,.5,.713,.95,.999,1.)
        for q in (1.,1.37,2.4,7.,33.,100.,6e4,1e5)]
    ref = build["references"]["DSS14"]
    dss_reference,dss_log = reference_rows(ref["executable"],ref["cwd"],dss_queries,9)
    native = hcat([CollinsFit.pion_components(weighted_components(dss,z,q),c) for (c,z,q) in dss_queries]...)
    check("DSS14 all pion charges versus unchanged Fortran",scaled(native,dss_reference),policy["legacy_max_scaled"])
    flavors = (2,-2,1,-1,3,-3,4,5,21)
    conversion = [scaled([pion_ff(dss,p,c,z,q) for p in flavors],dss_reference[:,i]./z)
        for (i,(c,z,q)) in enumerate(dss_queries)]
    check("DSS14 z removal and PDG flavor mapping",maximum(conversion),policy["legacy_max_scaled"])
    heavy = [abs(pion_ff(dss,4,1,.2,100.)-pion_ff(dss,-4,1,.2,100.)),
             abs(pion_ff(dss,5,-1,.2,100.)-pion_ff(dss,-5,-1,.2,100.))]
    check("DSS14 per-heavy-flavor equality",maximum(heavy),0.)
    check("DSS14 strict range rejection",all((
        rejected(() -> pion_ff(dss,2,1,.049,2.4)),
        rejected(() -> pion_ff(dss,2,1,.2,.9)),
        rejected(() -> pion_ff(dss,2,1,.2,1.01e5)),
        rejected(() -> pion_ff(dss,2,2,.2,2.4)))) ? 0. : 1.,0.)
    low_queries = [(c,.23,.6) for c in (-1,0,1)]
    low_reference,low_log = reference_rows(ref["executable"],ref["cwd"],low_queries,9)
    low = hcat([CollinsFit.pion_components(weighted_components(dss,z,q;policy=:reference),c)
        for (c,z,q) in low_queries]...)
    check("DSS14 explicit reference scale extrapolation",scaled(low,low_reference),policy["legacy_max_scaled"])
    domain = parton_domain(dss,.23,.6;policy=:reference)
    check("DSS14 extrapolation remains labeled",!domain["within_grid"] &&
        !domain["production_ready"] && domain["policy"]=="reference" && !isempty(low_log) ? 0. : 1.,0.)

    # Parse the raw official table independently of both library APIs.
    grid = readlines(CollinsFit.verified_parton_asset(root,"CT10_central"))
    findall(==("---"),strip.(grid)) == [3,length(grid)] || error("unexpected CT10 subgrid structure")
    xs = parse.(Float64,split(grid[4]))
    qs = parse.(Float64,split(grid[5])).^2
    pids = parse.(Int,split(grid[6]))
    length(grid)-7 == length(xs)*length(qs) || error("unexpected CT10 raw grid dimensions")
    raw = hcat([parse.(Float64,split(line)) for line in grid[7:end-1]]...)
    check("CT10 flavors match the source contract",pids==contract["CT10"]["flavors"] ? 0. : 1.,0.)
    ct10 = CT10Input(buildfile)
    reference_input = nothing
    try
        # Every node/flavor, including heavy thresholds and zero/negative entries.
        xx = repeat(xs;inner=length(qs))
        qq = repeat(qs;outer=length(xs))
        values = ct10_batch(ct10,pids,xx,qq)
        check("CT10 raw grid nodes with explicit ForcePositive=1",scaled(values,max.(0.,raw)),policy["ct10_grid_max_scaled"])
        queries = [(p,x,q) for p in pids for x in (1e-8,.00013,.013,.137,.61,.999,1.)
            for q in (1.3^2,2.4,4.75^2,31.,1e5,1e10)]
        reference,ct10_log = reference_rows(build["ct10_reference"],root,queries,1;args=[root,"strict"])
        batch = ct10_batch(ct10,pids,[x for x in (1e-8,.00013,.013,.137,.61,.999,1.) for q in (1.3^2,2.4,4.75^2,31.,1e5,1e10)],
            [q for x in (1e-8,.00013,.013,.137,.61,.999,1.) for q in (1.3^2,2.4,4.75^2,31.,1e5,1e10)])
        # Factory executable is pid-major; the ABI is query-major.
        actual = [batch[f,j] for f in eachindex(pids) for j in axes(batch,2)]
        check("CT10 C ABI ordering versus separate factory executable",scaled(actual,vec(reference)),policy["ct10_abi_max_scaled"])
        ratios = [abs(parton_pdf(ct10,p,x,q)*x-batch[f,j])/max(abs(batch[f,j]),floor) for (f,p) in enumerate(pids)
            for (j,(x,q)) in enumerate([(x,q) for x in (1e-8,.00013,.013,.137,.61,.999,1.) for q in (1.3^2,2.4,4.75^2,31.,1e5,1e10)])]
        check("CT10 public unweighted PDF removes x",maximum(ratios),policy["conversion_max_scaled"])
        check("CT10 strict range and flavor rejection",all((
            rejected(() -> parton_pdf(ct10,2,.1,.6)),
            rejected(() -> parton_pdf(ct10,2,1e-9,2.4)),
            rejected(() -> parton_pdf(ct10,6,.1,2.4)),
            rejected(() -> ct10_batch(ct10,[2],[.1,.2],[2.4])))) ? 0. : 1.,0.)
        check("CT10 kinematic AD cannot silently lose derivatives",
            rejected(() -> ForwardDiff.derivative(x -> parton_pdf(ct10,2,x,2.4),.1)) ? 0. : 1.,0.)
        reference_input = CT10Input(buildfile;policy=:reference)
        low_queries = [(p,.23,.6) for p in pids]
        low_reference,low_ct10_log = reference_rows(build["ct10_reference"],root,low_queries,1;args=[root,"reference"])
        low = vec(ct10_batch(reference_input,pids,[.23],[.6]))
        check("CT10 explicit continuation matches separate executable",scaled(low,vec(low_reference)),policy["ct10_abi_max_scaled"])
        domain = parton_domain(reference_input,.23,.6)
        check("CT10 extrapolation remains labeled",!domain["within_grid"] &&
            !domain["production_ready"] && domain["policy"]=="reference" ? 0. : 1.,0.)
        close(ct10)
        close(ct10)
        check("closed CT10 context rejects reuse",rejected(() -> parton_pdf(ct10,2,.1,2.4)) ? 0. : 1.,0.)
        gridfile = CollinsFit.verified_parton_asset(root,"DSSV08_grid")
        original = read(gridfile)
        tamper_rejected = try
            write(gridfile,vcat(original,UInt8[0x0a]))
            rejected(() -> load_dssv08(root)) && rejected(() -> CollinsFit.verify_parton_build(buildfile))
        finally
            write(gridfile,original)
        end
        check("changed source bytes invalidate loader and build evidence",tamper_rejected ? 0. : 1.,0.)
        check("source and build evidence stable during job",source_hashes()==sources &&
            bytes2hex(sha256(read(buildfile)))==build_hash &&
            CollinsFit.verify_parton_build(buildfile)==build ? 0. : 1.,0.)
        Dict("scope"=>"reconstructed_theory","production_ready"=>false,"original_KPSY_reproduction"=>false,
            "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
            "contract_id"=>contract["contract_id"],"build_receipt"=>abspath(buildfile),"build_sha256"=>build_hash,
            "build"=>build,"completed_at"=>string(now(UTC)),"julia_version"=>string(VERSION),
            "comparison_counts"=>Dict("DSSV08_queries"=>length(dssv_queries),"DSS14_queries"=>length(dss_queries),
                "CT10_raw_node_flavor_values"=>length(raw),"CT10_factory_queries"=>length(queries)),
            "reference_messages"=>Dict("DSSV08"=>dssv_log,"DSS14"=>dss_log,"DSS14_low_scale"=>low_log,
                "CT10"=>ct10_log,"CT10_low_scale"=>low_ct10_log),
            "limitations"=>contract["limitations"],
            "reference_independence"=>"DSS ports compare against unchanged Fortran; CT10 uses the official library in both ABI/factory paths plus independent raw-node values, not an independent interpolation algorithm")
    finally
        close(ct10)
        !isnothing(reference_input) && close(reference_input)
    end
end

buildfile = isempty(ARGS) ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[1]
output = length(ARGS)<2 ? joinpath(ROOT,"runs/parton-inputs/validation.json") : ARGS[2]
try
    report = parton_checks(buildfile)
    writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "counts"=>report["comparison_counts"],"failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
