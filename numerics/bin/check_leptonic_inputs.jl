# Pinned-input unpolarized reconstruction pilot; never a fitted-data or production receipt.
using CollinsFit, JSON, SHA, Dates

const ROOT = normpath(joinpath(@__DIR__,"../.."))
const CONTRACT = joinpath(ROOT,"config/leptonic-contract.json")
const SOURCES = vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/parton-contract.json"),
    joinpath(ROOT,"config/parton-sources.json"),joinpath(ROOT,"config/css-contract.json"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p -> endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes() = Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)

function input_edges(model,dss,root,B)
    grid = readlines(CollinsFit.verified_parton_asset(root,"CT10_central"))
    q2knots = vcat(parse.(Float64,split(grid[5])).^2,dss.q2)
    c02 = 4exp(-2Float64(Base.MathConstants.eulergamma))
    bknots = [inv(sqrt(q/c02-inv(model.bmax^2))) for q in q2knots if q>c02/model.bmax^2]
    base = [0.,1e-7,1e-6,1e-5,1e-4,.001,.005,.02,.05,.1,.2,.4,.7,1.,1.5,2.,3.,4.,6.,8.,12.,16.,20.,24.]
    sort(unique(vcat(filter(<(B),base),filter(b -> 0<b<B,bknots),B)))
end

function input_checks(buildfile)
    sources = source_hashes()
    build = CollinsFit.verify_parton_build(buildfile)
    build_hash = bytes2hex(sha256(read(buildfile)))
    policy = JSON.parsefile(CONTRACT)["validation_policy"]
    checks = Dict{String,Any}[]
    check(name,value,tolerance) = push!(checks,Dict("name"=>name,"error"=>value,
        "tolerance"=>tolerance,"passed"=>isfinite(value) && value<=tolerance))
    scaled(a,b) = maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
    function rejected(f,expected)
        try
            f()
            false
        catch exception
            occursin(expected,sprint(showerror,exception))
        end
    end
    model,dss = KPSYCSS(),load_dss14(build["input_root"])
    strict = CT10Input(buildfile)
    pdf = nothing
    try
        small = RadialRule([0.,.1,1.,4.];order=4)
        check("strict SIDIS preparation rejects unsupported scales",
            rejected(() -> prepare_unpolarized_sidis(strict,dss,model,small,2.4,.3,.1,.3,1),"outside") ? 0. : 1.,0.)
        check("strict SIA preparation rejects unsupported scales",
            rejected(() -> prepare_unpolarized_sia(dss,model,small,110.,.3,.3,.5,1,-1),"outside strict range") ? 0. : 1.,0.)
        pdf = CT10Input(buildfile;policy=:reference)
        check("mixed implicit range policies are rejected",
            rejected(() -> prepare_unpolarized_sidis(pdf,dss,model,small,2.4,.3,.1,.3,1),"policies must agree") ? 0. : 1.,0.)
        ps = [0.,.1,.3,.5,.7]
        refinements = []
        for setting in policy["input_refinements"]
            radial = RadialRule(input_edges(model,dss,build["input_root"],24.);order=setting["radial_order"])
            sidis = prepare_unpolarized_sidis(pdf,dss,model,radial,2.4,0.,.1,.3,1;
                matching_order=setting["matching_order"],range_policy=:reference)
            sia = prepare_unpolarized_sia(dss,model,radial,110.,0.,.3,.5,1,-1;
                matching_order=setting["matching_order"],range_policy=:reference)
            values = []
            for (process,Q2,z,z2,prepared) in ((:SIDIS,2.4,.3,nothing,sidis),(:SIA,110.,.3,.5,sia))
                predictions = [begin
                    plan = leptonic_plan(process,model,radial,Q2,P,z;z2=z2)
                    leptonic_unpolarized(plan,prepared.first,prepared.second)
                end for P in ps]
                push!(values,Dict("process"=>string(process),"P_GeV"=>ps,"unpolarized"=>predictions,
                    "provenance"=>prepared.provenance))
            end
            push!(refinements,Dict("matching_order"=>setting["matching_order"],
                "radial_order"=>setting["radial_order"],"radial_nodes"=>length(radial.b),"predictions"=>values))
        end
        final,previous = last(refinements)["predictions"],refinements[end-1]["predictions"]
        for (f,p) in zip(final,previous)
            process = f["process"]
            check("$process pinned-input quadrature refinement",scaled(f["unpolarized"],p["unpolarized"]),policy["input_prediction_max_scaled"])
            check("$process pinned-input finite positive zero-P denominator",
                all(isfinite,f["unpolarized"]) && first(f["unpolarized"])>0 ? 0. : 1.,0.)
            meta = f["provenance"]
            audit = meta["source_range_audits"]
            labeled = meta["scope"]=="reconstructed_theory" && !meta["production_ready"] &&
                meta["range_policy"]=="reference" && audit["DSS14_outside_strict_nodes"]>0
            process=="SIDIS" && (labeled &= audit["CT10_outside_strict_nodes"]>0)
            check("$process unsupported scales remain explicitly labeled",labeled ? 0. : 1.,0.)
        end
        check("source and build hashes stable during input pilot",source_hashes()==sources &&
            bytes2hex(sha256(read(buildfile)))==build_hash && CollinsFit.verify_parton_build(buildfile)==build ? 0. : 1.,0.)
        Dict("scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
            "status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","checks"=>checks,"sources"=>sources,
            "contract_id"=>"KPSY-LEPTONIC-001","build_sha256"=>build_hash,"build"=>build,
            "refinements"=>refinements,"completed_at"=>string(now(UTC)),"julia_version"=>string(VERSION),
            "input_scope"=>"pinned CT10/DSS14; explicit reference extrapolation, no experimental fit",
            "limitations"=>["Unpolarized components only; no initial transversity/Collins fit model",
                "Reference extrapolation is not a production prescription", "No independent physical-data benchmark or likelihood admission"])
    finally
        close(strict)
        !isnothing(pdf) && close(pdf)
    end
end

buildfile = isempty(ARGS) ? joinpath(ROOT,".runtime/parton-inputs/build.json") : ARGS[1]
output = length(ARGS)<2 ? joinpath(ROOT,"runs/leptonic-inputs/validation.json") : ARGS[2]
try
    report = input_checks(buildfile)
    writejson(output,report)
    println(JSON.json(Dict("scope"=>report["scope"],"status"=>report["status"],"checks"=>length(report["checks"]),
        "failed"=>filter(c -> !c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("scope"=>"reconstructed_theory","production_ready"=>false,"status"=>"FAIL",
        "error"=>sprint(showerror,exception)))
    rethrow()
end
