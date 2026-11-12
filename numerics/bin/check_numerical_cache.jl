# Real locked-source cache checks, separate from scientific or fit admission.
using CollinsFit, JSON, SHA, LinearAlgebra

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/numerical-cache-contract.json")
const SOURCES=vcat([@__FILE__,CONTRACT,joinpath(ROOT,"config/polarized-process-contract.json"),
    joinpath(ROOT,"config/leptonic-observable-contract.json"),joinpath(ROOT,"config/parton-sources.json"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))

function cache_checks(nativefile,inputfile)
    sources=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)
    native=CollinsFit.verified_hoppet_build(nativefile);inputs=CollinsFit.verify_parton_build(inputfile)
    builds=Dict(p=>bytes2hex(sha256(read(p))) for p in (nativefile,inputfile))
    checks=Dict{String,Any}[];mutations=Dict{String,Any}[];processes=Dict{String,Any}()
    check(name,passed)=push!(checks,Dict("name"=>name,"error"=>(passed ? 0. : 1.),"tolerance"=>0.,"passed"=>passed))
    function mutation(name,array,index,prediction,expected)
        baseline=prediction();old=array[index];changed=false;rejected=false;stage=nothing;reason=nothing
        try
            try
                array[index]=old+(old isa Integer ? 1 : 1e-8)
                changed=true
            catch exception
                reason=sprint(showerror,exception)
                rejected=occursin(expected,reason);stage="write"
            end
            if changed
                try
                    prediction()
                catch exception
                    reason=sprint(showerror,exception)
                    rejected=occursin(expected,reason);stage="application"
                end
            end
        finally
            changed && (array[index]=old)
        end
        restored=prediction()==baseline
        check(name*" rejected",rejected);check(name*" restored",restored)
        push!(mutations,Dict("name"=>name,"rejected"=>rejected,"rejected_at"=>stage,"reason"=>reason,
            "restored_exactly"=>restored))
    end
    model=KPSYCSS();ff=load_dss14(inputs["input_root"]);hel=load_dssv08(inputs["input_root"])
    pdf=CT10Input(inputfile;policy=:reference)
    radial=RadialRule(parton_radial_edges(model,ff,inputs["input_root"],24.);order=8)
    theta=collect(kpsy_published_parameters().values)
    try
        for process in (:SIDIS,:SIA)
            label=string(process)
            context=process==:SIDIS ?
                prepare_polarized_sidis(pdf,hel,ff,nativefile,model,radial,2.4,.3,.1,.3,1;
                    matching_order=4,intervals=512,range_policy=:reference) :
                prepare_polarized_sia(ff,nativefile,model,radial,110.,.3,.3,.5,1,-1;
                    matching_order=4,intervals=512,range_policy=:reference)
            momentum=MomentumRule(.05,.35;order=8)
            observable=process==:SIDIS ?
                sidis_observable(context;definition=:structure_function_ratio,momentum=momentum,source="software cache-integrity diagnostic") :
                sia_observable(context;ratio=:UL,momentum=momentum,
                    angular=SIAMeanOfRatio(.7;source="diagnostic fixed factor, not experiment",averaging_measure="fixed diagnostic measure"),
                    source="software cache-integrity diagnostic")
            prediction()=predict_observable(observable,theta)
            value=prediction();operator=context.first_operator;block=operator.blocks[1]
            check(label*" finite diagnostic observable",isfinite(value))
            check(label*" actual radial evolution scales",operator.Q2==css_mu_b2.(Ref(model),radial.b))
            check(label*" explicit reference source policy",context.provenance["range_policy"]=="reference" &&
                context.provenance["initial_source_kind"]=="pinned_inputs")
            mutation(label*" row coefficient",block.rows,min(1001,length(block.rows)),prediction,"read-only")
            mutation(label*" row write-through view",view(block.rows,1:2,1),1,prediction,"read-only")
            for (name,array) in (("block index",block.indices),("operator grid",operator.grid),
                                 ("operator query",operator.queries),("operator scales",operator.Q2),
                                 ("operator times",operator.tau),("operator physical alpha",operator.alpha))
                mutation(label*" "*name,array,max(1,length(array)÷2),prediction,"mutated HOPPET point settings/grid")
            end
            samples=context.first_samples
            for (name,array) in (("initial grid",samples.grid),("initial logarithm",samples.logarithm),
                                 ("initial complement logarithm",samples.complement_logarithm),("initial bases",samples.bases))
                mutation(label*" "*name,array,max(1,length(array)÷2),prediction,"mutated polarized initial samples")
            end
            plan=context.leptonic
            for (name,array) in (("first provider matrix",context.unpolarized_first),("second provider matrix",context.unpolarized_second),
                                 ("radial nodes",plan.radial.b),("radial weights",plan.radial.db),("radial edges",plan.radial.edges),
                                 ("CSS nodes",plan.model.nodes),("CSS weights",plan.model.weights),
                                 ("flavors",plan.flavors),("antiparticle map",plan.antiparticle),("charges squared",plan.charge2),
                                 ("UU weights",plan.uu_weights),("Collins weights",plan.collins_weights),
                                 ("UU damping",plan.uu_damping),("fixed exponent",plan.fixed_exponent),
                                 ("fragmentation shape",plan.fragmentation_shape))
                mutation(label*" "*name,array,max(1,length(array)÷2),prediction,"mutated polarized process context")
            end
            for (name,array) in (("momentum nodes",momentum.nodes),("momentum weights",momentum.weights))
                mutation(label*" "*name,array,max(1,length(array)÷2),prediction,"mutated momentum quadrature")
            end
            for (name,array) in (("observable UU weights",observable.uu_weights),("observable Collins weights",observable.collins_weights),
                                 ("observable denominators",observable.denominators))
                mutation(label*" "*name,array,max(1,length(array)÷2),prediction,"mutated observable settings, weights or metadata")
            end
            check(label*" unchanged after all mutations",prediction()==value)
            processes[label]=Dict("diagnostic_observable"=>value,"radial_nodes"=>length(radial.b),
                "intervals"=>512,"matching_order"=>4,"momentum_order"=>8,
                "source_policy"=>context.provenance["range_policy"],
                "row_payload_bytes"=>operator.coefficient_payload_bytes+
                    (context.second_operator===nothing ? 0 : context.second_operator.coefficient_payload_bytes))
        end
        check("source/build stability",all(bytes2hex(sha256(read(joinpath(ROOT,p))))==h for (p,h) in sources) &&
            all(bytes2hex(sha256(read(p)))==h for (p,h) in builds) &&
            CollinsFit.verified_hoppet_build(nativefile)==native && CollinsFit.verify_parton_build(inputfile)==inputs)
        contract=JSON.parsefile(CONTRACT)
        Dict("status"=>(all(c["passed"] for c in checks) ? "PASS" : "FAIL"),
            "scope"=>contract["scope"],"contract_id"=>contract["id"],"production_ready"=>false,"physical_extraction"=>false,
            "G2PERF_pass"=>false,"checks"=>checks,"sources"=>sources,"build_sha256"=>builds,
            "native_build"=>native,"parton_build"=>inputs,"processes"=>processes,"mutations"=>mutations,
            "Julia_version"=>string(VERSION),
            "limitations"=>contract["limitations"],
            "input_scope"=>"locked physical parton/evolution inputs; explicit reference mu_b scale policy; diagnostic angular factor and finite momentum rule, not an experimental observable or fit")
    finally
        close(pdf)
    end
end

nativefile=length(ARGS)>=1 ? ARGS[1] : joinpath(ROOT,".runtime/hoppet/build.json")
inputfile=length(ARGS)>=2 ? ARGS[2] : joinpath(ROOT,".runtime/parton-inputs/build.json")
output=length(ARGS)>=3 ? ARGS[3] : joinpath(ROOT,"runs/numerical-cache/integrity.json")
try
    report=cache_checks(nativefile,inputfile);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),
        "mutation_cases"=>length(report["mutations"]),"failed"=>filter(c->!c["passed"],report["checks"]))))
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"fixed_numerical_cache_integrity",
        "production_ready"=>false,"physical_extraction"=>false,"error"=>sprint(showerror,exception)))
    rethrow()
end
