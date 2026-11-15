module CouplingMomentReference

using CollinsFit, JSON, SHA, LinearAlgebra
module Producer
include("check_native_coupling.jl")
end
using .Producer.NativeCouplingReference
include("CollinearReference.jl")
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/coupling-moment-contract.json")
export verify_coupling_report, inspect_coupling_report, moment_profiles, study_sources,
       initial_value_gradient, verify_analytic_report

require(ok,label)=ok ? nothing : error(label)
flatten(x)=Float64.(reduce(vcat,x;init=Float64[]))
near(a,b)=length(a)==length(b) && all(isfinite,a) && all(isfinite,b) && all(isapprox.(a,b;rtol=2e-13,atol=2e-15))

function inspect_coupling_report(d,reportfile,buildfile;deep=true)
    C=JSON.parsefile(Producer.CONTRACT);P=C["validation_policy"]
    RUN=dirname(abspath(reportfile));REF_CACHE=Dict{Any,Any}()
    build=CollinsFit.verified_hoppet_build(buildfile)
    require(abspath(d["build_file"])==abspath(buildfile) && d["native_build"]==build &&
        d["build_sha256"]==bytes2hex(sha256(read(buildfile))),"current build binding")
    for path in ("src/collinsflow/hoppet.py","src/collinsflow/provenance.py","fortran/hoppet_bridge.f90",
        "fortran/native_coupling_reference.f90","scripts/prepare_hoppet.py","scripts/prepare_native_coupling.py",
        "config/native-coupling-contract.json","config/hoppet-contract.json","config/evolution-contract.json",
        "config/matched-rows-contract.json","config/hoppet-batch-contract.json")
        matches=filter(p->endswith(replace(p,'\\'=>'/'),"/"*path),collect(keys(build["source_hashes"])))
        require(length(matches)==1 && build["source_hashes"][only(matches)]==
            bytes2hex(sha256(read(joinpath(ROOT,path)))),"current compiled source binding: $path")
    end
    reference=build["coupling_reference"]
    require(reference["contract_sha256"]==bytes2hex(sha256(read(Producer.CONTRACT))) &&
        reference["paper_sha256"]==C["primary_source"]["sha256"],"comparator source binding")
    require(bytes2hex(sha256(read(joinpath(dirname(buildfile),"paper.html"))))==C["primary_source"]["sha256"],"current locked primary source")
    require(d["scope"]=="reconstructed_theory" && d["status"]=="PASS" && !haskey(d,"error"),"invalid/incomplete report")
    for key in ("production_ready","physical_extraction","G2PERF_pass","historical_coupling_recovered",
        "physical_partial_moments_evaluated","failed_paper_benchmark_replaced")
        require(d[key]===false,"forged admission claim: $key")
    end
    require(d["contract_id"]==C["id"] && d["contract_sha256"]==bytes2hex(sha256(read(Producer.CONTRACT))),"contract binding")
    require(d["sources"]==Producer.source_hashes(),"source binding")
    require(d["primary_source"]==C["primary_source"] && d["policy"]==P && d["settings"]==C["settings"] && d["limitations"]==C["limitations"],"policy/source/limit mismatch")
    expected=Dict{String,Any}();references=Float64[]
    add(name,a,b,tol;absolute=false)=(expected[name]=(a=Float64.(vec(collect(a))),b=Float64.(vec(collect(b))),tol=tol,absolute=absolute))
    qs=Float64.(d["Q2_probes_GeV2"]);q0=C["settings"]["Q02_GeV2"]
    require(qs==vcat(Float64.(C["settings"]["Q2_probes_GeV2"]),C["settings"]["MZ_probe_GeV"]^2),"probe scales")
    batches=d["native_batches"]
    expected_batches=sum(length(level["quadrature_trials"]) for model in d["native_models"] for level in model["dt_trials"])+
        length(d["native_models"])+length(d["one_loop_controls"])+length(d["protocol_controls"])
    require(length(batches)==expected_batches,"missing/extra native batches")
    batch_values=Dict{Int,Any}()
    for (i,batch) in enumerate(batches)
        for key in ("request","response","stderr")
            path=batch[key]
            require(isfile(path) && !islink(path) && startswith(realpath(path),realpath(joinpath(RUN,"batches"))*"/"),"escaping batch path")
            require(bytes2hex(sha256(read(path)))==batch["sha256"][path],"native batch evidence hash")
        end
        require(count("TOO HIGH ACCURACY REQUIRED",read(batch["stderr"],String))==batch["precision_errors"],"diagnostic count")
        lines=split(chomp(read(batch["request"],String)),'\n')
        if batch["exit_code"]==0
            header=parse.(Int,split(lines[2]));parameters=parse.(Float64,split(lines[3]));points=parse.(Float64,lines[4:end])
            require(header[1]==length(points),"native request count")
            values=parse_native_response(read(batch["response"],String),header...)
            batch_values[i]=(nf=header[2],loops=header[3],alpha0=parameters[1],q0=parameters[2],dt=parameters[3],qs=points,values=values)
        else
            require(isempty(read(batch["response"],String)) && occursin("REFERENCE_INPUT_OR_RESULT_REJECTED",read(batch["stderr"],String)),"invalid-request control not rejected by comparator")
        end
    end
    function rowcheck(row,probes)
        batch=batch_values[row["batch_index"]]
        require(batch.qs==row["query_Q2_GeV2"] && batch.values==row["native_alpha"] && batch.alpha0==row["alpha0"] &&
            batch.q0==row["Q02_GeV2"] && batch.dt==row["dt"] && batch.nf==row["nf"] && batch.loops==row["nloop"],"native operands disagree with protocol bytes")
        positions=searchsortedfirst.(Ref(batch.qs),probes)
        require(batch.values[positions]==row["probe_alpha"],"probe/native binding")
        reference=Any[]
        if deep
            for q in batch.qs
                key=(batch.alpha0,q,batch.q0,batch.nf,batch.loops)
                push!(reference,get!(REF_CACHE,key) do
                    rg_reference(key...;bits=C["independent_reference"]["precision_bits"])
                end)
            end
            require([r.alpha for r in reference]==row["reference_alpha"],"independent alpha reference changed")
            require([reference[i].tau for i in positions]==row["reference_tau"],"independent time reference changed")
            append!(references,[r.implicit_residual for r in reference])
        end
        require(row["reference_alpha"][positions]==row["probe_reference_alpha"],"probe/reference binding")
        times=Float64[]
        for (q,rule) in zip(probes,row["rules"])
            if q==batch.q0
                require(isempty(rule["nodes"]) && isempty(rule["weights"]),"identity integral rule")
                push!(times,0.)
            else
                nodes,weights=Float64.(rule["nodes"]),Float64.(rule["weights"])
                expected_rule=log_rule(batch.q0,q,row["quadrature_order"])
                require(near(nodes,expected_rule.nodes) && near(weights,expected_rule.weights),"changed quadrature rule")
                require(length(nodes)==length(weights)==row["quadrature_order"] && all(v->batch.q0<v<q,nodes) && all(>(0),weights),"quadrature node/measure domain")
                require(isapprox(sum(weights),log(q/batch.q0)/(2pi);rtol=2e-14,atol=2e-15),"missing log-Q2/2pi measure")
                push!(times,dot(weights,batch.values[searchsortedfirst.(Ref(batch.qs),nodes)]))
            end
        end
        require(near(times,row["tau"]),"quadrature payload/time mismatch")
        require(isapprox(scaled_error(row["native_alpha"],row["reference_alpha"],P["relative_floor"]),row["sampled_alpha_error"];rtol=2e-13,atol=1e-15),"sampled alpha error")
        require(isapprox(scaled_error(row["probe_alpha"][1:1],[row["alpha0"]],P["relative_floor"]),row["anchor_error"];rtol=2e-13,atol=1e-15),"anchor error")
        require(isapprox(scaled_error(row["tau"],row["reference_tau"],P["relative_floor"]),row["time_reference_error"];rtol=2e-13,atol=1e-15),"time/reference error")
    end
    base=d["asymptotic_baseline"]
    lambda=C["settings"]["Lambda_GeV"];nf=C["settings"]["nf"];beta0=11-2nf/3;beta1=102-38nf/3
    alpha(q)=begin L=log(q/lambda^2);4pi/(beta0*L)*(1-beta1/beta0^2*log(L)/L) end
    tau(q)=begin L,L0=log(q/lambda^2),log(q0/lambda^2);(2/beta0)*(log(L/L0)+beta1/beta0^2*((log(L)+1)/L-(log(L0)+1)/L0)) end
    require(near(alpha.(qs),base["alpha"]) && near(tau.(qs),base["tau"]) && isapprox(base["alpha0"],alpha(q0);rtol=3e-15),"asymptotic definition binding")
    require(length(d["native_models"])==2,"missing native models")
    for (mi,model) in enumerate(d["native_models"])
        id=C["models"][mi+1]["id"]
        require(model["id"]==id && model["role"]==C["models"][mi+1]["role"],"native model identity")
        require(model["alpha0"]==(mi==1 ? C["settings"]["alpha0_literal"] : base["alpha0"]),"native initial condition")
        levels=model["dt_trials"]
        require(3<=length(levels)<=length(P["native_dt_refinements"]) && [v["dt"] for v in levels]==P["native_dt_refinements"][1:length(levels)],"dt ladder")
        for (i,level) in enumerate(levels)
            rows=level["quadrature_trials"]
            require(3<=length(rows)<=length(P["quadrature_orders"]) && [r["quadrature_order"] for r in rows]==P["quadrature_orders"][1:length(rows)],"quadrature ladder")
            for (j,row) in enumerate(rows)
                require(row["dt"]==level["dt"] && row["alpha0"]==model["alpha0"] && row["nf"]==nf && row["nloop"]==2,"trial setting identity")
                rowcheck(row,qs)
                j>1 && require(isapprox(row["quadrature_refinement_error"],scaled_error(row["tau"],rows[j-1]["tau"],P["relative_floor"]);rtol=2e-13,atol=1e-15),"quadrature refinement error")
            end
            final=rows[end]
            if i>1
                require(isapprox(level["dt_refinement_error"],scaled_error(final["probe_alpha"],levels[i-1]["quadrature_trials"][end]["probe_alpha"],P["relative_floor"]);rtol=2e-13,atol=1e-15),"dt refinement error")
            end
            accepted=i>=3 && final["sampled_alpha_error"]<=P["alpha_max_scaled"] && final["anchor_error"]<=P["anchor_max_scaled"] &&
                level["dt_refinement_error"]<=P["dt_refinement_max_scaled"] && final["quadrature_refinement_error"]<=P["quadrature_refinement_max_scaled"] && final["time_reference_error"]<=P["time_max_scaled"]
            require(level["accepted"]===accepted && level["status"]==(accepted ? "SELECTED_NUMERICAL_RESOLUTION" : "UNACCEPTED_REFINEMENT_TRIAL"),"forged coarse/final qualification")
            require(i==length(levels) || !accepted,"continued after accepted resolution")
        end
        selected=levels[end];f=selected["quadrature_trials"][end];previous=levels[end-1]["quadrature_trials"][end]
        require(model["selected_dt"]==selected["dt"] && model["selected_order"]==f["quadrature_order"] && model["alpha"]==f["probe_alpha"] && model["tau"]==f["tau"],"selected result binding")
        add("$id native alpha/reference",f["native_alpha"],f["reference_alpha"],P["alpha_max_scaled"])
        add("$id initial anchor",f["probe_alpha"][1:1],[model["alpha0"]],P["anchor_max_scaled"])
        add("$id final dt refinement",f["probe_alpha"],previous["probe_alpha"],P["dt_refinement_max_scaled"])
        add("$id native time/reference",f["tau"],f["reference_tau"],P["time_max_scaled"])
        add("$id final quadrature refinement",f["tau"],selected["quadrature_trials"][end-1]["tau"],P["quadrature_refinement_max_scaled"])
        add("$id analytic Mellin factors",[exp(transversity_moment(N)*t) for t in f["tau"],N in (1,2,3)],
            [exp(PN*t) for t in f["reference_tau"],PN in (-2/3,-2.,-26/9)],P["mellin_factor_max_scaled"])
        chain=model["composition_control"];rowcheck(chain,qs[2:3])
        add("$id native alpha composition",chain["probe_alpha"][2:2],f["probe_alpha"][3:3],P["composition_max_scaled"])
        add("$id native time composition",[f["tau"][2]+chain["tau"][2]],f["tau"][3:3],P["composition_max_scaled"])
        require(near(model["difference_from_asymptotic_alpha"],f["probe_alpha"]-base["alpha"]) && near(model["difference_from_asymptotic_tau"],f["tau"]-base["tau"]),"model difference binding")
    end
    require(length(d["one_loop_controls"])==4,"missing one-loop control")
    for (flavors,row) in zip(P["one_loop_nf_controls"],d["one_loop_controls"])
        require(row["nf"]==flavors && row["nloop"]==1 && row["dt"]==P["one_loop_control_dt"],"one-loop setting identity")
        rowcheck(row,qs)
        add("one-loop nf$flavors alpha",row["native_alpha"],row["reference_alpha"],P["one_loop_max_scaled"])
        add("one-loop nf$flavors time",row["tau"],row["reference_tau"],P["one_loop_max_scaled"])
    end
    residual_check=only(filter(c->c["name"]=="independent implicit-root residuals",d["checks"]))
    if deep
        require(residual_check["actual"]==[maximum(references)],"full high-precision residual operands")
        require(length(REF_CACHE)==d["independent_reference_cache_entries"],"independent reference point count")
    end
    add(residual_check["name"],residual_check["actual"],[0.],P["implicit_residual_max_abs"];absolute=true)
    require(Set(c["name"] for c in d["protocol_controls"])==Set(["wrong protocol","zero query count","unsupported nf","unsupported loops","zero dt","nonfinite alpha","backward query","truncated queries","unexpected trailer"]),"protocol control identities")
    require(length(d["protocol_controls"])==9 && length(unique(c["name"] for c in d["protocol_controls"]))==9,"missing/duplicate protocol controls")
    for control in d["protocol_controls"]
        require(control["rejected"]===true && batches[control["batch_index"]]["exit_code"]!=0,"invalid request admitted")
        add("native invalid request: "*control["name"],[0.],[0.],0.;absolute=true)
    end
    add("zero native precision diagnostics",[sum(b["precision_errors"] for b in batches)],[0.],0.;absolute=true)
    add("source/build stability",[0.],[0.],0.;absolute=true)
    require(length(expected)==length(d["checks"])==36 && length(unique(c["name"] for c in d["checks"]))==36,"missing/duplicate/extra check")
    for check in d["checks"]
        e=expected[check["name"]]
        require(near(e.a,check["actual"]) && near(e.b,check["reference"]),"check operands")
        require(check["tolerance"]==e.tol && check["relative_floor"]==P["relative_floor"] && check["metric"]==(e.absolute ? "absolute" : "scaled"),"altered metric/tolerance")
        difference=e.absolute ? maximum(abs.(e.a-e.b)) : scaled_error(e.a,e.b,P["relative_floor"])
        require(isapprox(check["error"],difference;rtol=2e-13,atol=1e-15) && check["passed"]===(isfinite(difference) && difference<=e.tol),"forged check result")
        require(check["passed"]===true,"failed component gate")
    end
    true
end

function verify_coupling_report(path,buildfile)
    d=JSON.parsefile(path)
    inspect_coupling_report(d,path,buildfile)
    d
end

function study_sources()
    paths=vcat(["config/coupling-moment-contract.json","config/kpsy-paper-moment-contract.json",
        "config/initial-contract.json","config/evolution-contract.json","config/parton-contract.json",
        "config/parton-sources.json","config/hoppet-contract.json","config/native-coupling-contract.json",
        "numerics/Project.toml","numerics/Manifest.toml","fortran/hoppet_bridge.f90",
        "fortran/native_coupling_reference.f90","src/collinsflow/hoppet.py","src/collinsflow/parton.py",
        "src/collinsflow/provenance.py","bridges/lhapdf.cpp","scripts/prepare_parton_inputs.py",
        "scripts/prepare_hoppet.py","scripts/prepare_native_coupling.py","scripts/prepare_kpsy_paper.py",
        "workflow/coupling-moments.smk"],
        ["numerics/bin/"*name for name in ("CouplingMomentReference.jl","PaperMoments.jl",
         "NativeCouplingReference.jl","CollinearReference.jl","check_native_coupling.jl",
         "check_prepared_times.jl","check_coupling_moments.jl")],
        ["numerics/src/"*name for name in readdir(joinpath(ROOT,"numerics/src")) if endswith(name,".jl")])
    Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths)
end

function moment_profiles(d,path)
    q0=d["settings"]["Q02_GeV2"];nf=d["settings"]["nf"];qs=d["Q2_probes_GeV2"]
    digest=bytes2hex(sha256(read(path)))
    study=JSON.parsefile(CONTRACT);source=JSON.parsefile(Producer.CONTRACT)
    models=vcat([d["asymptotic_baseline"]],d["native_models"])
    require([m["id"] for m in models]==[m["id"] for m in source["models"]],"coupling source model order/identity")
    # The study's descriptive baseline name maps only this exact source model.
    # Preserve the original identifier in source_id as well as in its receipt.
    [PreparedEvolutionTimes(qs,m["tau"],q0,nf,m["alpha0"],study["models"][i],
        "native-coupling-validation/"*m["id"],digest) for (i,m) in enumerate(models)]
end

"Independent direct powers and analytic first derivatives in high precision."
function initial_value_gradient(grid,soffer,theta,functional;bits=256)
    size(soffer)==(2,length(grid)) && length(functional)==length(grid) && length(theta)==13 || error("invalid independent derivative inputs")
    setprecision(bits) do
        gradient=zeros(BigFloat,13);value=BigFloat(0)
        for (f,weight) in enumerate((big"0.7",big"-0.4"))
            n,a,b=BigFloat(theta[f]),BigFloat(theta[f+2]),BigFloat(theta[f+4])
            normalization=(a+b)^(a+b)/(a^a*b^b)
            da,db=log((a+b)/a),log((a+b)/b)
            for (j,xvalue) in enumerate(grid)
                xvalue==1 && continue
                x=BigFloat(xvalue)
                shape=x^a*(1-x)^b*normalization*BigFloat(soffer[f,j])
                measure=weight*BigFloat(functional[j])
                value+=measure*n*shape
                gradient[f]+=measure*shape
                gradient[f+2]+=measure*n*shape*(log(x)+da)
                gradient[f+4]+=measure*n*shape*(log(1-x)+db)
            end
        end
        (value=value,gradient=gradient)
    end
end

function inspect_analytic_report(d,couplingfile,buildfile)
    contract=JSON.parsefile(CONTRACT)
    require(d["status"]=="PASS" && d["scope"]=="reconstructed_theory" && !haskey(d,"error"),"failed/incomplete prepared-time controls")
    for key in ("production_ready","physical_extraction","G2PERF_pass","historical_coupling_recovered","failed_paper_benchmark_replaced")
        require(d[key]===false,"forged analytic admission: $key")
    end
    require(d["contract_id"]==contract["prepared_time_interface"]["id"] &&
        d["contract_sha256"]==bytes2hex(sha256(read(CONTRACT))) && d["sources"]==study_sources(),"stale analytic sources/contract")
    require(d["build_sha256"]==bytes2hex(sha256(read(buildfile))) &&
        d["coupling_report_sha256"]==bytes2hex(sha256(read(couplingfile))),"analytic dependency binding")
    require(!isempty(d["checks"]) && length(d["checks"])==length(unique(c["name"] for c in d["checks"])),"missing/duplicate analytic controls")
    expected=analytic_check_policy(contract)
    coupling=JSON.parsefile(couplingfile);profiles=moment_profiles(coupling,couplingfile)
    require(Set(keys(expected))==Set(c["name"] for c in d["checks"]),"missing/extra analytic check")
    for check in d["checks"]
        metric,tolerance,floor=expected[check["name"]]
        if isnan(floor)
            profile=only(filter(p->startswith(check["name"],p.model_id*" native Mellin Q2="),profiles))
            q=parse(Float64,last(split(check["name"],"Q2=")))
            refs=[CollinearReference.mellin_polynomial_tau(x,evolution_time(profile,q,profile.Q02))
                for x in contract["analytic_policy"]["mellin_probes"]]
            require(near(refs,check["reference"]),"changed analytic Mellin reference")
            floor=contract["analytic_policy"]["mellin_floor_fraction_of_max"]*maximum(abs,refs)
        end
        require(check["metric"]==metric && check["tolerance"]==tolerance && check["relative_floor"]==floor,"altered analytic tolerance")
        a,b=Float64.(check["actual"]),Float64.(check["reference"])
        require(length(a)==length(b)>0 && all(isfinite,a) && all(isfinite,b),"absent/nonfinite analytic operands")
        difference=metric=="absolute" ? maximum(abs.(a-b)) : scaled_error(a,b,floor)
        require(isfinite(difference) && difference<=tolerance && check["passed"]===true &&
            isapprox(check["error"],difference;rtol=2e-13,atol=1e-15),"failed/forged analytic check")
    end
    d
end

verify_analytic_report(path,couplingfile,buildfile)=inspect_analytic_report(JSON.parsefile(path),couplingfile,buildfile)

function analytic_check_policy(contract)
    p=contract["analytic_policy"];f=p["relative_floor"];out=Dict{String,Tuple{String,Float64,Float64}}()
    add(name,metric,tol;floor=f)=(out[name]=(metric,tol,floor))
    for name in ("profile validation and unsupported API controls","source/build/dependency stability")
        add(name,"absolute",0.)
    end
    for model in contract["models"]
        for dy in p["native_dy"]
            add("$model native dy=$dy identity","absolute",p["identity_max_abs"])
            add("$model native dy=$dy endpoint","absolute",p["endpoint_max_abs"])
        end
        for q in p["Q2_GeV2"][2:end]
            add("$model native Mellin Q2=$q","scaled",p["mellin_max_scaled"];floor=NaN)
        end
        add("$model native composition","scaled",p["native_composition_max_scaled"])
        add("$model dense composition","scaled",p["dense_composition_max_scaled"])
        add("$model dense identity","absolute",p["identity_max_abs"])
        for backend in ("native","dense"), kind in ("gradient","hessian")
            add("$model $backend $kind","scaled",p["derivative_max_scaled"])
        end
    end
    for backend in ("native","dense")
        add("$backend unchanged default values","absolute",p["default_compatibility_max_abs"])
    end
    out
end

end
