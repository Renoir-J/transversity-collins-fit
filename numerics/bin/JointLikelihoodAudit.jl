module JointLikelihoodAudit
using CollinsFit, JSON, SHA, LinearAlgebra
include("JointLikelihoodReference.jl")
using .JointLikelihoodReference: augmented_reference

const ROOT=normpath(joinpath(@__DIR__,"../.."))
require(condition,message)=condition || error(message)
matrix(value)=hcat((Float64.(column) for column in value)...)
function close(actual,expected,tolerance=1e-11;floor=1e-8)
    a=vec(Float64.(actual));b=vec(Float64.(expected))
    require(length(a)==length(b) && all(isfinite,a) && all(isfinite,b),"invalid numerical operand")
    require(maximum(abs.(a-b)./max.(abs.(b),floor))<=tolerance,"numerical operand differs")
end

function audit_joint_report(report)
    contract=JSON.parsefile(joinpath(ROOT,"config/joint-likelihood-contract.json"))
    input_contract=JSON.parsefile(joinpath(ROOT,"config/joint-pinned-input-contract.json"))
    policy=contract["validation_policy"]
    require(report["status"]=="PASS" && !haskey(report,"error"),"failed or incomplete producer")
    require(report["scope"]=="synthetic_closure" && report["theory_scope"]=="reconstructed_theory","wrong scientific scope")
    for key in ("production_ready","physical_extraction","G2PERF_pass","real_data_fit_implemented","experimental_data_admitted")
        require(report[key]===false,"invalid admission flag")
    end
    require(report["contract_id"]==contract["id"] && report["policy"]==policy,"changed contract or acceptance policy")
    require(report["contract_sha256"]==bytes2hex(sha256(read(joinpath(ROOT,"config/joint-likelihood-contract.json")))) &&
        report["input_contract_sha256"]==bytes2hex(sha256(read(joinpath(ROOT,"config/joint-pinned-input-contract.json")))),"stale contract")
    expected=vcat(["numerics/bin/check_joint_likelihood.jl","numerics/bin/JointLikelihoodReference.jl","numerics/bin/PaperMoments.jl",
        "config/joint-likelihood-contract.json","config/joint-pinned-input-contract.json","numerics/Project.toml","numerics/Manifest.toml"],
        ["config/"*n*"-contract.json" for n in ("initial","evolution","hoppet","parton","polarized-process",
            "leptonic-observable","leptonic-bin","sidis-convention","star-evolved","star-born")],
        ["numerics/src/"*f for f in readdir(joinpath(ROOT,"numerics/src")) if endswith(f,".jl")])
    require(Set(keys(report["sources"]))==Set(expected),"incomplete computational source inventory")
    for path in expected
        require(report["sources"][path]==bytes2hex(sha256(read(joinpath(ROOT,path)))),"stale computational source: $path")
    end
    kind=report["kind"]
    require(kind in ("analytic","inputs"),"unknown producer kind")
    dependencies=report["dependency_sha256"]
    if kind=="inputs"
        require(report["input_kind"]=="pinned_CT10_DSSV08_DSS14_native_rows" &&
            report["input_preparation_policy"]==input_contract && length(dependencies)==2,"wrong physical input policy")
        for (path,sha) in dependencies
            require(bytes2hex(sha256(read(path)))==sha,"changed native/input receipt")
            receipt=JSON.parsefile(path)
            if receipt==report["native_build"]
                require(CollinsFit.verified_hoppet_build(path)==receipt,"invalid native build")
            elseif receipt==report["parton_build"]
                require(CollinsFit.verify_parton_build(path)==receipt,"invalid parton build")
            else
                error("unbound input receipt")
            end
        end
        require(report["native_build"]!=report["parton_build"],"distinct native and parton inputs required")
    else
        require(report["input_kind"]=="explicitly_synthetic_source_rows" && isempty(dependencies),"wrong synthetic input policy")
    end
    require(report["matrix_encoding"]=="outer arrays are Julia matrix columns","unknown matrix encoding")
    d=report["data"];p=report["profile"];derivative=report["derivatives"]
    ids=["synthetic:$i" for i in 1:8]
    require(d["point_ids"]==ids && report["model"]["point_ids"]==ids && report["likelihood"]["point_ids"]==ids,"point order changed")
    require(report["model"]["processes"]==["SIDIS","SIDIS","SIA","SIA","SIDIS","STAR_200","STAR_510","SIA"],"process order changed")
    require(report["model"]["scope"]=="reconstructed_theory" && report["likelihood"]["scope"]=="synthetic_closure","child scope changed")
    for metadata in (report["model"],report["likelihood"]),key in ("production_ready","physical_extraction","G2PERF_pass")
        require(metadata[key]===false,"child admission changed")
    end
    t=Float64.(p["predictions"]);y=Float64.(d["observations"])
    C=matrix(d["covariance"]);A=matrix(d["additive"]);R=matrix(d["relative"]);S=matrix(d["nuisance_prior"])
    require(length(t)==8 && size(C)==(8,8) && size(A)==size(R)==(8,2) && size(S)==(2,2),"wrong fixture dimensions")
    close(y,t+[.01cos(i) for i in 1:8],0.)
    sigma=[.02+.001i for i in 1:8]
    close(C,[sigma[i]*sigma[j]*.3^abs(i-j) for i in 1:8,j in 1:8],0.)
    close(A,[.001sin(i*j) for i in 1:8,j in 1:2],0.)
    close(R,[.02cos(.3i+j) for i in 1:8,j in 1:2],0.)
    close(S,[1. .35;.35 1.4],0.)
    require(d["nuisance_ids"]==["synthetic shared 1","synthetic shared 2"]==report["likelihood"]["nuisance_ids"],"nuisance order changed")
    theta=collect(kpsy_published_parameters().values);theta[2]=-.8
    for i in (3,4,5,6,9,10,11,12);theta[i]+=.2;end
    close(report["parameters"],theta,0.)
    reference=setprecision(256) do
        augmented_reference(BigFloat.(t),y,C,A,R,S)
    end
    close([p["chi2"]],[reference.chi2]);close(p["nuisance"],reference.nuisance)
    close([p["independent_QR_chi2"]],[reference.chi2]);close(p["independent_QR_nuisance"],reference.nuisance)
    residual=y-t;design=A+reshape(t,:,1).*R;V=C+design*S*transpose(design)
    nuisance=Float64.(p["nuisance"]);shifted=residual-design*nuisance
    close(p["residual"],residual,0.);close(p["shifted_residual"],shifted,0.)
    close([p["chi2"]],[dot(shifted,C\shifted)+dot(nuisance,S\nuisance)])
    counter=report["counterfactual_values"]
    data_design=A+reshape(y,:,1).*R
    close([counter["profile"]],[p["chi2"]],0.)
    close([counter["data_scaled"]],[dot(residual,(C+data_design*S*transpose(data_design))\residual)])
    diagonal=augmented_reference(t,y,Matrix(Diagonal(diag(C))),A,R,S).chi2
    close([counter["diagonalized_C"]],[diagonal])
    logdet_difference=logdet(Symmetric(V))-logdet(Symmetric(C))
    close([counter["marginal_relative_logdet"]],[logdet_difference])
    direction=Float64.(derivative["direction"])
    close(direction,Float64.(input_contract["derivative_direction"]),0.)
    trials=derivative["finite_difference_trials"]
    require(length(trials)==2 && [r["step"] for r in trials]==policy["finite_difference_steps"],"changed difference steps")
    gradient,hvp=setprecision(policy["reference_bits"]) do
        gradients=Vector{BigFloat}[];products=Vector{BigFloat}[]
        for row in trials
            h=BigFloat(row["step"]);gv=row["gradient_function_values_decimal"];hv=row["hvp_function_values_decimal"]
            require(length(gv)==length(hv)==13 && all(length(v)==2 for v in gv) && all(length(v)==4 for v in hv),"missing raw FD operands")
            parse_values(v)=parse.(BigFloat,v)
            g=[begin v=parse_values(values);require(all(isfinite,v),"nonfinite FD operand");(v[1]-v[2])/(2h) end for values in gv]
            q=[begin v=parse_values(values);require(all(isfinite,v),"nonfinite FD operand");(v[1]-v[2]-v[3]+v[4])/(4h^2) end for values in hv]
            close(g,parse_values(row["gradient_decimal"]),1e-60)
            close(q,parse_values(row["hvp_decimal"]),1e-60)
            push!(gradients,g);push!(products,q)
        end
        (Float64.((4gradients[2]-gradients[1])/3),Float64.((4products[2]-products[1])/3))
    end
    close(derivative["reference_gradient"],gradient,0.);close(derivative["reference_hvp"],hvp,0.)
    H=matrix(derivative["hessian"])
    require(size(H)==(13,13),"incomplete Hessian")
    close(derivative["hvp"],H*direction,0.)
    controls=report["negative_controls"]
    control_names=vcat(["reject production","reject experimental_candidate","reject historical_benchmark_candidate",
        "wrong point order","unsupported lattice recipe","duplicate source phase"],["mutated $f" for f in
            (:observations,:covariance,:additive,:relative,:nuisance_prior,:nuisance_factor)])
    require([c["name"] for c in controls]==control_names && all(c["rejected"]===true && !isempty(c["error"]) for c in controls),"missing rejection controls")
    # Recompute the comparison ledger independently; reported PASS bits are insufficient.
    checks=report["checks"];seen=Set{String}()
    absolute=Set(["subset retains covariance","subset retains correlated nuisance prior","data-scaled nuisance is different",
        "discarded correlations are different","marginal logdet is not profile chi2","all scope identity and mutation controls reject",
        "in-place gradient agreement","in-place value agreement","sources and build dependencies stable"])
    names=Set(vcat(collect(absolute),["independent augmented QR profile","independent augmented QR nuisance",
        "global transversity Collins sign ambiguity preserved","global sign ambiguity in joint profile","effective covariance identity",
        "original-coordinate nuisance stationarity","zero nuisance GLS reduction","nuisance-coordinate invariant profile",
        "nuisance-coordinate invariant physical pulls","reordered marginal subset profile","process holdout marginal profile",
        "complete 13-coordinate gradient","complete 13-coordinate HVP"]))
    require(length(checks)==length(names),"incomplete comparison ledger")
    indices=[7,1,3];keep=[1,2,5,6,7]
    sub=augmented_reference(t[indices],y[indices],C[indices,indices],A[indices,:],R[indices,:],S)
    hold=augmented_reference(t[keep],y[keep],C[keep,keep],A[keep,:],R[keep,:],S)
    operands=Dict{String,Any}(
        "independent augmented QR profile"=>[reference.chi2],"independent augmented QR nuisance"=>reference.nuisance,
        "global transversity Collins sign ambiguity preserved"=>t,"global sign ambiguity in joint profile"=>[p["chi2"]],
        "effective covariance identity"=>[dot(residual,V\residual)],
        "original-coordinate nuisance stationarity"=>S\nuisance,"zero nuisance GLS reduction"=>[dot(residual,C\residual)],
        "nuisance-coordinate invariant profile"=>[p["chi2"]],"nuisance-coordinate invariant physical pulls"=>nuisance,
        "reordered marginal subset profile"=>[sub.chi2],"process holdout marginal profile"=>[hold.chi2],
        "subset retains covariance"=>C[indices,indices],"subset retains correlated nuisance prior"=>S,
        "complete 13-coordinate gradient"=>gradient,"complete 13-coordinate HVP"=>hvp,
        "in-place gradient agreement"=>derivative["gradient"],"in-place value agreement"=>[p["chi2"]])
    require(abs(counter["data_scaled"]-p["chi2"])>policy["relative_floor"] &&
        abs(diagonal-p["chi2"])>policy["relative_floor"] && abs(logdet_difference)>policy["relative_floor"],"ineffective counterfactual fixture")
    for name in setdiff(absolute,Set(keys(operands)));operands[name]=[0.];end
    for check in checks
        name=check["name"];require(name in names && !(name in seen),"unknown or duplicate comparison");push!(seen,name)
        is_absolute=name in absolute;tol=is_absolute ? 0. : startswith(name,"complete 13-coordinate") ? policy["derivative_max_scaled"] : policy["algebra_max_scaled"]
        require(check["metric"]==(is_absolute ? "absolute" : "scaled") && check["tolerance"]==tol &&
            check["relative_floor"]==policy["relative_floor"],"comparison policy changed")
        a=vec(Float64.(check["actual"]));b=vec(Float64.(check["reference"]))
        require(length(a)==length(b) && !isempty(a) && all(isfinite,a) && all(isfinite,b),"invalid comparison values")
        error=maximum(is_absolute ? abs.(a-b) : abs.(a-b)./max.(abs.(b),policy["relative_floor"]))
        require(check["error"]==error && check["passed"]===true && error<=tol,"failed comparison or forged PASS")
        close(b,operands[name],max(tol,1e-11));close(a,operands[name],max(tol,1e-11))
    end
    close(only(c for c in checks if c["name"]=="complete 13-coordinate gradient")["actual"],derivative["gradient"],0.)
    close(only(c for c in checks if c["name"]=="complete 13-coordinate gradient")["reference"],gradient,0.)
    close(only(c for c in checks if c["name"]=="complete 13-coordinate HVP")["actual"],derivative["hvp"],0.)
    close(only(c for c in checks if c["name"]=="complete 13-coordinate HVP")["reference"],hvp,0.)
    Dict("status"=>"PASS","scope"=>"synthetic_closure_report_audit","theory_scope"=>"reconstructed_theory",
        "kind"=>kind,"production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,
        "comparison_count"=>length(checks),"native_input_runtime_revalidated"=>(kind=="inputs"),
        "process_predictions_replayed"=>false,"finite_difference_function_values_replayed"=>false,
        "limitations"=>["Recomputes retained profile/FD algebra and checks current sources/builds; does not replay process predictions or FD function calls",
            "Synthetic observations remain synthetic with pinned physical inputs; no accuracy, fit, coverage or production admission"])
end
end
