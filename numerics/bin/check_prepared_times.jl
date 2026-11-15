# Analytic/native controls only. No physical input data or production admission.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff
include("CouplingMomentReference.jl")
include("PaperMoments.jl")
include("CollinearReference.jl")
using .CouplingMomentReference
using .PaperMoments: compare_check
using .CollinearReference: validation_grid, mellin_polynomial_tau

function prepared_checks!(report,buildfile,couplingfile)
    contract=JSON.parsefile(CouplingMomentReference.CONTRACT);p=contract["analytic_policy"]
    sources=study_sources();buildsha=bytes2hex(sha256(read(buildfile)));couplingsha=bytes2hex(sha256(read(couplingfile)))
    coupling=verify_coupling_report(couplingfile,buildfile);profiles=moment_profiles(coupling,couplingfile)
    report["contract_id"]=contract["prepared_time_interface"]["id"]
    report["contract_sha256"]=bytes2hex(sha256(read(CouplingMomentReference.CONTRACT)))
    report["sources"]=sources;report["build_sha256"]=buildsha;report["coupling_report_sha256"]=couplingsha
    report["policy"]=p;report["limitations"]=contract["limitations"]
    report["profiles"]=prepared_time_metadata.(profiles)
    qs=Float64.(p["Q2_GeV2"]);q0=first(qs);theta=[.2,-.3];checks=report["checks"]
    check(name,a,b,tol;metric=:scaled,floor=p["relative_floor"])=push!(checks,compare_check(name,a,b,tol;metric=metric,floor=floor))
    function reject(label,f)
        message="";rejected=false
        try;f();catch exception;rejected=true;message=sprint(showerror,exception);end
        push!(report["negative_controls"],Dict("name"=>label,"rejected"=>rejected,"error"=>message))
        rejected
    end
    controls=Bool[]
    for (label,mutate) in [
        ("coupling receipt status only",r->(r["status"]="IN_PROGRESS")),
        ("coupling missing check",r->pop!(r["checks"])),
        ("coupling altered tolerance",r->(r["checks"][1]["tolerance"]*=2)),
        ("coupling forged operand",r->(r["checks"][1]["actual"][1]=999.)),
        ("coupling stale source",r->(r["sources"][first(keys(r["sources"]))]=repeat("0",64))),
        ("coupling altered selected time",r->(r["native_models"][1]["tau"][2]+=.01)),
        ("coupling production claim",r->(r["production_ready"]=true)),
        ("coupling promoted coarse result",r->(r["native_models"][1]["dt_trials"][1]["accepted"]=true))]
        forged=deepcopy(coupling);mutate(forged)
        push!(controls,reject(label,()->inspect_coupling_report(forged,couplingfile,buildfile;deep=false)))
    end
    profile=first(profiles);meta=prepared_time_metadata(profile)
    for (label,f) in [
        ("unavailable exact scale",()->prepared_times(profile,[nextfloat(10.)],q0)),
        ("wrong initial scale",()->prepared_times(profile,[10.],3.)),
        ("wrong nf",()->prepared_times(profile,[10.],q0;nf=4)),
        ("empty request",()->prepared_times(profile,Float64[],q0)),
        ("scale AD",()->ForwardDiff.derivative(q->evolution_time(profile,q,q0),10.)),
        ("time AD",()->ForwardDiff.derivative(t->PreparedEvolutionTimes(qs,[0.,t,.2],q0,5,.327,"test","test",meta.source_sha256).alpha0,.06)),
        ("duplicate labels",()->PreparedEvolutionTimes([q0,q0],[0.,.1],q0,5,.327,"test","test",meta.source_sha256)),
        ("missing zero anchor",()->PreparedEvolutionTimes([q0,10.],[.001,.1],q0,5,.327,"test","test",meta.source_sha256)),
        ("incorrect sign",()->PreparedEvolutionTimes([q0,10.],[0.,-.1],q0,5,.327,"test","test",meta.source_sha256)),
        ("nonfinite time",()->PreparedEvolutionTimes([q0,10.],[0.,NaN],q0,5,.327,"test","test",meta.source_sha256)),
        ("invalid source digest",()->PreparedEvolutionTimes([q0,10.],[0.,.1],q0,5,.327,"test","test","missing")),
        ("point row unsupported",()->HoppetPointPlan("absent",[.1],[10.],q0,profile)),
        ("matching unsupported",()->HoppetMatchingPlan("absent",[.1],[10.],q0,profile;process=:SIDIS)),
        ("no alpha reconstruction",()->alpha_s(profile,10.))]
        push!(controls,reject(label,f))
    end
    ownedq,ownedt=collect(profile.Q2),collect(profile.tau)
    owned=PreparedEvolutionTimes(ownedq,ownedt,q0,5,profile.alpha0,profile.model_id,profile.source_id,profile.source_sha256)
    ownedq[2]=11.;ownedt[2]=.07
    push!(controls,owned.Q2==profile.Q2 && owned.tau==profile.tau && !ismutabletype(typeof(owned)))
    other=PreparedEvolutionTimes(profile.Q2,profile.tau,q0,5,profile.alpha0,profile.model_id,"other-source",profile.source_sha256)
    push!(controls,CollinsFit.numerical_snapshot(other)!=CollinsFit.numerical_snapshot(profile))
    native_default=nothing;dense_default=nothing
    for (mi,profile) in enumerate(profiles)
        id=profile.model_id;fine=nothing;finegrid=nothing
        for dy in p["native_dy"]
            plan=HoppetPlan(buildfile,qs,q0,profile;nf=profile.nf,xmin=p["native_xmin"],dy=dy,
                layout=:nested,order=p["native_order"],eps=p["native_eps"],du=p["native_du"])
            try
                x=plan.grid;h=x.*(1 .-x).^3
                result=apply_operator(plan,h)
                values=dropdims(hoppet_evaluate(plan,reshape(result,length(x),1,length(qs)),p["mellin_probes"]);dims=2)
                check("$id native dy=$dy identity",result[:,1],h,p["identity_max_abs"];metric=:absolute)
                check("$id native dy=$dy endpoint",result[end,:],zeros(length(qs)),p["endpoint_max_abs"];metric=:absolute)
                push!(report["native_refinements"],Dict("model_id"=>id,"dy"=>dy,"grid"=>x,"initial"=>h,
                    "evolved"=>result,"probe_values"=>values,"tau"=>plan.tau,"signature"=>plan.signature))
                if dy==last(p["native_dy"])
                    fine=result;finegrid=copy(x)
                    for (k,q) in enumerate(qs[2:end])
                        reference=[mellin_polynomial_tau(xx,plan.tau[k+1]) for xx in p["mellin_probes"]]
                        floor=p["mellin_floor_fraction_of_max"]*maximum(abs,reference)
                        check("$id native Mellin Q2=$q",values[:,k+1],reference,p["mellin_max_scaled"];floor=floor)
                    end
                    amplitude=t->exp(t[1])+t[2]^2+t[1]*t[2]
                    objective=t->sum(apply_operator(plan,amplitude(t).*h)[:,end])
                    coefficient=sum(result[:,end])
                    check("$id native gradient",ForwardDiff.gradient(objective,theta),coefficient.*[exp(theta[1])+theta[2],2theta[2]+theta[1]],p["derivative_max_scaled"])
                    check("$id native hessian",ForwardDiff.hessian(objective,theta),coefficient.*[exp(theta[1]) 1.;1. 2.],p["derivative_max_scaled"])
                    saved=plan.tau[2];plan.tau[2]=nextfloat(saved)
                    push!(controls,reject("$id mutated native tau",()->apply_operator(plan,h)));plan.tau[2]=saved
                    saved=plan.grid[2];plan.grid[2]=nextfloat(saved)
                    push!(controls,reject("$id mutated native grid",()->apply_operator(plan,h)));plan.grid[2]=saved
                    mi==1 && (native_default=copy(result))
                end
            finally;close(plan);end
        end
        shifted=PreparedEvolutionTimes(qs[2:3],[0.,profile.tau[3]-profile.tau[2]],qs[2],profile.nf,
            (mi==1 ? coupling["asymptotic_baseline"] : coupling["native_models"][mi-1])["alpha"][2],
            id*"-shifted",profile.source_id,profile.source_sha256)
        plan=HoppetPlan(buildfile,[qs[3]],qs[2],shifted;nf=profile.nf,xmin=p["native_xmin"],dy=last(p["native_dy"]),
            layout=:nested,order=p["native_order"],eps=p["native_eps"],du=p["native_du"])
        try
            plan.grid==finegrid || error("composition grid mismatch")
            composed=vec(apply_operator(plan,fine[:,2]))
            check("$id native composition",composed,fine[:,3],p["native_composition_max_scaled"])
        finally;close(plan);end
        x=validation_grid(p["dense_nodes"],p["native_xmin"]);h=x.*(1 .-x).^3
        identity=collinear_lo_operator(x,q0,q0,profile)
        op=collinear_lo_operator(x,qs[3],q0,profile)
        direct=apply_operator(op,h)
        composed=apply_operator(collinear_lo_operator(x,qs[3],qs[2],shifted),
            apply_operator(collinear_lo_operator(x,qs[2],q0,profile),h))
        check("$id dense composition",composed,direct,p["dense_composition_max_scaled"])
        check("$id dense identity",apply_operator(identity,h),h,p["identity_max_abs"];metric=:absolute)
        amplitude=t->exp(t[1])+t[2]^2+t[1]*t[2]
        objective=t->sum(apply_operator(op,amplitude(t).*h))
        coefficient=sum(direct)
        check("$id dense gradient",ForwardDiff.gradient(objective,theta),coefficient.*[exp(theta[1])+theta[2],2theta[2]+theta[1]],p["derivative_max_scaled"])
        check("$id dense hessian",ForwardDiff.hessian(objective,theta),coefficient.*[exp(theta[1]) 1.;1. 2.],p["derivative_max_scaled"])
        push!(report["dense_controls"],Dict("model_id"=>id,"grid"=>x,"initial"=>h,"direct"=>direct,"composed"=>composed))
        mi==1 && (dense_default=copy(direct))
        saved=op.matrix[1,2];op.matrix[1,2]=nextfloat(saved)
        push!(controls,reject("$id mutated dense matrix",()->apply_operator(op,h)))
        op.matrix[1,2]=saved;op.grid[2]=nextfloat(op.grid[2])
        push!(controls,reject("$id mutated dense grid",()->apply_operator(op,h)))
        println(JSON.json(Dict("stage"=>"prepared-time-controls","model"=>id,"checks"=>length(checks))));flush(stdout)
    end
    original=AsymptoticTwoLoopCoupling()
    plan=HoppetPlan(buildfile,qs,q0,original;xmin=p["native_xmin"],dy=last(p["native_dy"]),
        layout=:nested,order=p["native_order"],eps=p["native_eps"],du=p["native_du"])
    try
        h=plan.grid.*(1 .-plan.grid).^3
        check("native unchanged default values",apply_operator(plan,h),native_default,0.;metric=:absolute)
        push!(controls,plan.time_profile===nothing && plan.signature==CollinsFit.numerical_snapshot((plan.grid,plan.Q2,plan.tau)))
    finally;close(plan);end
    x=validation_grid(p["dense_nodes"],p["native_xmin"])
    check("dense unchanged default values",apply_operator(collinear_lo_operator(x,qs[3],q0,original),x.*(1 .-x).^3),dense_default,0.;metric=:absolute)
    check("profile validation and unsupported API controls",[all(controls) ? 0. : 1.],[0.],0.;metric=:absolute)
    stable=study_sources()==sources && bytes2hex(sha256(read(buildfile)))==buildsha &&
        bytes2hex(sha256(read(couplingfile)))==couplingsha && verify_coupling_report(couplingfile,buildfile)==coupling
    check("source/build/dependency stability",[stable ? 0. : 1.],[0.],0.;metric=:absolute)
    report["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
    CouplingMomentReference.inspect_analytic_report(report,couplingfile,buildfile)
    report
end

function main()
    length(ARGS)==3 || error("usage: check_prepared_times.jl NATIVE_BUILD COUPLING_VALIDATION OUTPUT_JSON")
    report=Dict{String,Any}("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
        "physical_extraction"=>false,"G2PERF_pass"=>false,"historical_coupling_recovered"=>false,
        "failed_paper_benchmark_replaced"=>false,"checks"=>Dict{String,Any}[],"negative_controls"=>Dict{String,Any}[],
        "native_refinements"=>Dict{String,Any}[],"dense_controls"=>Dict{String,Any}[],"julia_version"=>string(VERSION))
    try;prepared_checks!(report,ARGS[1],ARGS[2]);catch exception;report["status"]="FAIL";report["error"]=sprint(showerror,exception,catch_backtrace());
    finally;report["completed_at"]=string(now(UTC));writejson(ARGS[3],report);end
    println(JSON.json(Dict("status"=>report["status"],"checks"=>length(report["checks"]),"error"=>get(report,"error",nothing),
        "failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
end

if abspath(PROGRAM_FILE)==(@__FILE__);main();end
