# Analytic-input native observable reconstruction. No experimental admission.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, QuadGK
include("ObservableChecks.jl")

analytic_edges(B)=vcat([0.,1e-7,1e-6,1e-5,1e-4,.001,.005,.02,.05,.1,.2,.4,.7,1.,1.5,2.,3.,4.,6.,8.],filter(<(B),[12.,16.,20.,24.]),B)

function analytic_observable_checks(buildfile,progressfile)
    sources=observable_source_hashes(@__FILE__)
    build=CollinsFit.verified_hoppet_build(buildfile);build_hash=bytes2hex(sha256(read(buildfile)))
    contract=JSON.parsefile(OBS_CONTRACT);policy=contract["validation_policy"]
    parent=JSON.parsefile(joinpath(OBS_ROOT,contract["parent_contract"]))["validation_policy"]
    checks=Dict{String,Any}[]
    check(name,e,t)=observable_check!(checks,name,e,t)
    scaled(a,b)=observable_scaled(a,b,policy)
    theta=[.5,-.4,1.,1.,2.,2.,.3,-.2,1.,1.,0.,0.,.01];model=KPSYCSS()
    levels=Dict{String,Any}();fine=Dict{Symbol,Any}();coarse=Dict{Symbol,Any}()
    oracles=[];momentum=Dict{String,Any}();tails=Dict{String,Any}();derivatives=[]
    for process in (:SIDIS,:SIA)
        label=string(process);levels[label]=Any[]
        for (level,setting) in enumerate(parent["refinements"])
            radial=RadialRule(analytic_edges(24.);order=setting["radial_order"])
            context=analytic_context(buildfile,process,model,radial;intervals=setting["intervals"])
            values=all_observable_values(context,theta)
            source_values=angular_source_level(context,theta,values,level,policy,checks)
            error=isempty(levels[label]) ? nothing : scaled(values,last(levels[label])["values"])
            accepted=level>=3 && error!==nothing && error<=policy["physical_value_max_scaled"]
            push!(levels[label],Dict("setting"=>setting,"radial_nodes"=>length(radial.b),"values"=>values,
                "previous_level_max_error"=>error,"refinement_accepted"=>accepted,"mean_of_angular_ratio_values"=>source_values))
            coarse[process]=get(coarse,process,context);fine[process]=context
            writejson(progressfile,Dict("status"=>"IN_PROGRESS","scope"=>"reconstructed_theory","production_ready"=>false,
                "sources"=>sources,"refinements"=>levels))
            println(JSON.json(Dict("completed_process"=>label,"level"=>level,"radial_nodes"=>length(radial.b),"error"=>error,"accepted"=>accepted)));flush(stdout)
            accepted && break
        end
        check("$process complete analytic observable refinement admitted",last(levels[label])["refinement_accepted"] ? 0. : 1.,0.)
        check("$process adjacent radial/evolution observable difference",last(levels[label])["previous_level_max_error"],policy["physical_value_max_scaled"])
        context=fine[process]
        amplitude_indices=[1,2,7,8,13];a=theta[amplitude_indices]
        for interval in (nothing,OBS_BINS...)
            pair=observable_pair(context;interval=interval)
            oracle=mellin_observable_reference(context,theta;momentum_interval=interval,y=OBS_Y,angular=OBS_ANGULAR)
            function amplitude_prediction(aa)
                t=[(i in amplitude_indices ? aa[something(findfirst(==(i),amplitude_indices))] : zero(aa[1])+theta[i]) for i in 1:13]
                observable_values(pair,t)
            end
            actual=observable_values(pair,theta)
            J=ForwardDiff.jacobian(amplitude_prediction,a)
            H=cat([ForwardDiff.hessian(x->amplitude_prediction(x)[i],a) for i in 1:2]...;dims=3)
            errors=[scaled(actual,oracle.values),scaled(J,oracle.gradient),scaled(H,oracle.hessian)]
            check("$process interval $interval independent Mellin/adaptive observable",errors[1],policy["physical_value_max_scaled"])
            check("$process interval $interval amplitude/gc gradient moments",errors[2],policy["derivative_max_scaled"])
            check("$process interval $interval amplitude/gc Hessian moments",errors[3],policy["derivative_max_scaled"])
            push!(oracles,Dict("process"=>label,"momentum_bounds_GeV"=>interval,"definitions"=>OBS_NAMES[process],
                "values"=>actual,"reference_values"=>oracle.values,"gradient"=>J,"reference_gradient"=>oracle.gradient,
                "hessian"=>H,"reference_hessian"=>oracle.hessian,"amplitude_indices"=>amplitude_indices,
                "reference_charge_moments"=>oracle.charge_moments,"errors"=>errors))
        end
        momentum[label]=momentum_refinements(context,theta,policy,checks)
        # Gaussian Hankel identity is a mathematical integration check only.
        for interval in OBS_BINS
            weights=CollinsFit.observable_weights(context,MomentumRule(interval...;order=128))
            A=.7;plan=context.leptonic;z=plan.z1;n=process==:SIDIS ? 1 : 2
            prefactor=process==:SIDIS ? -inv(2z^3) : inv(4z^3*plan.z2)
            gaussian=exp.(-A.*plan.radial.b.^2)
            expectedUU=quadgk(p->p*PolarizedReference.LeptonicReference.gaussian_hankel(0,p/z,A)/z^2,interval...;rtol=1e-12,atol=1e-14)[1]
            expectedPol=quadgk(p->p*prefactor*PolarizedReference.LeptonicReference.gaussian_hankel(n,p/z,A),interval...;rtol=1e-12,atol=1e-14)[1]
            check("$process interval $interval Gaussian projected UU and Collins integrals",scaled([dot(weights[1],gaussian),dot(weights[2],gaussian)],[expectedUU,expectedPol]),policy["analytic_integral_max_scaled"])
        end
        zero_pair=observable_pair(with_momentum(context,0.))
        check("$process zero-P observable exact zero",maximum(abs.(observable_values(zero_pair,theta))),0.)
        tail=[]
        setting=last(levels[label])["setting"]
        for B in parent["radial_bounds"]
            tailcontext=B==24 ? context : analytic_context(buildfile,process,model,
                RadialRule(analytic_edges(B);order=setting["radial_order"]);intervals=setting["intervals"])
            push!(tail,Dict("B_inverse_GeV"=>B,"values"=>all_observable_values(tailcontext,theta)))
        end
        tails[label]=tail
        check("$process finite-bound observable extension 24 to 32",scaled(tail[end]["values"],tail[end-1]["values"]),policy["physical_value_max_scaled"])
        append!(derivatives,observable_derivatives(coarse[process],theta,policy,checks))
    end
    source_checks=angular_source_checks(fine[:SIA],coarse[:SIA],theta,oracles,derivatives,contract,checks;
        reference_tolerance=policy["physical_value_max_scaled"],reference_scope="independent Mellin/adaptive values at the same fixed effective angular factor")
    check("source/build stability",observable_source_hashes(@__FILE__)==sources && bytes2hex(sha256(read(buildfile)))==build_hash &&
        CollinsFit.verified_hoppet_build(buildfile)==build ? 0. : 1.,0.)
    observable_result(checks,sources,contract;build=build,build_sha256=build_hash,refinements=levels,tails=tails,
        momentum_refinements=momentum,independent_mellin_references=oracles,derivative_checks=derivatives,
        angular_source_checks=source_checks,parameters=theta,input_scope="analytic polynomial initial and unpolarized inputs; native matched evolution; no physical source replacement")
end

buildfile=isempty(ARGS) ? joinpath(OBS_ROOT,".runtime/hoppet/build.json") : ARGS[1]
output=length(ARGS)<2 ? joinpath(OBS_ROOT,"runs/observables/analytic.json") : ARGS[2]
try
    report=analytic_observable_checks(buildfile,output*".progress.json");writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"scope"=>report["scope"],"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
        "physical_extraction"=>false,"sources"=>observable_source_hashes(@__FILE__),"error"=>sprint(showerror,exception)))
    rethrow()
end
