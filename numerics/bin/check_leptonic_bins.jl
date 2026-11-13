# Complete Julia bin-composition validation job; never a fit or admission gate.
using CollinsFit, JSON, SHA, Dates, ForwardDiff, LinearAlgebra, QuadGK
include("ObservableChecks.jl")
include("BinReference.jl")

const BIN_CONTRACT=joinpath(OBS_ROOT,"config/leptonic-bin-contract.json")
bin_sources()=merge(observable_source_hashes(@__FILE__),Dict(
    relpath(p,OBS_ROOT)=>bytes2hex(sha256(read(p))) for p in (BIN_CONTRACT,joinpath(@__DIR__,"BinReference.jl"))))

function validation_bin(nodes,weights,measure)
    leptonic_bin(nodes,weights;source="explicit reconstructed bin-composition diagnostic; no experimental admission",
        measure=measure,weight_definition="fixed quadrature, cross-section prefactor and explicitly synthetic acceptance",
        acceptance="supplied synthetic conditional weights; not recovered experimental efficiency")
end

function bin_analytic_context(native,process,model,radial,Q2,P,x,z1,z2)
    plan=leptonic_plan(process,model,radial,Q2,P,z1;z2=process==:SIA ? z2 : nothing)
    scales=css_mu_b2.(Ref(model),radial.b)
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    first=HoppetMatchingPlan(native,process==:SIDIS ? [x] : [z1,z2],scales,model.Q02,model.coupling;
        process=process,quantity=kind,intervals=512)
    bases=process==:SIDIS ? [(i==1 ? 1.3 : .8)*(1-v) for i in 1:2,v in first.grid] :
        [PolarizedReference.FF_BASES[i]*(1-v)^3 for i in 1:4,v in first.grid]
    samples=PolarizedInitialSamples(kind,first.grid,bases;input_kind=:analytic)
    second=process==:SIDIS ? HoppetMatchingPlan(native,[z1],scales,model.Q02,model.coupling;process=process,intervals=512) : nothing
    samples2=second===nothing ? nothing : PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,second.grid,
        [PolarizedReference.FF_BASES[i]*(1-v)^3 for i in 1:4,v in second.grid];input_kind=:analytic)
    U=[1+.1pid for pid in plan.flavors,b in radial.b]
    D=[.8-.05pid for pid in plan.flavors,b in radial.b]
    polarized_leptonic_plan(plan,U,D,first,samples;second_operator=second,second_samples=samples2,
        target=process==:SIDIS ? NucleonTarget(:proton) : nothing,charges=process==:SIDIS ? (1,0) : (1,-1),
        provenance=Dict{String,Any}("validation_kinematics"=>Dict("x"=>x,"z1"=>z1,"z2"=>z2,"Q2"=>Q2,"P"=>P)))
end

function synthetic_acceptance_checks(contexts,theta,policy,checks)
    rows=[];settings=policy["synthetic_joint_acceptance"]
    uv,uw=CollinsFit.gauss_rule(settings["auxiliary_order"]);v=(uv .+1)./2;vw=uw./2
    context=first(contexts[:SIDIS]);base=prepare_generic_reference(context)
    R=base.polarized_integrals(theta)[1]/base.denominators[1]
    lower,upper=Float64.(settings["SIDIS_y_interval"])
    # Integral of (1+1.5y)*2(1-y)/y divided by
    # integral of (1+1.5y)*(1+(1-y)^2)/y, from the Eq.3 prefactor.
    numerator(y)=2log(y)+y-1.5y^2
    denominator(y)=2log(y)+y-y^2+.5y^3
    exact=R*(numerator(upper)-numerator(lower))/(denominator(upper)-denominator(lower))
    for order in settings["SIDIS_y_orders"]
        q,w=CollinsFit.gauss_rule(order);ys=(lower+upper)/2 .+(upper-lower)/2 .*q;yw=(upper-lower)/2 .*w
        nodes=LeptonicObservablePlan[];weights=Float64[]
        for (y,dy) in zip(ys,yw),(u,du) in zip(v,vw)
            push!(nodes,sidis_observable(context;definition=:sin_phi_sum_amplitude,y=y,source="synthetic y/u acceptance"))
            push!(weights,dy*du*(1+3y*u)*(1+(1-y)^2)/y)
        end
        value=predict_bin(validation_bin(nodes,weights,"synthetic dy du, not physical fixed-beam phase space"),theta)
        error=observable_scaled([value],[exact],policy)
        push!(rows,Dict("process"=>"SIDIS","y_order"=>order,"node_count"=>length(nodes),"value"=>value,"closed_reference"=>exact,"scaled_error"=>error))
    end
    observable_check!(checks,"synthetic correlated y/u quadrature final closed integral",last(rows)["scaled_error"],policy["analytic_integral_max_scaled"])
    context=first(contexts[:SIA]);reference=prepare_generic_reference(context)
    values=reference.polarized_integrals(theta);D=reference.denominators
    NU,NL=values[2]+values[3],values[1]+values[4];DU,DL=D[2]+D[3],D[1]+D[4]
    q,w=CollinsFit.gauss_rule(settings["SIA_cos_theta_order"])
    for ratio in (:UL,:UC)
        nodes=LeptonicObservablePlan[];weights=Vector{Float64}[]
        for (c,dc) in zip(q,w),(u,du) in zip(v,vw)
            a=(1-c^2)/(1+c^2)
            push!(nodes,sia_observable(context;ratio=ratio,angular=SIAAngularFactor(a;source="synthetic point theta"),source="synthetic theta/u acceptance"))
            push!(weights,dc*du*(1+c^2).*[1+3c^2*u,1+.5c^2*(1-u)])
        end
        W=reduce(vcat,permutedims.(weights))
        value=predict_bin(validation_bin(nodes,W,"synthetic dcos(theta) du"),theta)
        exact=NU*(26/15)/(DU*(64/15))-(ratio==:UL ? NL*(21/15)/(DL*(44/15)) :
            (NU*(26/15)+NL*(21/15))/(DU*(64/15)+DL*(44/15)))
        error=observable_scaled([value],[exact],policy)
        observable_check!(checks,"synthetic joint theta/u $ratio closed integral",error,policy["analytic_integral_max_scaled"])
        push!(rows,Dict("process"=>"SIA","ratio"=>string(ratio),"value"=>value,"closed_reference"=>exact,"scaled_error"=>error,
            "exact_U_sin_cos_integrals"=>[26/15,64/15],"exact_L_sin_cos_integrals"=>[21/15,44/15]))
    end
    rows
end

function validate_bins(mode,nativefile,inputfile,progressfile)
    mode in ("analytic","inputs") || error("mode must be analytic or inputs")
    sources=bin_sources();native=CollinsFit.verified_hoppet_build(nativefile)
    build_hashes=Dict(nativefile=>bytes2hex(sha256(read(nativefile))))
    contract=JSON.parsefile(BIN_CONTRACT);policy=contract["validation_policy"]
    checks=Dict{String,Any}[];cases=[];derivatives=[];contexts=Dict{Symbol,Any}()
    model=KPSYCSS();theta=collect(kpsy_published_parameters().values)
    # Coarse composition fixture, not a replacement for parent physical refinement.
    radial=RadialRule([0.,.25,.75,2.,8.];order=4)
    inputbuild=nothing;pdf=nothing;hel=nothing;ff=nothing
    if mode=="inputs"
        inputbuild=CollinsFit.verify_parton_build(inputfile);build_hashes[inputfile]=bytes2hex(sha256(read(inputfile)))
        pdf=CT10Input(inputfile;policy=:reference);hel=load_dssv08(inputbuild["input_root"]);ff=load_dss14(inputbuild["input_root"])
    end
    try
        for process in (:SIDIS,:SIA)
            settings=process==:SIDIS ? [(4.,.15,.1,.3,nothing),(6.,.3,.2,.4,nothing)] :
                [(110.,.15,0.,.3,.5),(110.,.3,0.,.4,.6)]
            contexts[process]=[mode=="analytic" ? bin_analytic_context(nativefile,process,model,radial,Q2,P,x,z1,z2) :
                process==:SIDIS ? prepare_polarized_sidis(pdf,hel,ff,nativefile,model,radial,Q2,P,x,z1,1;
                    matching_order=4,intervals=512,range_policy=:reference) :
                prepare_polarized_sia(ff,nativefile,model,radial,Q2,P,z1,z2,1,-1;
                    matching_order=4,intervals=512,range_policy=:reference) for (Q2,P,x,z1,z2) in settings]
            for integrated in (false,true)
                intervals=integrated ? [(.1,.3),(.2,.4)] : [nothing,nothing]
                momenta=[interval===nothing ? nothing : MomentumRule(interval...;order=128) for interval in intervals]
                definitions=process==:SIDIS ? (:structure_function_ratio,:sin_phi_sum_amplitude) : (:UL,:UC)
                for definition in definitions
                    ys=[.4,.3] # Q2/(x*S) for the explicit synthetic S=100 GeV2.
                    angles=[.2,.8]
                    nodes=process==:SIDIS ? [sidis_observable(c;definition=definition,y=ys[i],momentum=momenta[i],source="synthetic fixed-beam bin node") for (i,c) in enumerate(contexts[process])] :
                        [sia_observable(c;ratio=definition,angular=SIAAngularFactor(angles[i];source="synthetic conditional point-angle factor"),momentum=momenta[i],source="synthetic conditional bin node") for (i,c) in enumerate(contexts[process])]
                    weights=process==:SIDIS ? [w*(1+(1-y)^2)/(y*c.leptonic.Q2) for (w,y,c) in zip([.3,.7],ys,contexts[process])] :
                        [(.3+.2i+.1j)/(c.leptonic.Q2*(1+angles[i])) for (i,c) in enumerate(contexts[process]),j in 1:2]
                    factors=process==:SIDIS ? (definition==:structure_function_ratio ? ones(2) : [2(1-y)/(1+(1-y)^2) for y in ys]) : angles
                    bin=validation_bin(nodes,weights,"explicit fixed two-node quadrature in correlated kinematics; finite P included only when declared")
                    reference=prepare_bin_reference(nodes,weights,intervals,factors)
                    value=predict_bin(bin,theta);expected=reference(theta)
                    label="$process $definition integrated=$integrated"
                    error=observable_scaled([value],[expected],policy)
                    observable_check!(checks,label*" independent scalar/adaptive composition",error,policy["algebra_max_scaled"])
                    append!(derivatives,bin_derivative_checks(bin,reference,theta,policy,checks,label))
                    push!(cases,Dict("label"=>label,"parameters"=>theta,"value"=>value,"reference"=>expected,"scaled_error"=>error,
                        "weights"=>weights,"independent_angular_factors"=>factors,"momentum_intervals"=>intervals,
                        "integrated_denominators"=>bin.denominators,"node_provenance"=>bin.provenance["node_provenance"]))
                    writejson(progressfile,Dict("status"=>"RUNNING","scope"=>"reconstructed_theory","mode"=>mode,"completed_cases"=>cases,"checks"=>checks))
                    println(JSON.json(Dict("completed"=>label,"scaled_error"=>error)));flush(stdout)
                end
            end
        end
        joint=synthetic_acceptance_checks(contexts,theta,policy,checks)
        stable=bin_sources()==sources && all(bytes2hex(sha256(read(p)))==h for (p,h) in build_hashes) && CollinsFit.verified_hoppet_build(nativefile)==native
        mode=="inputs" && (stable=stable && CollinsFit.verify_parton_build(inputfile)==inputbuild)
        observable_check!(checks,"source and build stability",stable ? 0. : 1.,0.)
        Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
            "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,"historical_fit_ready"=>false,
            "contract_id"=>contract["id"],"mode"=>mode,"checks"=>checks,"sources"=>sources,"build_sha256"=>build_hashes,
            "parameters"=>theta,"cases"=>cases,"derivative_checks"=>derivatives,"synthetic_joint_acceptance"=>joint,
            "composition_grid"=>Dict("radial_edges"=>radial.edges,"radial_order"=>radial.order,"native_intervals"=>512,"matching_order"=>4,
                "scope"=>"coarse fixed-row composition only; parent independent evolution/refinement remains required"),
            "input_scope"=>(mode=="inputs" ? "strict pinned initial CT10/DSSV08/DSS14; explicit reference scale extrapolation; synthetic acceptance" : "analytic polynomial inputs; retained native rows; synthetic acceptance"),
            "limitations"=>contract["limitations"],"julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
    finally
        pdf!==nothing && close(pdf)
    end
end

length(ARGS)==4 || error("usage: check_leptonic_bins.jl analytic|inputs native-build.json input-build.json|- output.json")
mode,nativefile,inputfile,output=ARGS
try
    report=validate_bins(mode,nativefile,inputfile,output*".progress.json");writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"mode"=>mode,"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "sources"=>bin_sources(),"error"=>sprint(showerror,exception)))
    rethrow()
end
