# One complete Julia validation job. No experiment or production job is admitted.
using CollinsFit, JSON, SHA, Dates, LinearAlgebra, ForwardDiff, QuadGK
using SpecialFunctions: besselj
include("StarEvolvedReference.jl")
include("StarEvolvedFixture.jl")
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/star-evolved-contract.json")
const SOURCEFILES=vcat([@__FILE__,joinpath(@__DIR__,"StarEvolvedReference.jl"),joinpath(@__DIR__,"StarEvolvedFixture.jl"),
    joinpath(@__DIR__,"LeptonicReference.jl"),joinpath(@__DIR__,"MatchedReference.jl"),joinpath(@__DIR__,"CollinearReference.jl"),
    joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml"),joinpath(ROOT,"fortran/hoppet_bridge.f90")],
    filter(p->endswith(p,"-contract.json") || endswith(p,"parton-sources.json"),readdir(joinpath(ROOT,"config");join=true)),
    filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
source_hashes()=Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in SOURCEFILES)

function central_gradient(f,t,h)
    [begin
        e=zeros(eltype(t),length(t));e[i]=h
        (-f(t+2e)+8f(t+e)-8f(t-e)+f(t-2e))/(12h)
    end for i in eachindex(t)]
end
function full_derivatives(p,theta,policy,check,progress)
    reference=StarEvolvedReference.scalar_reference(p)
    f(t)=predict_star(p,t;phase=:sin_phiS_minus_phiH)
    r(t)=reference(t).UT/reference(t).UU
    J=ForwardDiff.gradient(f,theta);d=collect(1.:13.);d[13]=.01;d./=norm(d)
    H=ForwardDiff.jacobian(s->ForwardDiff.gradient(f,theta+s[1]*d),[0.])[:,1]
    check("all13 active finite gradient and directional HVP",all(isfinite,J)&&all(isfinite,H)&&all(abs.(J).>1e-12) ? 0. : 1.,0.)
    rows=[]
    scaled(a,b)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
    for step in policy["finite_difference_steps"]
        gradient,hvp=setprecision(policy["finite_difference_reference_bits"]) do
            t=BigFloat.(theta);h=BigFloat(step);dd=BigFloat.(d)
            G(x)=central_gradient(r,x,h)
            (Float64.(G(t)),Float64.((-G(t+2h*dd)+8G(t+h*dd)-8G(t-h*dd)+G(t-2h*dd))/(12h)))
        end
        eg,eh=scaled(J,gradient),scaled(H,hvp)
        check("full13 scalar256bit gradient step$step",eg,policy["derivative_max_scaled"])
        check("full13 scalar256bit HVP step$step",eh,policy["derivative_max_scaled"])
        push!(rows,Dict("step"=>step,"reference_bits"=>policy["finite_difference_reference_bits"],
            "gradient"=>J,"reference_gradient"=>gradient,"HVP"=>H,"reference_HVP"=>hvp,
            "gradient_max_scaled"=>eg,"HVP_max_scaled"=>eh))
        progress("full13 derivatives step$step")
    end
    Dict("scope"=>"scalar_assembly_controls_shared_native_rows_and_provider_samples",
        "parameters"=>theta,"direction"=>d,"rows"=>rows,"point_asymmetry"=>f(theta),"reference_asymmetry"=>r(theta))
end

function validate(kind,nativefile,inputfile,output)
    kind in ("analytic","inputs") || error("explicit analytic or inputs validation job required")
    sources=source_hashes();contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    native=CollinsFit.verified_hoppet_build(nativefile)
    builds=Dict(nativefile=>bytes2hex(sha256(read(nativefile))))
    checks=Dict{String,Any}[];details=Dict{String,Any}()
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,
        "passed"=>error isa Real && isfinite(error)&&error<=tolerance))
    function progress(stage)
        writejson(output*".progress.json",Dict("status"=>"RUNNING","stage"=>stage,"scope"=>"reconstructed_theory",
            "production_ready"=>false,"sources"=>sources,"checks"=>checks,"details"=>details))
        println(JSON.json(Dict("completed"=>stage,"checks"=>length(checks))));flush(stdout)
    end
    scaled(a,b)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
    model=KPSYCSS();rs,pt,eta,z,jt=200.,12.9,.3,.3,.3
    polynomial=[.5,-.4,1.,1.,2.,2.,.3,-.2,1.,1.,0.,0.,.01]
    theta=collect(kpsy_published_parameters().values);theta[2]=-.91;theta[11]=.4;theta[12]=.7
    pdf=hel=ff=inputs=nothing
    if kind=="inputs"
        inputs=CollinsFit.verify_parton_build(inputfile);builds[inputfile]=bytes2hex(sha256(read(inputfile)))
        pdf=CT10Input(inputfile;policy=:reference);hel=load_dssv08(inputs["input_root"]);ff=load_dss14(inputs["input_root"])
    end
    normalized(p,t)=collect(star_structure(p,t)).*(rs^4/alpha_s(model.coupling,pt^2)^2)
    contexts=Dict{Symbol,StarPointPlan}();alllevels=Dict{String,Any}()
    try
        if kind=="analytic"
            # Standalone exact-model controls of the complete C11 Fourier
            # prefactor. They are explicitly not evolved physical predictions.
            gaussians=[];r=RadialRule([0.,.1,.4,1.,2.,4.,8.,16.,32.];order=32)
            for zz in (.23,.48,.7),width in (.08,.16,.3),j in (0.,.19,.7)
                H=-.07;A=width/(4zz^2)
                actual=dot(-CollinsFit.hankel_weights(r,1,j/zz)./(2zz^3),H.*exp.(-A.*r.b.^2))
                expected=-H*j*exp(-j^2/width)/(pi*width^2)
                check("synthetic Gaussian C11 J1 z$zz width$width j$j",scaled([actual],[expected]),policy["algebra_max_scaled"])
                push!(gaussians,Dict("z"=>zz,"width"=>width,"jT"=>j,"actual"=>actual,"reference"=>expected))
            end
            details["synthetic_gaussian_normalization_controls"]=gaussians
            r=RadialRule([0.,.1,1.,4.];order=4);qs=css_mu_b2.(Ref(model),r.b)
            star=HoppetMatchingPlan(nativefile,[z],qs,model.Q02,model.coupling;process=:STAR_C11,intervals=256)
            sia=HoppetMatchingPlan(nativefile,[z],qs,model.Q02,model.coupling;process=:SIA,intervals=256)
            plain=HoppetPointPlan(nativefile,[z],qs,model.Q02,model.coupling;quantity=:KPSY_Hhat3_homogeneous,intervals=256)
            source=star.grid.*(1 .-star.grid).^3
            actual=apply_operator(star,source)
            expected=apply_operator(sia,source).-permutedims(sia.alpha./pi.*((pi^2-8)/3)).*apply_operator(plain,source)
            check("native matched C11 equals SIA with its hard delta removed",scaled(actual,expected),policy["algebra_max_scaled"])
            progress("analytic Gaussian and native C11 delta controls")
        end
        for fragmentation in (:quark_only,:quark_and_gluon)
            label=string(fragmentation);levels=[];lastp=nothing
            t=kind=="analytic" ? polynomial : theta
            for (order,intervals) in zip(policy["radial_orders"],policy["evolution_intervals"])
                edges=kind=="analytic" ? star_edges(24.) : parton_radial_edges(model,ff,inputs["input_root"],24.)
                radial=RadialRule(edges;order=order)
                cold=@elapsed lastp=kind=="analytic" ?
                    analytic_star_context(nativefile,model,radial;fragmentation=fragmentation,outer_order=64,intervals=intervals) :
                    prepare_star_point(pdf,hel,ff,nativefile,model,radial,rs,pt,eta,z,jt,1;
                        fragmentation=fragmentation,gluon_model=star_gluon_model(fragmentation),jet_radius=.6,
                        axis=:standard,scale_prescription=:historical_pT_C11,outer_order=64,matching_order=8,
                        intervals=intervals,range_policy=:reference)
                values=normalized(lastp,t);previous=isempty(levels) ? nothing : last(levels)["normalized_UU_UT"]
                err=previous===nothing ? nothing : scaled(values,previous)
                push!(levels,Dict("radial_order"=>order,"evolution_intervals"=>intervals,"outer_order"=>64,
                    "radial_nodes"=>length(radial.b),"normalized_UU_UT"=>values,"previous_level_max_scaled"=>err,
                    "meets_final_tolerance"=>err!==nothing && err<=policy["physical_value_max_scaled"],
                    "preparation_seconds_informational"=>cold,"provenance"=>lastp.provenance))
                alllevels[label]=levels;details["coupled_refinements"]=alllevels
                progress("$label radial$order evolution$intervals")
            end
            check("$label final combined radial/evolution refinement",last(levels)["previous_level_max_scaled"],policy["physical_value_max_scaled"])
            ffplan=lastp.fragmentation;outerlevels=[]
            for order in policy["outer_orders"]
                p=order==64 ? lastp : star_on_outer_grid(nativefile,ffplan,rs,pt,eta;
                    outer_order=order,intervals=last(policy["evolution_intervals"]),pdf=pdf,hel=hel)
                value=normalized(p,t);previous=isempty(outerlevels) ? nothing : last(outerlevels)["normalized_UU_UT"]
                err=previous===nothing ? nothing : scaled(value,previous)
                push!(outerlevels,Dict("order"=>order,"normalized_UU_UT"=>value,"previous_level_max_scaled"=>err,
                    "meets_final_tolerance"=>err!==nothing && err<=policy["physical_value_max_scaled"]))
                contexts[fragmentation]=p;details["$(label)_outer_refinements"]=outerlevels
                progress("$label outer$order")
            end
            check("$label declared final outer refinement",last(outerlevels)["previous_level_max_scaled"],policy["physical_value_max_scaled"])
            p=contexts[fragmentation];scalar=StarEvolvedReference.scalar_reference(p)
            factor=rs^4/alpha_s(model.coupling,pt^2)^2
            check("$label direct scalar full assembly",scaled(factor.*collect(star_structure(p,t)),factor.*collect(scalar(t))),policy["algebra_max_scaled"])
            check("$label canonical/opposite phase",predict_star(p,t;phase=:sin_phiS_minus_phiH)==-predict_star(p,t;phase=:sin_phiH_minus_phiS) ? 0. : 1.,0.)
            check("$label jT0 differential numerator",abs(star_structure(with_momentum(p,0.),t).UT),0.)
            # The finer B=32 leg retains exactly the final native/outer policies.
            radial=RadialRule(kind=="analytic" ? star_edges(32.) : parton_radial_edges(model,ff,inputs["input_root"],32.);order=last(policy["radial_orders"]))
            extension=kind=="analytic" ? analytic_star_fragmentation(nativefile,model,radial,pt,z,jt,1,fragmentation;
                intervals=last(policy["evolution_intervals"])) : begin
                    scales=css_mu_b2.(Ref(model),radial.b)
                    H=CollinsFit.star_native_checked(()->HoppetMatchingPlan(nativefile,[z],scales,model.Q02,model.coupling;
                        process=:STAR_C11,intervals=last(policy["evolution_intervals"])))
                    zz,bases=CollinsFit.prepare_initial_ff(ff,H.grid)
                    sample=PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,zz,bases;input_kind=:pinned_inputs)
                    D=CollinsFit.matched_star_ff_inputs(ff,model,radial,z,1,p.fragmentation.pids,8,:reference;
                        fragmentation=fragmentation,gluon_model=star_gluon_model(fragmentation))
                    audit=Dict("beam_nodes"=>2length(p.born),"radial_nodes"=>length(scales),
                        "CT10_beam_nodes_outside_strict"=>count(x->!parton_domain(pdf,x,pt^2)["within_grid"],
                            vcat([b.kinematics.x_a for b in p.born],[b.kinematics.x_b for b in p.born])),
                        "DSS14_radial_nodes_outside_strict"=>count(q->!parton_domain(ff,z,q;policy=:reference)["within_grid"],scales))
                    meta=CollinsFit.prepared_input_metadata(:STAR_C11,model,radial,8,:reference,audit)
                    meta["outer_order"]=length(p.born);meta["input_kind"]="pinned_CT10_DSSV08_DSS14"
                    star_fragmentation_plan(model,radial,pt^2,z,jt,p.fragmentation.pids,D,H,sample;
                        fragmentation=fragmentation,gluon_model=star_gluon_model(fragmentation),charge=1,provenance=meta)
                end
            p32=star_point_plan(p.born,p.weights,p.couplings,p.fA,p.fB,p.operator,p.samples,extension;
                jet_radius=p.jet_radius,axis=p.axis,scale_prescription=p.scale_prescription,measure_source=p.measure_source,
                provenance=merge(deepcopy(p.provenance),extension.provenance))
            check("$label finite radial bound24to32",scaled(normalized(p32,t),normalized(p,t)),policy["physical_value_max_scaled"])
            details["$(label)_bound_extension"]=Dict("bounds"=>[24.,32.],"normalized_UU_UT"=>[normalized(p,t),normalized(p32,t)])
            if kind=="analytic"
                hard(a,b,c,s,tt,u)=begin
                    kinds=CollinsFit.jet_born_channel(a,b,c)
                    kinds===nothing ? (0.,0.) : (CollinsFit.jet_born_hard(kinds[1],s,tt,u),CollinsFit.jet_born_transfer(kinds[2],s,tt,u))
                end
                oracle=StarEvolvedReference.analytic_reference(model,rs,pt,eta,z,jt,1,fragmentation,polynomial,hard)
                factor=rs^4/alpha_s(model.coupling,pt^2)^2
                expected=collect(oracle.values).*factor
                check("$label independent Mellin/Fourier/recoil-rapidity prediction",scaled(normalized(p,polynomial),expected),policy["physical_value_max_scaled"])
                check("$label oracle radial estimated relative error",oracle.radial_estimated_error/max(norm(oracle.radial_integrals),policy["relative_floor"]),policy["physical_value_max_scaled"])
                check("$label oracle recoil estimated relative error",oracle.recoil_estimated_error/max(norm(oracle.recoil_integrals),policy["relative_floor"]),policy["physical_value_max_scaled"])
                details["$(label)_independent_oracle"]=Dict("normalized_UU_UT"=>expected,"oracle"=>oracle)
            end
            progress("$label complete point and reference checks")
        end
        # Different charge/z/jet cases exercise preparation beyond one baseline;
        # full13 derivatives use smaller complete contexts, with explicit scope.
        small=RadialRule(kind=="analytic" ? star_edges(24.) : parton_radial_edges(model,ff,inputs["input_root"],24.);order=8)
        derivatives=Dict{String,Any}();cases=[]
        for (rs2,pt2,eta2,z2,j2,c2) in ((200.,12.9,.3,.3,.3,1),(510.,31.,-.2,.5,.5,-1))
            p=kind=="analytic" ? analytic_star_context(nativefile,model,small;pt=pt2,sqrt_s=rs2,eta=eta2,z=z2,jT=j2,charge=c2,
                fragmentation=:quark_and_gluon,outer_order=16,intervals=256) :
                prepare_star_point(pdf,hel,ff,nativefile,model,small,rs2,pt2,eta2,z2,j2,c2;
                    fragmentation=:quark_and_gluon,gluon_model=star_gluon_model(:quark_and_gluon),jet_radius=.6,axis=:standard,
                    scale_prescription=:historical_pT_C11,outer_order=16,matching_order=8,intervals=256,range_policy=:reference)
            key="$(Int(rs2))GeV"
            derivatives[key]=full_derivatives(p,theta,policy,check,progress)
            push!(cases,p)
        end
        details["full13_derivatives"]=derivatives
        # One beam energy per bin, as required by the physical composition API.
        base=cases[1];high=with_momentum(base,.5);measure=[.3,1.7]
        bin=star_bin([base,high],measure;phase=:sin_phiS_minus_phiH,measure_source="synthetic two-node jT response, not measured acceptance")
        values=[star_structure(x,theta) for x in (base,high)]
        target=sum(measure[i]*values[i].UT for i in 1:2)/sum(measure[i]*values[i].UU for i in 1:2)
        check("supplied-bin absolute-yield ratio",scaled([predict_star(bin,theta)],[target]),policy["algebra_max_scaled"])
        averaged=sum(measure[i]*predict_star(x,theta;phase=:sin_phiS_minus_phiH) for (i,x) in enumerate((base,high)))/sum(measure)
        check("averaged point ratios negative control",abs(averaged-target)>1e-8 ? 0. : 1.,0.)
        mr=MomentumRule(.05,.6;order=16)
        mb=star_momentum_bin(base,mr;phase=:sin_phiS_minus_phiH,measure=:jT_djT,measure_source="synthetic finite momentum interval")
        ab=star_momentum_bin(base,mr;phase=:sin_phiS_minus_phiH,measure=:d2jT,measure_source="synthetic finite momentum interval with azimuth")
        check("explicit d2jT normalization",scaled(collect(star_structure(ab,theta)),2pi.*collect(star_structure(mb,theta))),policy["algebra_max_scaled"])
        check("finite momentum integral not jT0 sentinel",abs(star_structure(mb,theta).UT)>0 ? 0. : 1.,0.)
        integrated_reference=StarEvolvedReference.scalar_reference(base;momentum_interval=(.05,.6))(theta)
        factor=rs^4/alpha_s(model.coupling,pt^2)^2
        check("finite jT integral versus independent adaptive kernels",scaled(factor.*collect(star_structure(mb,theta)),factor.*collect(integrated_reference)),policy["physical_value_max_scaled"])
        details["supplied_bin_controls"]=Dict("scope"=>"synthetic_measure_not_experiment_acceptance","bin_asymmetry"=>target,
            "incorrect_average_of_ratios"=>averaged,"momentum_integrals"=>star_structure(mb,theta),"momentum_bounds"=>[.05,.6],"momentum_order"=>16)
        if kind=="inputs"
            p=contexts[:quark_and_gluon];f=p.fragmentation
            check("pinned reference extrapolation provenance",p.provenance["range_policy"]=="reference" &&
                p.provenance["initial_source_kind"]=="pinned_inputs" && !p.provenance["production_ready"] &&
                p.provenance["source_range_audits"]["DSS14_radial_nodes_outside_strict"]>0 ? 0. : 1.,0.)
            details["input_provenance"]=p.provenance
            # Source-input checks use direct interpolation, not the cached FF matrices.
            matches=[]
            for bindex in unique(round.(Int,range(1,length(f.radial.b);length=7)))
                q2=css_mu_b2(model,f.radial.b[bindex]);a=alpha_s(model.coupling,q2)
                q(y)=pion_ff(ff,2,1,y,q2;policy=:reference);g(y)=pion_ff(ff,21,1,y,q2;policy=:reference)
                singlet(y)=sum(pion_ff(ff,pid,1,y,q2;policy=:reference) for pid in f.pids if pid!=21)
                # Split in v at every transformed input knot. This is independent
                # adaptive integration, not matching_rule's log(y) Gauss weights.
                edges=sort!(unique(vcat(z,filter(v->z<v<1,z./ff.x),1.)))
                cf(v)=(4/3)*(((1-v)/2+(1+v^2)*log(v)/(1-v))*q(z/v)+(v/2+(2/v-2+v)*log(v))*g(z/v))/v
                gf(v)=(6*(v/(1-v)+(1-v)/v+v*(1-v))*log(v)*g(z/v)+.5*(v*(1-v)+(v^2+(1-v)^2)*log(v))*singlet(z/v))/v
                expected=[q(z)+a/pi*quadgk(cf,edges...;rtol=1e-11,atol=1e-13)[1],
                    g(z)+a/pi*quadgk(gf,edges...;rtol=1e-11,atol=1e-13)[1]]
                actual=f.matched[[findfirst(==(2),f.pids),findfirst(==(21),f.pids)],bindex]
                err=scaled(actual,expected)
                check("source-input quark/gluon C11 matching bindex$bindex",err,policy["physical_value_max_scaled"])
                push!(matches,Dict("bindex"=>bindex,"Q2"=>q2,"actual"=>actual,"reference"=>expected,"error"=>err))
            end
            details["independent_input_matching"]=matches
        end
        check("source/build stability",sources==source_hashes() && all(bytes2hex(sha256(read(p)))==h for (p,h) in builds) ? 0. : 1.,0.)
        CollinsFit.verified_hoppet_build(nativefile)==native || error("native build changed during STAR validation")
        kind=="inputs" && CollinsFit.verify_parton_build(inputfile)!=inputs && error("parton build changed during STAR validation")
        progress("all declared checks completed")
        Dict("status"=>all(x["passed"] for x in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
            "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,
            "coupled_STAR_TMD_reconstruction_implemented"=>true,"historical_author_code_recovered"=>false,
            "finite_R_matching_implemented"=>false,"native_FeynFacet_calculation_executed"=>false,
            "contract_id"=>contract["id"],"kind"=>kind,"sources"=>sources,"builds"=>builds,
            "checks"=>checks,"details"=>details,"limitations"=>contract["remaining_production_requirements"],
            "julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
    catch exception
        writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,
            "kind"=>kind,"sources"=>sources,"checks"=>checks,"details"=>details,"error"=>sprint(showerror,exception)))
        rethrow()
    end
end

length(ARGS)==4 || error("usage: check_star_evolved.jl analytic|inputs HOPPET_BUILD INPUT_BUILD OUTPUT")
mkpath(dirname(ARGS[4]))
diagnostic_path=ARGS[4]*".native-stderr.log"
report=open(diagnostic_path,"w") do stream
    redirect_stderr(stream) do
        validate(ARGS...)
    end
end
diagnostic=read(diagnostic_path,String)
precision_errors=count("TOO HIGH ACCURACY REQUIRED",diagnostic)
report["native_diagnostics"]=Dict("path"=>diagnostic_path,"sha256"=>bytes2hex(sha256(read(diagnostic_path))),
    "precision_error_count"=>precision_errors,"guard"=>"native precision failures reject even if ABI status is zero")
push!(report["checks"],Dict("name"=>"no native precision-failure diagnostics","error"=>precision_errors,"tolerance"=>0,"passed"=>precision_errors==0))
precision_errors>0 && (report["status"]="FAIL")
writejson(ARGS[4],report)
println(JSON.json(Dict("status"=>report["status"],"kind"=>report["kind"],"checks"=>length(report["checks"]))))
report["status"]=="PASS" || error("STAR reconstruction checks failed; retain the failed report")
