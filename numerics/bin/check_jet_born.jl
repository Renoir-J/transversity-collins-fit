# Complete Julia reconstruction check. No data fit, evolution or production admission.
using CollinsFit, JSON, SHA, Dates, ForwardDiff, LinearAlgebra, QuadGK
using SpecialFunctions: besselj
include("JetBornReference.jl")
const JR=JetBornReference
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONTRACT=joinpath(ROOT,"config/star-born-contract.json")
const IDS=(21,2,-2,1,-1,3,-3)
const YS=(1//32,1//8,1//4,1//2,3//4,7//8,31//32)

function source_hashes()
    files=vcat([@__FILE__,CONTRACT,joinpath(@__DIR__,"JetBornReference.jl"),
        joinpath(ROOT,"scripts/prepare_star_born_references.py"),
        joinpath(ROOT,"config/initial-contract.json"),joinpath(ROOT,"numerics/Project.toml"),joinpath(ROOT,"numerics/Manifest.toml")],
        filter(p->endswith(p,".jl"),readdir(joinpath(ROOT,"numerics/src");join=true)))
    Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in files)
end
scaled(a,b,policy)=maximum(abs.(a.-b)./max.(abs.(b),policy["relative_floor"]))
function check!(checks,name,error,tolerance)
    push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error)&&error<=tolerance))
end

function source_algebra(references,policy,checks)
    rows=[]
    for y in YS,scale in (1,17)
        s=big(scale)//1;t=-s*y;u=-s*(1-y)
        env=Dict{Symbol,Any}(:s=>s,:t=>t,:u=>u,:D=>4//1,:CA=>3//1,:CF=>4//3,:Nc=>3//1,:v=>1-y,:alphaS=>1//5)
        for (name,expression) in sort(collect(references.FeynUU);by=first)
            ids=JR.tagged_ids(name);uk,_=CollinsFit.jet_born_channel(ids...)
            actual=CollinsFit.jet_born_hard(uk,s,t,u);expected=JR.arithmetic_expression(expression,env)
            check!(checks,"released UU $name y=$y s=$scale",actual==expected ? 0. : 1.,0.)
            push!(rows,Dict("family"=>"FeynFacet_UU","channel"=>name,"y"=>string(y),"s"=>scale,"actual"=>string(actual),"reference"=>string(expected)))
        end
        for (class,row) in sort(collect(references.FeynUT);by=first)
            ids=JR.tagged_ids(row.process);kinds=CollinsFit.jet_born_channel(ids...)
            actual=kinds===nothing ? zero(s) : CollinsFit.jet_born_transfer(kinds[2],s,t,u)
            for component in (:XX,:NN)
                expected=JR.arithmetic_expression(getproperty(row,component),env)*s*(1-y)/env[:alphaS]^2
                check!(checks,"released transfer class$class $component y=$y s=$scale",actual==expected ? 0. : 1.,0.)
                push!(rows,Dict("family"=>"FeynFacet_transfer","class"=>class,"channel"=>row.process,"component"=>string(component),
                    "y"=>string(y),"s"=>scale,"actual"=>string(actual),"reference"=>string(expected)))
            end
        end
    end
    five=(21,1,-1,2,-2,3,-3,4,-4,5,-5)
    routing=[]
    for y in YS
        s=big(17)//1;t=-s*y;u=-s*(1-y)
        table=JR.exported_channel_table(references,five,s,t,u)
        failures=[];maxU=0.;maxT=0.;nonzeroU=0;structuralT=0
        for a in five,b in five,c in five
            kinds=CollinsFit.jet_born_channel(a,b,c)
            U=kinds===nothing ? zero(s) : CollinsFit.jet_born_hard(kinds[1],s,t,u)
            T=kinds===nothing ? zero(s) : CollinsFit.jet_born_transfer(kinds[2],s,t,u)
            rU=get(table.UU,(a,b,c),zero(s));rT=get(table.UT,(a,b,c),zero(s))
            maxU=max(maxU,Float64(abs(U-rU)));maxT=max(maxT,Float64(abs(T-rT)))
            U==rU && T==rT || push!(failures,Dict("ids"=>[a,b,c],"actual"=>string.([U,T]),"reference"=>string.([rU,rT])))
            nonzeroU+=!iszero(U);structuralT+=kinds!==nothing && kinds[2]!=:zero
        end
        check!(checks,"all 1331 five-flavor routes against relabeled released channels y=$y",isempty(failures) ? 0. : 1.,0.)
        check!(checks,"five-flavor route counts y=$y",nonzeroU==331 && structuralT==120 ? 0. : 1.,0.)
        push!(routing,Dict("y"=>string(y),"s"=>17,"routes"=>11^3,"nonzero_UU"=>nonzeroU,"structural_UT"=>structuralT,
            "maximum_absolute_UU_error"=>maxU,"maximum_absolute_UT_error"=>maxT,"mismatches"=>failures))
    end
    (coefficients=rows,routing=routing)
end

function legacy_plan_checks(references,policy,checks)
    result=[]
    for args in ((200.,20.,0.4,0.4),(510.,31.,0.2,0.6)),fragmentation in (:quark_only,:quark_and_gluon)
        plan=JetBornPlan(IDS,JetBornKinematics(args...);fragmentation=fragmentation)
        k=plan.kinematics;actuals=[];expected=[]
        for a in 1:7,b in 1:7,c in 1:7
            f=zeros(7);ft=zeros(7);D=zeros(7);f[a]=ft[b]=D[c]=1.
            actual=jet_born_contractions(plan,f,ft,D,f,D)
            referenceD=copy(D);fragmentation==:quark_only && (referenceD[1]=0.)
            reference=JR.legacy_density(references,k.s,k.t,k.u,f,ft,referenceD,f,referenceD)
            push!(actuals,[actual.UU,actual.UT]);push!(expected,[Float64(reference.UU),Float64(reference.UT)])
        end
        errors=[scaled(a,b,policy) for (a,b) in zip(actuals,expected)]
        check!(checks,"legacy 343 one-hot flavors $(args) $fragmentation",maximum(errors),policy["algebra_max_scaled"])
        push!(result,Dict("kinematics"=>collect(args),"fragmentation"=>string(fragmentation),"flavors"=>collect(IDS),
            "iteration_order"=>"a outermost, b middle, c innermost; each1:7","actual_UU_UT"=>actuals,"reference_UU_UT"=>expected,"maximum_scaled_error"=>maximum(errors)))
    end
    result
end

# This family is explicitly synthetic. Fixed input bases are shared; source
# reference contractions and direct Gaussian C are assembled separately.
function response_arrays(theta,plans)
    p=KPSYParameters(theta);n=length(plans);z=(.23,.48);jt=(.19,.34);charges=(1,-1)
    f=[(1+.04i)*(1-plans[j].kinematics.x_a)^2 for i in 1:7,j in 1:n]
    ft=[(1+.03i)*(1-plans[j].kinematics.x_b)^3 for i in 1:7,j in 1:n]
    D=[(.8+.06i)*(1-z[j])^2*exp(-jt[j]^2/(4p.gh))/(pi*4p.gh) for i in 1:7,j in 1:n]
    h=[pid==21 ? zero(theta[1]) : initial_transversity(p,pid,plans[j].kinematics.x_a,(1+.1i)*(1-plans[j].kinematics.x_a)^2) for (i,pid) in enumerate(IDS),j in 1:n]
    H=[pid==21 ? zero(theta[1]) : initial_collins(p,pid,charges[j],z[j],[1.2,.7,.4,.3].*(1-z[j])^3) for pid in IDS,j in 1:n]
    w=4*(p.gh-p.values[13])
    C=[2z[j]*jet_collins_fourier_moment(H[i,j],z[j])*jt[j]*exp(-jt[j]^2/w)/(pi*w^2) for i in 1:7,j in 1:n]
    direct=[-H[i,j]*jt[j]*exp(-jt[j]^2/w)/(pi*w^2) for i in 1:7,j in 1:n]
    (f=f,ft=ft,D=D,h=h,C=C,direct=direct)
end

function synthetic_predictions(references,plans,theta;reference=false)
    weights=(.3,.7);couplings=(.2,.17);out=[]
    for group in plans
        arrays=response_arrays(theta,group)
        if reference
            U=zero(theta[1]);T=zero(theta[1])
            for j in eachindex(group)
                k=group[j].kinematics;D=copy(arrays.D[:,j]);C=copy(arrays.direct[:,j])
                group[j].fragmentation==:quark_only && (D[1]=0.;C[1]=zero(theta[1]))
                value=JR.legacy_density(references,k.s,k.t,k.u,arrays.f[:,j],arrays.ft[:,j],D,arrays.h[:,j],C)
                pref=couplings[j]^2/(k.S*k.x_a*k.x_b*(k.x_a*k.S+k.U))*weights[j]
                U+=pref*value.UU;T+=pref*value.UT
            end
            append!(out,[T/U,-T/U])
        else
            values=jet_born_integrals(group,weights,couplings,arrays.f,arrays.ft,arrays.D,arrays.h,arrays.C;
                measure_source="synthetic fixed two-node response; no experimental acceptance")
            append!(out,[jet_born_asymmetry(values;phase=:sin_phiS_minus_phiH),jet_born_asymmetry(values;phase=:sin_phiH_minus_phiS)])
        end
    end
    collect(promote(out...))
end

function central_gradient(f,t,h)
    columns=map(eachindex(t)) do i
        e=zeros(eltype(t),length(t));e[i]=h
        (-f(t+2e)+8f(t+e)-8f(t-e)+f(t-2e))/(12h)
    end
    reduce(hcat,columns)
end

function derivative_checks(references,policy,checks,progressfile)
    plans=[[JetBornPlan(IDS,JetBornKinematics(args...);fragmentation=model) for args in
        ((200.,20.,0.,.3),(200.,30.,0.,.5))] for model in (:quark_only,:quark_and_gluon)]
    theta=collect(kpsy_published_parameters().values);theta[2]=-.91;theta[11]=.4;theta[12]=.7
    actual(t)=synthetic_predictions(references,plans,t)
    reference(t)=synthetic_predictions(references,plans,t;reference=true)
    value=actual(theta);expected=reference(theta)
    check!(checks,"synthetic two-node response vs independent legacy/Gaussian assembly",scaled(value,expected,policy),policy["algebra_max_scaled"])
    J=ForwardDiff.jacobian(actual,theta)
    direction=collect(1.:13.);direction[13]=.01;direction./=norm(direction)
    H=reshape(ForwardDiff.jacobian(s->vec(ForwardDiff.jacobian(actual,theta+s[1]*direction)),[0.])[:,1],size(J))
    check!(checks,"synthetic all 13 parameter columns active and finite",all(isfinite,J)&&all(isfinite,H)&&all(sum(abs.(J);dims=1).>1e-12) ? 0. : 1.,0.)
    rows=[]
    for step in policy["finite_difference_steps"]
        g,hv=setprecision(policy["finite_difference_reference_bits"]) do
            t=BigFloat.(theta);h=BigFloat(step);d=BigFloat.(direction)
            gradient(x)=central_gradient(reference,x,h)
            (Float64.(gradient(t)),Float64.((-gradient(t+2h*d)+8gradient(t+h*d)-8gradient(t-h*d)+gradient(t-2h*d))/(12h)))
        end
        eg=scaled(J,g,policy);eh=scaled(H,hv,policy)
        check!(checks,"synthetic full-13 fourth-order FD gradient step$step",eg,policy["derivative_max_scaled"])
        check!(checks,"synthetic full-13 fourth-order FD HVP step$step",eh,policy["derivative_max_scaled"])
        push!(rows,Dict("step"=>step,"parameters"=>theta,"direction"=>direction,"actual_gradient"=>J,"reference_gradient"=>g,
            "actual_HVP"=>H,"reference_HVP"=>hv,"gradient_error"=>eg,"HVP_error"=>eh,"reference_bits"=>policy["finite_difference_reference_bits"]))
        writejson(progressfile,Dict("status"=>"RUNNING","scope"=>"reconstructed_theory","checks"=>checks,"completed_derivative_checks"=>rows))
        println(JSON.json(Dict("completed"=>"synthetic derivatives","step"=>step,"gradient_error"=>eg,"HVP_error"=>eh)));flush(stdout)
    end
    Dict("scope"=>"synthetic_response_controls_not_evolved_STAR_predictions","values"=>value,"reference_values"=>expected,
        "definitions"=>["quark_only phase S-H","quark_only phase H-S","quark_and_gluon phase S-H","quark_and_gluon phase H-S"],
        "checks"=>rows,"shared_inputs"=>"analytic Soffer/FF bases and existing initial-family functions; no independent evolution claimed")
end

# Independent recoil-rapidity measure; do not call production kinematics or prefactor.
function rapidity_reference(references,rs,pt,eta)
    ratio=pt/rs;S=rs^2
    lower=-log(rs/pt-exp(-eta));upper=log(rs/pt-exp(eta))
    function density(yd)
        xa=ratio*(exp(eta)+exp(yd));xb=ratio*(exp(-eta)+exp(-yd))
        s=xa*xb*S;t=-xa*rs*pt*exp(-eta);u=-xb*rs*pt*exp(eta)
        f=[(1+.04i)*(1-xa)^2 for i in 1:7];ft=[(1+.03i)*(1-xb)^3 for i in 1:7]
        D=[(.8+.06i) for i in 1:7];h=f.*[0.,.3,-.1,-.2,.07,.05,-.04];C=D.*[0.,-.2,.1,.12,-.05,.03,-.02]
        value=JR.legacy_density(references,s,t,u,f,ft,D,h,C)
        [value.UU,value.UT]./(xa*xb)
    end
    # Scalar component integrals also bound the small signed UT reference;
    # a vector norm dominated by UU would not establish that accuracy.
    parts=[quadgk(yd->density(yd)[i],lower,upper;rtol=1e-13,atol=1e-13) for i in 1:2]
    (values=first.(parts),errors=last.(parts),bounds=[lower,upper])
end

function integration_checks(references,policy,checks)
    rows=[]
    for (rs,pt,eta) in ((200.,20.,.4),(510.,31.,.2))
        reference=rapidity_reference(references,rs,pt,eta);xmin=jet_born_xmin(rs,pt,eta);alpha=.2
        reference_error=maximum(reference.errors./max.(abs.(reference.values),policy["relative_floor"]))
        check!(checks,"independent recoil-rapidity component error estimates rs=$rs",reference_error,policy["algebra_max_scaled"])
        refinements=[]
        for order in policy["integration_orders"]
            q,w=CollinsFit.gauss_rule(order);xs=(1+xmin)/2 .+(1-xmin)/2 .*q;weights=(1-xmin)/2 .*w
            plans=[JetBornPlan(IDS,JetBornKinematics(rs,pt,eta,x);fragmentation=:quark_and_gluon) for x in xs]
            f=[(1+.04i)*(1-k.kinematics.x_a)^2 for i in 1:7,k in plans]
            ft=[(1+.03i)*(1-k.kinematics.x_b)^3 for i in 1:7,k in plans]
            D=[.8+.06i for i in 1:7,j in eachindex(plans)]
            h=f.*[0.,.3,-.1,-.2,.07,.05,-.04];C=D.*[0.,-.2,.1,.12,-.05,.03,-.02]
            value=jet_born_integrals(plans,weights,fill(alpha,length(plans)),f,ft,D,h,C;
                measure_source="synthetic dx_a integration, compared with independent recoil rapidity")
            normalized=[value.UU,value.UT].*(rs^4/alpha^2)
            error=scaled(normalized,reference.values,policy)
            push!(refinements,Dict("order"=>order,"normalized_UU_UT"=>normalized,"reference_normalized_UU_UT"=>reference.values,
                "scaled_error"=>error,"meets_final_tolerance"=>error<=policy["integration_max_scaled"]))
        end
        check!(checks,"delta-Jacobian x_a vs recoil-rapidity final integral rs=$rs",last(refinements)["scaled_error"],policy["integration_max_scaled"])
        push!(rows,Dict("sqrt_s"=>rs,"pT"=>pt,"eta"=>eta,"x_a_min"=>xmin,"normalization"=>"S^2/alpha_s^2 times both yields before relative-floor comparison",
            "reference_adaptive_component_errors"=>reference.errors,"reference_max_scaled_error_estimate"=>reference_error,
            "recoil_rapidity_bounds"=>reference.bounds,"refinements"=>refinements))
    end
    rows
end

function gaussian_checks(policy,checks)
    rows=[]
    for z in (.23,.48,.7),width in (.08,.16,.3),j in (0.,.19,.7)
        H=-.07;moment=jet_collins_fourier_moment(H,z)
        # STAR Eq.6b with a Gaussian in b conjugate to j/z. z factors
        # are part of this check, rather than absorbed into an arbitrary norm.
        transform=quadgk(b->moment*b^2*besselj(1,j*b/z)*exp(-width*b^2/(4z^2))/(2pi*z^2),0.,Inf;rtol=1e-12,atol=1e-14)[1]
        direct=-H*j*exp(-j^2/width)/(pi*width^2)
        error=scaled([transform],[direct],policy)
        check!(checks,"synthetic Gaussian J1 mass/z normalization z=$z w=$width j=$j",error,policy["integration_max_scaled"])
        push!(rows,Dict("z"=>z,"width"=>width,"jperp"=>j,"Hhat3"=>H,"fourier_J1"=>transform,"direct_C"=>direct,"scaled_error"=>error))
    end
    rows
end

function validate(referencefile,output)
    sources=source_hashes();contract=JSON.parsefile(CONTRACT);policy=contract["validation_policy"]
    reference_sha=bytes2hex(sha256(read(referencefile)));references=JR.load_source_references(referencefile)
    checks=Dict{String,Any}[]
    algebra=source_algebra(references,policy,checks);legacy=legacy_plan_checks(references,policy,checks)
    writejson(output*".progress.json",Dict("status"=>"RUNNING","scope"=>"reconstructed_theory","production_ready"=>false,
        "checks"=>checks,"released_coefficient_checks"=>algebra.coefficients,"five_flavor_routing"=>algebra.routing,"legacy_one_hot_checks"=>legacy))
    println(JSON.json(Dict("completed"=>"source algebra and one-hot plans","checks"=>length(checks))));flush(stdout)
    derivative=derivative_checks(references,policy,checks,output*".progress.json")
    integrals=integration_checks(references,policy,checks);gaussians=gaussian_checks(policy,checks)
    check!(checks,"source and reference stability",source_hashes()==sources && bytes2hex(sha256(read(referencefile)))==reference_sha && JR.load_source_references(referencefile).data==references.data ? 0. : 1.,0.)
    Dict("status"=>all(c["passed"] for c in checks) ? "PASS" : "FAIL","scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,"full_STAR_TMD_observable_implemented"=>false,
        "native_FeynFacet_calculation_executed"=>false,"legacy_imports_or_pickle_deserialization_executed"=>false,
        "contract_id"=>contract["id"],"sources"=>sources,"reference_data"=>Dict("path"=>referencefile,"sha256"=>reference_sha,"sources"=>references.data["sources"]),
        "checks"=>checks,"released_coefficient_checks"=>algebra.coefficients,"five_flavor_routing"=>algebra.routing,
        "legacy_one_hot_checks"=>legacy,"synthetic_response_controls"=>derivative,"synthetic_delta_measure_controls"=>integrals,
        "synthetic_gaussian_controls"=>gaussians,"limitations"=>contract["not_implemented"],"julia_version"=>string(VERSION),"completed_at"=>string(now(UTC)))
end

length(ARGS)==2 || error("usage: check_jet_born.jl reference-data.json output.json")
referencefile,output=ARGS
try
    report=validate(referencefile,output);writejson(output,report)
    println(JSON.json(Dict("status"=>report["status"],"checks"=>length(report["checks"]),"failed"=>filter(c->!c["passed"],report["checks"]))));flush(stdout)
    report["status"]=="PASS" || exit(2)
catch exception
    writejson(output,Dict("status"=>"FAIL","scope"=>"reconstructed_theory","production_ready"=>false,"physical_extraction"=>false,
        "sources"=>source_hashes(),"error"=>sprint(showerror,exception),
        "partial_progress"=>(isfile(output*".progress.json") ? Dict("path"=>output*".progress.json","sha256"=>bytes2hex(sha256(read(output*".progress.json")))) : nothing)))
    rethrow()
end
