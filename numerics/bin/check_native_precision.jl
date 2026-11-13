# Fixed-grid regression and independent Mellin checks; no production admission.
using CollinsFit,JSON,SHA,Dates,LinearAlgebra
const PACKAGE_ROOT=normpath(joinpath(dirname(pathof(CollinsFit)),"../.."))
const DRIVER_ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(PACKAGE_ROOT,"numerics/bin/CollinearReference.jl"))
using .CollinearReference

function precision_sources()
    files=vcat([joinpath(PACKAGE_ROOT,p) for p in ("fortran/hoppet_bridge.f90","numerics/Project.toml","numerics/Manifest.toml",
        "config/hoppet-contract.json","numerics/bin/CollinearReference.jl")],
        filter(p->endswith(p,".jl"),readdir(joinpath(PACKAGE_ROOT,"numerics/src");join=true)))
    result=Dict("package/"*relpath(p,PACKAGE_ROOT)=>bytes2hex(sha256(read(p))) for p in files)
    result["driver/check_native_precision.jl"]=bytes2hex(sha256(read(@__FILE__)))
    result["driver/hoppet-contract.json"]=bytes2hex(sha256(read(joinpath(DRIVER_ROOT,"config/hoppet-contract.json"))))
    result
end

function capture_precision(f,path)
    value=open(path,"w") do stream
        redirect_stderr(stream) do
            try f() finally flush(stderr);Base.Libc.flush_cstdio() end
        end
    end
    text=read(path,String)
    value,Dict("path"=>path,"sha256"=>bytes2hex(sha256(read(path))),
        "precision_error_count"=>count("TOO HIGH ACCURACY REQUIRED",text))
end

function native_case(build,setting,construction)
    x=Float64(setting["xmin"]);dy=haskey(setting,"dy") ? Float64(setting["dy"]) : -log(x)/setting["intervals"]
    native=HoppetPlan(build,[1.2,2.4,10.,1000.],2.4,AsymptoticTwoLoopCoupling();
        xmin=x,dy=dy,layout=Symbol(setting["layout"]),construction=construction)
    try
        grid=native.grid;h=hcat(grid.*(1 .-grid).^3,-2grid.^2 .* (1 .-grid).^2,grid.^3 .* (1 .-grid).^3)
        values=apply_operator(native,h);rows=CollinsFit.hoppet_rows(native,1)
        generator=CollinsFit.hoppet_generator(native,h)
        queries=[x,(x+1)/2,1.];observed=hoppet_evaluate(native,values,queries)
        cotangent=reshape(cos.(1:length(values)),size(values))
        adjoint=operator_adjoint(native,cotangent)
        reference=[q==1 ? 0. : mellin_polynomial_tau(q,t) for q in queries,t in native.tau]
        floor=1e-8
        mellin_error=maximum(abs.(observed[:,1,:].-reference)./max.(abs.(reference),floor))
        forward=permutedims(values[1,:,:]);row_error=maximum(abs.(rows*h.-forward)./max.(abs.(forward),.01maximum(abs,forward),1e-15))
        payload=Dict("grid"=>copy(grid),"tau"=>copy(native.tau),"forward"=>vec(values),"rows"=>vec(rows),
            "generator"=>vec(generator),"interpolated"=>vec(observed),"adjoint"=>vec(adjoint))
        regular_error=nothing
        if native.layout==:single
            regular_rows=CollinsFit.hoppet_rows(native,1;collins=true)
            regular=CollinsFit.hoppet_collins_regular(native,values)
            truth=permutedims(regular[1,:,:])
            regular_error=maximum(abs.(regular_rows*h.-truth)./max.(abs.(truth),.01maximum(abs,truth),1e-15))
            payload["Collins_rows"]=vec(regular_rows);payload["Collins_forward"]=vec(regular)
        end
        Dict("label"=>setting["label"],"construction"=>string(construction),"xmin"=>x,"dy"=>dy,"layout"=>setting["layout"],
            "Q2_GeV2"=>copy(native.Q2),"Q02_GeV2"=>native.Q02,"nf"=>native.nf,"quantity"=>string(native.quantity),
            "payload"=>payload,"mellin_max_scaled"=>mellin_error,"row_native_max_scaled"=>row_error,
            "Collins_row_native_max_scaled"=>regular_error,"queries"=>queries,"reference"=>vec(reference))
    finally
        close(native)
    end
end

function check_precision(role,buildfile,output,parentfile=nothing)
    role in ("parent","candidate") || error("explicit parent or candidate role required")
    policy=JSON.parsefile(joinpath(DRIVER_ROOT,"config/hoppet-contract.json"))["native_precision"]["validation_policy"]
    policy["native_precision_errors_allowed"]==0 && policy["parent_retained_prediction_max_abs"]==0 || error("precision policy changed")
    parent=role=="candidate" ? JSON.parsefile(parentfile) : nothing
    if parent!==nothing
        parent["role"]=="parent" && parent["status"]=="REFERENCE_COMPLETED_WITH_DIAGNOSTICS_NOT_ADMITTED" || error("invalid parent comparison receipt")
        parent["sources"]["driver/check_native_precision.jl"]==precision_sources()["driver/check_native_precision.jl"] || error("different reference producer")
        parent["sources"]["driver/hoppet-contract.json"]==precision_sources()["driver/hoppet-contract.json"] || error("different frozen policy")
    end
    sources=precision_sources();build=CollinsFit.verified_hoppet_build(buildfile)
    buildhash=bytes2hex(sha256(read(buildfile)));cases=[];checks=[];nf_cases=[]
    check(name,error,tolerance)=push!(checks,Dict("name"=>name,"error"=>error,"tolerance"=>tolerance,"passed"=>isfinite(error)&&0<=error<=tolerance))
    result=Dict("status"=>"RUNNING","role"=>role,"scope"=>"parent_native_diagnostic_reference" ,"sources"=>sources,
        "contract_id"=>"KPSY-NATIVE-PRECISION-001","build_path"=>buildfile,"build_sha256"=>buildhash,
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,"cases"=>cases,"checks"=>checks,"nf_cases"=>nf_cases,
        "uniform_near_endpoint_h1_accuracy_certified"=>false,"synthetic_tests_grant_physics_admission"=>false)
    mkpath(dirname(output))
    try
        for construction in Symbol.(policy["constructions"]),setting in policy["constructor_cases"]
            key=string(construction)*"/"*setting["label"]
            row,diagnostic=capture_precision(output*".$(construction).$(setting["label"]).stderr.log") do
                native_case(buildfile,setting,construction)
            end
            row["native_diagnostics"]=diagnostic;push!(cases,row)
            if role=="candidate"
                baseline=filter(c->c["label"]==row["label"] && c["construction"]==row["construction"],parent["cases"])
                length(baseline)==1 || error("missing or duplicate parent case")
                prior=only(baseline);Set(keys(prior["payload"]))==Set(keys(row["payload"])) || error("payload fields changed")
                for name in sort!(collect(keys(row["payload"])))
                    a,b=row["payload"][name],Float64.(prior["payload"][name])
                    length(a)==length(b) || error("parent payload shape changed")
                    check(key*" parent exact "*name,reinterpret(UInt64,a)==reinterpret(UInt64,b) ? 0. : 1.,0.)
                end
                check(key*" zero native errors",diagnostic["precision_error_count"],policy["native_precision_errors_allowed"])
                check(key*" independent Mellin",row["mellin_max_scaled"],policy["point_Mellin_max_scaled"])
                check(key*" retained rows/forward",row["row_native_max_scaled"],policy["same_grid_row_native_max_scaled"])
                row["Collins_row_native_max_scaled"]===nothing || check(key*" retained Collins rows/forward",row["Collins_row_native_max_scaled"],policy["same_grid_row_native_max_scaled"])
            end
            CollinsFit.writejson(output*".progress.json",result)
            println(JSON.json(Dict("completed"=>key,"checks"=>length(checks),"native_precision_errors"=>diagnostic["precision_error_count"])));flush(stdout)
        end
        length(cases)==8 || error("incomplete constructor coverage")
        if role=="candidate"
            for nf in policy["nf_cases"]
                row,diagnostic=capture_precision(output*".nf$(nf).stderr.log") do
                    native=HoppetPlan(buildfile,[10.],2.4,AsymptoticTwoLoopCoupling(nf=nf);nf=nf,xmin=.1,dy=-log(.1)/512,layout=:single)
                    try
                        h=native.grid.*(1 .-native.grid).^3;v=apply_operator(native,h)[1,1]
                        reference=mellin_polynomial_tau(.1,native.tau[1])
                        Dict{String,Any}("nf"=>nf,"x"=>.1,"Q2_GeV2"=>10.,"Q02_GeV2"=>2.4,"tau"=>native.tau[1],
                            "actual"=>v,"reference"=>reference,"max_scaled"=>abs(v-reference)/max(abs(reference),policy["point_reference_floor"]))
                    finally close(native) end
                end
                row["native_diagnostics"]=diagnostic;push!(nf_cases,row)
                diagnostic["precision_error_count"]==0 || error("nf$(nf) native precision error")
                check("nf$(nf) independent Mellin",row["max_scaled"],policy["point_Mellin_max_scaled"])
            end
            check("source/build stability",sources==precision_sources() && bytes2hex(sha256(read(buildfile)))==buildhash && CollinsFit.verified_hoppet_build(buildfile)==build ? 0. : 1.,0.)
            result["status"]=all(c["passed"] for c in checks) ? "PASS" : "FAIL"
            result["scope"]="reconstructed_theory"
            result["parent_report_sha256"]=bytes2hex(sha256(read(parentfile)))
            result["parent_native_precision_errors_preserved"]=sum(c["native_diagnostics"]["precision_error_count"] for c in parent["cases"])
        else
            sources==precision_sources() && bytes2hex(sha256(read(buildfile)))==buildhash || error("parent sources/build changed")
            result["status"]="REFERENCE_COMPLETED_WITH_DIAGNOSTICS_NOT_ADMITTED"
        end
        result["completed_at"]=string(now(UTC));CollinsFit.writejson(output,result)
        result["status"]=="FAIL" && error("native precision comparison failed; preserve the failed report")
        result
    catch exception
        result["status"]="FAIL";result["error"]=sprint(showerror,exception);CollinsFit.writejson(output,result);rethrow()
    end
end

length(ARGS) in (3,4) || error("usage: parent|candidate BUILD OUTPUT [PARENT_REPORT]")
check_precision(ARGS...)
