using Test, CollinsFit, JSON

function captured_job_scope_error(f)
    try
        f()
        nothing
    catch error
        error
    end
end

function check_job_scope_error(error,label)
    @test error isa ErrorException
    @test error!==nothing && startswith(sprint(showerror,error),label*" requires explicit scope=\"synthetic_closure\"")
end

@testset "finite-kernel loader rejects unsupported or missing source scopes" begin
    root=normpath(joinpath(@__DIR__,"../.."))
    original_data=JSON.parsefile(joinpath(root,"fixtures/closure/data.json"))
    original_kernel=JSON.parsefile(joinpath(root,"fixtures/closure/kernel.json"))
    unsupported=("production","experimental_candidate","historical_benchmark_candidate","reconstructed_theory","",true,1,nothing)
    mktempdir() do directory
        datafile=joinpath(directory,"data.json");kernelfile=joinpath(directory,"kernel.json")
        for scope in unsupported
            data=deepcopy(original_data);data["scope"]=scope
            writejson(datafile,data);writejson(kernelfile,original_kernel)
            check_job_scope_error(captured_job_scope_error(()->load_context(datafile,kernelfile)),"data")
        end
        for scope in unsupported
            kernel=deepcopy(original_kernel);kernel["scope"]=scope
            writejson(datafile,original_data);writejson(kernelfile,kernel)
            check_job_scope_error(captured_job_scope_error(()->load_context(datafile,kernelfile)),"kernel")
        end
        for label in ("data","kernel")
            data=deepcopy(original_data);kernel=deepcopy(original_kernel)
            delete!(label=="data" ? data : kernel,"scope")
            writejson(datafile,data);writejson(kernelfile,kernel)
            check_job_scope_error(captured_job_scope_error(()->load_context(datafile,kernelfile)),label)
        end
        data=deepcopy(original_data);kernel=deepcopy(original_kernel)
        data["scope"]="production";kernel["scope"]="production"
        writejson(datafile,data);writejson(kernelfile,kernel)
        check_job_scope_error(captured_job_scope_error(()->load_context(datafile,kernelfile)),"data")
        missing_kernel=joinpath(directory,"must-not-open-kernel.json")
        check_job_scope_error(captured_job_scope_error(()->load_context(datafile,missing_kernel)),"data")
        @test !isfile(missing_kernel)
        writejson(datafile,original_data);writejson(kernelfile,Dict("scope"=>"production"))
        check_job_scope_error(captured_job_scope_error(()->load_context(datafile,kernelfile)),"kernel")
        for label in ("data","kernel")
            writejson(datafile,label=="data" ? [] : original_data)
            writejson(kernelfile,label=="kernel" ? [] : original_kernel)
            error=captured_job_scope_error(()->load_context(datafile,kernelfile))
            @test error isa ErrorException
            @test error!==nothing && startswith(sprint(showerror,error),label*" must be a JSON object")
        end
    end
end

@testset "complete numerical job requires scope before data IO or dispatch" begin
    root=normpath(joinpath(@__DIR__,"../.."))
    unsupported=("production","experimental_candidate","historical_benchmark_candidate","reconstructed_theory","",true,1,nothing)
    mktempdir() do directory
        output=joinpath(directory,"must-not-write.json")
        bad_data=joinpath(directory,"must-not-open-data.json")
        bad_kernel=joinpath(directory,"must-not-open-kernel.json")
        # Even the pre-fix implementation cannot fit: both inputs are absent.
        job=Dict{String,Any}("data"=>bad_data,"kernel"=>bad_kernel,"action"=>"fit","output"=>output,"seed"=>1729)
        for scope in unsupported
            request=copy(job);request["scope"]=scope
            check_job_scope_error(captured_job_scope_error(()->run_job(request)),"job")
            @test !isfile(output) && !isfile(output*".tmp")
        end
        check_job_scope_error(captured_job_scope_error(()->run_job(job)),"job")
        @test !isfile(output)
        for request in ([],nothing,true)
            error=captured_job_scope_error(()->run_job(request))
            @test error isa ErrorException
            @test error!==nothing && startswith(sprint(showerror,error),"numerical job must be a JSON object")
        end
        data=JSON.parsefile(joinpath(root,"fixtures/closure/data.json"));data["scope"]="production"
        datafile=joinpath(directory,"relabeled-data.json");writejson(datafile,data)
        request=Dict("scope"=>"synthetic_closure","data"=>datafile,"kernel"=>bad_kernel,
            "action"=>"scope_probe_no_numerics","output"=>output)
        check_job_scope_error(captured_job_scope_error(()->run_job(request)),"data")
        @test !isfile(output)
        valid=Dict("scope"=>"synthetic_closure","data"=>joinpath(root,"fixtures/closure/data.json"),
            "kernel"=>joinpath(root,"fixtures/closure/kernel.json"),"action"=>"benchmark","output"=>output)
        result=run_job(valid)
        @test result["scope"]=="synthetic_closure" && result["passed"]===true
        @test isfile(output) && !isfile(output*".tmp")
        @test JSON.parsefile(output)["scope"]=="synthetic_closure"
    end
end
