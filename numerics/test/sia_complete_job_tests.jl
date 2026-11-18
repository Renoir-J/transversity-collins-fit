using Test, JSON, CollinsFit
include("../bin/SIACompleteJobs.jl")
const C=SIACompleteJobs
module SyntheticCompleteOperands
include("sia_geometry_global_tests.jl")
end
const F=SyntheticCompleteOperands

# Load the actual CLI in an isolated namespace. This catches entrypoint
# top-level errors that loading only the numerical library would miss.
module IsolatedReadOnlyEntrypoint
include("../bin/inspect_sia_geometry_execution.jl")
end

@testset "Read-only source inspection entrypoint loads independently" begin
    @test isdefined(IsolatedReadOnlyEntrypoint,:main)
    @test_throws ErrorException IsolatedReadOnlyEntrypoint.main(String[])
end

@testset "Complete global worker arithmetic in an isolated module, synthetic only" begin
    mktempdir() do directory
        parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"))
        data=F.F.fixture_request();reports=F.native_reports(data,parent)
        c=(source=(request=(data=data,),parent=parent),m=Dict("run"=>directory,"limitations"=>["Explicitly synthetic operand fixture; not physical validation"]))
        records=Dict{Tuple{String,Int},Any}(("native123",i)=>Dict("leaf_report"=>reports[1][i]) for i in 1:36)
        coarse=C.produce_global(c,"native-select3",records)
        C.audit_global(c,"native-select3",records,coarse)
        @test coarse["selection"]["selected_level"]==3
        records[("native-select3",0)]=coarse
        path=joinpath(directory,"global/native-select3/report.json");mkpath(dirname(path));writejson(path,coarse)
        native=C.produce_global(c,"native-selection",records)
        C.audit_global(c,"native-selection",records,native)
        @test native["selection"]["selected_level"]==3
        records[("native-selection",0)]=native
        @test C.selected_native(c,records).selected==3
        bad=deepcopy(native);bad["selection"]["selected_level"]=4
        @test_throws ErrorException C.audit_global(c,"native-selection",records,bad)
        @test_throws ErrorException C.folder(c,"derivative",2)
        @test C.folder(c,"native123",1)==joinpath(directory,"native123/g0001")
        @test_throws ErrorException C.folder(c,"native123",true)
        @test C.dependencies(c,"native-selection")==vcat([("native123",i) for i in 1:36],[("native-select3",0)])
        @test length(C.dependencies(c,"diagnostic"))==78
        incomplete=copy(records);delete!(incomplete,("native123",36))
        @test_throws KeyError C.produce_global(c,"native-select3",incomplete)
    end
end

@testset "Saved complete GNU-time command pair is not a status or marker" begin
    mktempdir() do directory
        budget=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-job-contract.json"))["resources"]["global"]
        commands=[["explicitly-synthetic-command-not-executed",name] for name in ("producer","audit")]
        steps=Any[]
        for (name,command) in zip(("producer","audit"),commands)
            log=joinpath(directory,name*".log");resource=joinpath(directory,name*".resources.txt")
            write(log,"Explicitly synthetic saved measurement fixture; no physical execution\n")
            write(resource,"COLLINSFLOW_TIME_V1\t100\t1.0\t0.5\t0.1\t0\n")
            push!(steps,Dict{String,Any}("name"=>name,"command"=>command,"returncode"=>0,"process_returncode"=>0,
                "interruption"=>nothing,"termination_signal"=>nothing,"wrapper_elapsed_seconds"=>1.1,
                "minimum_available_memory_observed_before_command"=>2^40,"resource"=>C.resource(resource),
                "log_sha256"=>C.digest(log),"resource_sha256"=>C.digest(resource)))
        end
        record=Dict("resources"=>budget,"steps"=>steps)
        @test C.steps(record,budget,commands,directory)===nothing
        for mutation in ("effective","raw","bool","signal","rss","duration","memory","command","missing","budget","measurement")
            bad=deepcopy(record);step=first(bad["steps"])
            if mutation=="effective";step["returncode"]=124
            elseif mutation=="raw";step["process_returncode"]=2
            elseif mutation=="bool";step["returncode"]=false
            elseif mutation=="signal";step["interruption"]="parent_signal"
            elseif mutation=="rss";step["resource"]["maximum_RSS_KiB"]=2^40
            elseif mutation=="duration";step["wrapper_elapsed_seconds"]=1e9
            elseif mutation=="memory";step["minimum_available_memory_observed_before_command"]=1
            elseif mutation=="command";push!(step["command"],"substitution")
            elseif mutation=="missing";pop!(bad["steps"])
            elseif mutation=="budget";bad["resources"]["maximum_command_RSS_KiB"]*=2
            elseif mutation=="measurement";step["resource"]["child_exit_code"]=1
            end
            @test_throws ErrorException C.steps(bad,budget,commands,directory)
        end
        write(joinpath(directory,"audit.resources.txt"),"COLLINSFLOW_TIME_V1\t100\t1\t.5\t.1\t0\nCOLLINSFLOW_TIME_V1\t100\t1\t.5\t.1\t0\n")
        @test_throws ErrorException C.resource(joinpath(directory,"audit.resources.txt"))
    end
end

@testset "Worker entrypoints reject missing manifest identity without dispatch" begin
    @test_throws ErrorException C.main(String[])
    @test_throws ErrorException C.main(String[];audit=true)
    mktempdir() do directory
        file=joinpath(directory,"execution.json");write(file,"{}\n")
        @test_throws ErrorException C.load_manifest(file,"0"^64)
        @test_throws KeyError C.load_manifest(file,C.digest(file))
    end
end
