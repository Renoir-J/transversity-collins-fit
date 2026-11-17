using Test
include("../bin/SIDISContinuationPaths.jl")
const Paths=SIDISContinuationPaths
include("../bin/HistoricalSIDISJobRecords.jl")
const original_string_method=which(HistoricalSIDISJobRecords.point_folder,Tuple{String,String,Int})
include("../bin/SIDISContinuationDispatch.jl")

@testset "Explicitly synthetic continuation path dispatch" begin
    mktempdir() do temp
        root=joinpath(temp,"new");original=joinpath(temp,"old")
        mkpath(root);mkpath(original)
        run=Paths.ContinuationRun(root,original)
        for phase in Paths.POINT_PHASES,i in 1:132
            expected=phase=="native123" && i<=13 ? original : root
            @test Paths.point_folder(run,phase,i)==joinpath(expected,phase,"p"*lpad(string(i),4,'0'))
            @test HistoricalSIDISJobRecords.point_folder(run,phase,i)==Paths.point_folder(run,phase,i)
        end
        @test which(HistoricalSIDISJobRecords.point_folder,Tuple{String,String,Int})===original_string_method
        @test HistoricalSIDISJobRecords.point_folder(root,"native123",1)==joinpath(root,"native123/p0001")
        ctx=(requestfile="SYNTHETIC_REQUEST",nativefile="SYNTHETIC_NATIVE",inputfile="SYNTHETIC_PARTON")
        deps=HistoricalSIDISJobRecords.global_dependencies(ctx,run,"native-select3")
        @test length(deps)==135
        for i in 1:132
            @test deps[i+3]==joinpath(Paths.point_folder(run,"native123",i),"invocation.json")
        end
        fourth=HistoricalSIDISJobRecords.point_dependencies(ctx,run,"native4",13)
        @test fourth[end]==joinpath(original,"native123/p0013/report.json")
        @test joinpath(run,"diagnostic","report.json")==joinpath(root,"diagnostic/report.json")
        @test_throws ErrorException Paths.point_folder(run,"native123",true)
        @test_throws ErrorException Paths.point_folder(run,"native123",0)
        @test_throws ErrorException Paths.point_folder(run,"native123",133)
        @test_throws ErrorException Paths.point_folder(run,"../native123",14)
        @test_throws ErrorException Paths.ContinuationRun(root,root)
        nested=joinpath(original,"nested");mkpath(nested)
        @test_throws ErrorException Paths.ContinuationRun(nested,original)
        output=joinpath(root,"audit.json");m=(run=run,)
        @test Paths.require_new_output(m,output)==output
        write(output,"synthetic occupied output")
        @test_throws ErrorException Paths.require_new_output(m,output)
        @test_throws ErrorException Paths.require_new_output(m,joinpath(original,"audit.json"))
        partial=joinpath(root,"partial.json");write(partial*".progress.json","synthetic partial")
        @test_throws ErrorException Paths.require_new_output(m,partial)
        alias=joinpath(root,"alias");symlink(original,alias)
        @test_throws ErrorException Paths.require_new_output(m,joinpath(alias,"audit.json"))
    end
end
