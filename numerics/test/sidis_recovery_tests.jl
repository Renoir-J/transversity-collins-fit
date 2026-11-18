using Test, JSON
include("../bin/SIDISRecoveryBindings.jl")
const RB=SIDISRecoveryBindings
const RP=RB.RP
include("../bin/HistoricalSIDISJobRecords.jl")
const SP=HistoricalSIDISJobRecords
const old_point_method=which(SP.point_folder,Tuple{String,String,Int})
const old_dependencies_method=which(SP.point_dependencies,Tuple{Any,String,String,Int})
const old_audit_method=which(SP.audit_point,Tuple{Any,String,String,Int})
include("../bin/SIDISRecoveryDispatch.jl")

function synthetic_views(temp)
    original=joinpath(temp,"original/full");predecessor=joinpath(temp,"predecessor");current=joinpath(temp,"current")
    foreach(mkpath,(original,predecessor,current))
    first=joinpath(dirname(original),"invocation.json");second=joinpath(predecessor,"continuation.json")
    for path in (first,second);write(path,"{\"scope\":\"explicitly_synthetic_location_fixture_not_a_job\"}");end
    rows=Dict{String,Any}[]
    for index in 1:37
        producer,manifest,protocol=index<=13 ? (original,first,"frozen_original") : (predecessor,second,"continuation_v1")
        folder=joinpath(producer,"native123","p"*lpad(string(index),4,'0'));mkpath(folder)
        push!(rows,RB.entry("native123",index,folder,producer,manifest,RB.digest(manifest),protocol))
    end
    legacy=RB.LegacyRun(predecessor,original,Dict{String,Any}("scope"=>RB.CP.SCOPE,"manifest_sha256"=>RB.digest(second),
        "adapter_commit"=>"SYNTHETIC_NOT_A_COMMIT","contract_sha256"=>RB.CP.CONTRACT_SHA))
    run=RB.RecoveryRun(RP.RecoveryRunPaths(current,rows),(legacy,),Dict{String,Any}("scope"=>RB.SCOPE,"manifest_sha256"=>"SYNTHETIC_NOT_A_MANIFEST"))
    (;original,predecessor,current,rows,legacy,run)
end

@testset "Synthetic complete-source recovery dispatch preserves actual producing paths" begin
    mktempdir() do temp
        s=synthetic_views(temp)
        ctx=(requestfile="SYNTHETIC_REQUEST",nativefile="SYNTHETIC_NATIVE",inputfile="SYNTHETIC_PARTON")
        for phase in RP.POINT_PHASES,index in 1:132
            expected=phase=="native123" && index<=13 ? s.original : phase=="native123" && index<=37 ? s.predecessor : s.current
            @test SP.point_folder(s.run,phase,index)==joinpath(expected,phase,"p"*lpad(string(index),4,'0'))
        end
        deps=SP.global_dependencies(ctx,s.run,"native-select3")
        @test length(deps)==135
        for index in 1:132
            @test deps[index+3]==joinpath(SP.point_folder(s.run,"native123",index),"invocation.json")
        end
        @test SP.point_dependencies(ctx,s.run,"native123",14)==[ctx.requestfile,ctx.nativefile,ctx.inputfile]
        @test SP.point_dependencies(ctx,s.run,"native4",1)[end-1:end]==[
            joinpath(s.current,"native-select3/report.json"),joinpath(s.original,"native123/p0001/report.json")]
        @test RB.origin(s.run,"native123",1)==s.original
        @test RB.origin(s.run,"native123",14)===s.legacy
        @test RB.origin(s.run,"native123",38)===s.run
        @test which(SP.point_folder,Tuple{String,String,Int})===old_point_method
        @test which(SP.point_dependencies,Tuple{Any,String,String,Int})===old_dependencies_method
        @test which(SP.audit_point,Tuple{Any,String,String,Int})===old_audit_method
        @test SP.point_folder(s.current,"native123",1)==joinpath(s.current,"native123/p0001")
        @test check_recovery_transport(s.run,"native123",1,Dict())
        @test_throws ErrorException check_recovery_transport(s.run,"native123",1,Dict("recovery"=>s.run.transport))
        @test check_recovery_transport(s.run,"native123",14,Dict("continuation"=>s.legacy.transport))
        @test_throws ErrorException check_recovery_transport(s.run,"native123",14,Dict("continuation"=>Dict()))
        @test check_recovery_transport(s.run,"native123",38,Dict("recovery"=>s.run.transport))
        @test_throws ErrorException check_recovery_transport(s.run,"native123",38,Dict("continuation"=>s.legacy.transport))
        @test RP.metadata(s.run.paths)["saved_numerical_arithmetic_verified"]===false
    end
end

@testset "Synthetic later-phase recovery keeps old prerequisite evidence and new current decisions separate" begin
    mktempdir() do temp
        s=synthetic_views(temp)
        manifest=joinpath(s.current,"recovery.json");write(manifest,"{\"scope\":\"explicitly_synthetic_older_namespace\"}")
        folder=joinpath(s.current,"native4/p0001");mkpath(folder)
        rows=vcat(s.rows,[RB.entry("native4",1,folder,s.current,manifest,RB.digest(manifest),"recovery_v2")])
        latest=joinpath(temp,"latest");mkpath(latest)
        current=RB.RecoveryRun(RP.RecoveryRunPaths(latest,rows),(s.legacy,s.run),Dict{String,Any}("scope"=>RB.SCOPE))
        ctx=(requestfile="SYNTHETIC_REQUEST",nativefile="SYNTHETIC_NATIVE",inputfile="SYNTHETIC_PARTON")
        old_deps=SP.point_dependencies(ctx,current,"native4",1)
        @test old_deps[end-1]==joinpath(s.current,"native-select3/report.json")
        @test old_deps[end]==joinpath(s.original,"native123/p0001/report.json")
        new_deps=SP.point_dependencies(ctx,current,"native4",2)
        @test new_deps[end-1]==joinpath(latest,"native-select3/report.json")
        @test new_deps[end]==joinpath(s.original,"native123/p0002/report.json")
        @test RB.origin(current,"native4",1)===s.run
        @test RB.origin(current,"native4",2)===current
        @test_throws ErrorException RB.require_new_output((run=current,),joinpath(folder,"fresh.json"))
        @test RB.require_new_output((run=current,),joinpath(latest,"fresh.json"))==joinpath(latest,"fresh.json")
        partial=joinpath(latest,"partial.json");write(partial*".progress.json","SYNTHETIC_PARTIAL")
        @test_throws ErrorException RB.require_new_output((run=current,),partial)
        @test_throws ErrorException RP.fresh_folder(current.paths,"native4",1)
    end
end

struct SyntheticCurrentPrerequisites
    data::Dict{String,Any}
end
SP.selection(ctx,run::SyntheticCurrentPrerequisites,phase)=run.data[phase]

@testset "Synthetic retained packets must agree with fresh current selections" begin
    choices=Dict{String,Any}("native-select3"=>Dict("status"=>"NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED","selected_level"=>0),
        "native-selection"=>Dict("selected_level"=>4),
        "derivative-select3"=>Dict("status"=>"DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED","selected_reference_level"=>0))
    run=SyntheticCurrentPrerequisites(choices)
    @test check_current_recovery_prerequisites(nothing,run,"native4",Dict())
    @test check_current_recovery_prerequisites(nothing,run,"tailref123",Dict("native_level"=>4))
    @test_throws ErrorException check_current_recovery_prerequisites(nothing,run,"tailref123",Dict("native_level"=>3))
    @test check_current_recovery_prerequisites(nothing,run,"reference4",Dict("native_level"=>4))
    choices["derivative-select3"]["selected_reference_level"]=3
    @test_throws ErrorException check_current_recovery_prerequisites(nothing,run,"reference4",Dict("native_level"=>4))
    choices["native-select3"]["selected_level"]=3
    @test_throws ErrorException check_current_recovery_prerequisites(nothing,run,"native4",Dict())
    @test which(SP.point_dependencies,Tuple{Any,String,String,Int})===old_dependencies_method
end

@testset "Synthetic immutable inventories preserve links and reject stale or failed jobs" begin
    mktempdir() do temp
        old=joinpath(temp,"old");outside=joinpath(temp,"outside");mkpath(old);mkpath(outside)
        record=joinpath(old,"a.json");write(record,"SYNTHETIC_IMMUTABLE")
        write(joinpath(outside,"outside.json"),"OUTSIDE_NOT_READ")
        symlink(outside,joinpath(old,"old-link"))
        rows=RB.inventory(old)
        @test length(rows)==2
        @test RB.verify_inventory(old,rows,temp)==rows
        write(record,"CHANGED")
        @test_throws ErrorException RB.verify_inventory(old,rows,temp)
        @test read(joinpath(outside,"outside.json"),String)=="OUTSIDE_NOT_READ"
        folder=joinpath(old,"native123/p0001");mkpath(folder)
        inv=Dict{String,Any}("status"=>RB.COMPLETED,"phase"=>"native123","point_index"=>1,"completed_at"=>"SYNTHETIC_ONLY",
            "steps"=>[Dict{String,Any}("name"=>name,"returncode"=>0,"process_returncode"=>0,"interruption"=>nothing,"termination_signal"=>nothing) for name in ("producer","audit")])
        path=joinpath(folder,"invocation.json");write(path,JSON.json(inv));touch(joinpath(folder,"job.completed"))
        @test length(RB.complete_locations(old,"recovery_v2",joinpath(old,"recovery.json"),"0"^64))==1
        inv["steps"][1]["interruption"]="SYNTHETIC_TIMEOUT";write(path,JSON.json(inv))
        @test_throws ErrorException RB.complete_locations(old,"recovery_v2",joinpath(old,"recovery.json"),"0"^64)
        inv["status"]="FAIL";inv["completed_at"]=nothing;write(path,JSON.json(inv));rm(joinpath(folder,"job.completed"))
        @test isempty(RB.complete_locations(old,"recovery_v2",joinpath(old,"recovery.json"),"0"^64))
        inv["completed_at"]="SYNTHETIC_TERMINAL_FAILURE";write(path,JSON.json(inv))
        @test_throws ErrorException RB.complete_locations(old,"recovery_v2",joinpath(old,"recovery.json"),"0"^64)
        @test_throws ErrorException RB.ordered([Dict("phase"=>"native123","index"=>2)])
        @test_throws ErrorException RB.ordered([Dict("phase"=>"native123","index"=>1),Dict("phase"=>"native123","index"=>1)])
    end
end

@testset "Recovery entrypoints parse without numerical dispatch" begin
    for name in ("run_sidis_recovery_point.jl","run_sidis_recovery_global.jl","audit_sidis_recovery_job.jl","audit_sidis_recovery_retained.jl")
        @test Meta.parseall(read(joinpath(@__DIR__,"../bin",name),String)) isa Expr
    end
end
