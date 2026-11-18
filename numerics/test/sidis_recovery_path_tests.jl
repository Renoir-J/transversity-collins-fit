using Test, JSON, SHA
const original_string_join_method=which(joinpath,Tuple{String,String})
include("../bin/SIDISRecoveryPaths.jl")
const R=SIDISRecoveryPaths

function synthetic_recovery_locations(temp)
    original=joinpath(temp,"original/full");predecessor=joinpath(temp,"predecessor");current=joinpath(temp,"current")
    foreach(mkpath,(original,predecessor,current))
    first_manifest=joinpath(dirname(original),"invocation.json");second_manifest=joinpath(predecessor,"continuation.json")
    for path in (first_manifest,second_manifest);write(path,"{\"scope\":\"explicitly_synthetic_location_fixture_not_a_job\"}\n");end
    rows=Any[]
    for index in 1:37
        producer,manifest,protocol=index<=13 ? (original,first_manifest,"frozen_original") : (predecessor,second_manifest,"continuation_v1")
        folder=joinpath(producer,"native123","p"*lpad(string(index),4,'0'));mkpath(folder)
        push!(rows,Dict{String,Any}("phase"=>"native123","index"=>index,"folder"=>folder,"producer_root"=>producer,
            "producer_manifest"=>manifest,"producer_manifest_sha256"=>R.digest(manifest),"protocol"=>protocol))
    end
    partial=joinpath(predecessor,"native123/p0038/report.json");mkpath(dirname(partial));write(partial,"explicitly synthetic incomplete point38")
    current,original,predecessor,rows,partial
end

@testset "Synthetic owned recovery locations retain37 paths without admitting jobs" begin
    mktempdir() do temp
        current,original,predecessor,rows,partial=synthetic_recovery_locations(temp)
        partial_sha=R.digest(partial);run=R.RecoveryRunPaths(current,rows)
        for phase in R.POINT_PHASES,index in 1:132
            expected=phase=="native123" && index<=13 ? original : phase=="native123" && index<=37 ? predecessor : current
            @test R.folder(run,phase,index)==joinpath(expected,phase,"p"*lpad(string(index),4,'0'))
        end
        for phase in R.GLOBAL_PHASES;@test R.folder(run,phase)==joinpath(current,phase);end
        @test R.fresh_folder(run,"native123",38)==joinpath(current,"native123/p0038")
        @test R.retained(run,"native123",38)===nothing && R.digest(partial)==partial_sha
        @test which(joinpath,Tuple{String,String})===original_string_join_method
        @test joinpath(run,"diagnostic/report.json")==joinpath(current,"diagnostic/report.json")
        metadata=R.metadata(run)
        @test metadata["scope"]=="owned_whole_job_path_view_not_admission"
        @test all(metadata[key]===false for key in (R.FLAGS...,"numerical_dispatch_ready","complete_job_provenance_verified","saved_numerical_arithmetic_verified"))
        rows[1]["index"]=999
        @test endswith(R.folder(run,"native123",1),"p0001")
        returned=R.retained(run,"native123",1);returned["folder"]="changed"
        @test endswith(R.folder(run,"native123",1),"p0001")
        for (phase,index) in (("native123",true),("native123",0),("native123",133),("diagnostic",true),("diagnostic",1),("unknown",0))
            @test_throws ErrorException R.folder(run,phase,index)
        end
        baseline=deepcopy(run.entries)
        for mutation in ("boolean","float","phase","missing-field","admission-field","protocol","point-folder","producer-root","sha","sha-type","order","duplicate")
            bad=deepcopy(baseline)
            if mutation=="boolean";bad[1]["index"]=true
            elseif mutation=="float";bad[1]["index"]=1.
            elseif mutation=="phase";bad[1]["phase"]="../native123"
            elseif mutation=="missing-field";delete!(bad[1],"protocol")
            elseif mutation=="admission-field";bad[1]["accepted_fitting_gate"]=true
            elseif mutation=="protocol";bad[1]["protocol"]="unreviewed"
            elseif mutation=="point-folder";bad[end]["folder"]=bad[end-1]["folder"]
            elseif mutation=="producer-root";bad[end]["producer_root"]=current
            elseif mutation=="sha";bad[end]["producer_manifest_sha256"]="0"^64
            elseif mutation=="sha-type";bad[end]["producer_manifest_sha256"]=true
            elseif mutation=="order";reverse!(bad)
            elseif mutation=="duplicate";push!(bad,deepcopy(bad[end]))
            end
            @test_throws ErrorException R.RecoveryRunPaths(current,bad)
        end
        @test_throws ErrorException R.fresh_folder(run,"native123",14)
        nested=joinpath(dirname(original),"other-output");mkpath(nested)
        @test_throws ErrorException R.RecoveryRunPaths(nested,baseline)
        alias=joinpath(temp,"alias");symlink(predecessor,alias)
        bad=deepcopy(baseline);bad[end]["folder"]=joinpath(alias,"native123/p0037")
        @test_throws ErrorException R.RecoveryRunPaths(current,bad)
        @test_throws ErrorException joinpath(run,"..","escaped")
        symlink(predecessor,joinpath(current,"alias"))
        @test_throws ErrorException joinpath(run,"alias/new-output")
        failure=joinpath(current,"native123/p0038/report.json");mkpath(dirname(failure));write(failure,"explicitly synthetic occupied failure")
        failure_sha=R.digest(failure)
        @test_throws ErrorException R.fresh_folder(run,"native123",38)
        @test R.digest(failure)==failure_sha
        manifest=joinpath(predecessor,"continuation.json");bytes=read(manifest);write(manifest,vcat(bytes,UInt8[10]))
        @test_throws ErrorException R.valid(run)
        write(manifest,bytes);run.entries[end]["index"]=36
        @test_throws ErrorException R.valid(run)
    end
end

@testset "Synthetic later-phase locations preserve actual producing roots" begin
    mktempdir() do temp
        current,_,_,rows,_=synthetic_recovery_locations(temp)
        producer=joinpath(temp,"older-recovery");mkpath(producer)
        manifest=joinpath(producer,"recovery.json");write(manifest,"{\"scope\":\"explicitly_synthetic_recovery_location_fixture\"}\n")
        for phase in R.POINT_PHASES[2:end]
            folder=joinpath(producer,phase,"p0002");mkpath(folder)
            push!(rows,Dict{String,Any}("phase"=>phase,"index"=>2,"folder"=>folder,"producer_root"=>producer,
                "producer_manifest"=>manifest,"producer_manifest_sha256"=>R.digest(manifest),"protocol"=>"recovery_v2"))
        end
        run=R.RecoveryRunPaths(current,rows)
        for phase in R.POINT_PHASES[2:end]
            @test R.folder(run,phase,2)==joinpath(producer,phase,"p0002")
            @test_throws ErrorException R.fresh_folder(run,phase,2)
        end
        @test R.metadata(run)["saved_numerical_arithmetic_verified"]===false
    end
end
