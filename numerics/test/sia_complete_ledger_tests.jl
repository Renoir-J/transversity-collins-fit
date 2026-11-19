using Test, JSON
include("SIACompleteJobFixture.jl")
const F=SIACompleteJobFixture
const C=F.C
const ROOT=F.ROOT

function fresh_reader(directory,name)
    output=joinpath(directory,name*".json");log=joinpath(directory,name*".log")
    command=`$(Base.julia_cmd()) --startup-file=no --project=$(joinpath(ROOT,"numerics")) $(joinpath(@__DIR__,"sia_complete_fixture_probe.jl")) $directory $output`
    process=open(log,"w") do stream
        run(pipeline(ignorestatus(command),stdout=stream,stderr=stream))
    end
    process.exitcode,JSON.parsefile(output)
end

function reseal(c,jobs)
    for (phase,index) in jobs;F.seal(c,phase,index;verify=false);end
end

@testset "Complete79/115 saved Julia job fixtures and fresh full72 reader; no physical admission" begin
    for force4 in (false,true)
        base=joinpath(ROOT,"runs","sia-full-ledger-fixtures");mkpath(base)
        directory=mktempdir(base;prefix=force4 ? "native4-" : "native3-",cleanup=false)
        let directory=directory
            c,jobs=F.build(directory;force4=force4)
            @test length(jobs)==(force4 ? 115 : 79)
            @test count(first(job)=="native4" for job in jobs)==(force4 ? 36 : 0)
            @test C.needs4(c)===force4
            # The physical loader rejects this manifest even though its complete
            # mathematical and saved-job fixture can be reviewed below.
            @test_throws ErrorException C.load_manifest(c.path,c.sha)
            code,report=fresh_reader(directory,"fresh-positive")
            @test code==0
            @test report["status"]=="SYNTHETIC_FULL72_SAVED_LEDGER_REVIEWED_NOT_PHYSICAL_ADMISSION"
            @test report["theory_scope"]=="synthetic_closure" && report["fixture_scope"]==F.FIXTURE
            @test report["selected_native_level"]==(force4 ? 4 : 3)
            @test report["prerequisite_pairs"]==(force4 ? 114 : 78)
            @test report["source_point_ids"]==[p["id"] for p in c.source.request.data["points"]]
            @test all(report[key]===false for key in ("per_job_commands_measured","production_locks_acquired",
                "native_evolution_executed","physical_SIA_jobs_started","fit_performed","production_ready","accepted_fitting_gate"))
            legacy=C.load_library(ROOT)
            proof=Base.invokelatest(C.collect_dependencies,c,"diagnostic",legacy)
            raw=Base.invokelatest(C.collect_dependencies,c,"diagnostic",legacy;reuse_observed_file_hashes=false)
            @test isequal(proof.records,raw.records)
            @test isequal(proof.input_hashes,raw.input_hashes)
            @test isequal(proof.canonical_files,raw.canonical_files)
            observed=C.ObservedFileHashes(proof.input_hashes)
            @test_throws ErrorException C.observed_digest(joinpath(directory,"unobserved-file.json"),observed)
            if force4
                # A native4 packet directly names the native3 selector, whose
                # own native123 inputs must remain protected through the work.
                nested=Base.invokelatest(C.collect_dependencies,c,"native4",legacy)
                ancestor=joinpath(C.folder(c,"native123",1),"report.json")
                @test haskey(nested.input_hashes,ancestor)
                original=read(ancestor)
                open(ancestor,"a") do stream
                    write(stream,"\n")
                end
                @test_throws ErrorException C.hashes(nested.input_hashes)
                write(ancestor,original)
                @test C.hashes(nested.input_hashes)===true
            end
            if !force4
                file=joinpath(C.folder(c,"native123",1),"producer.log");old=read(file)
                write(file,"synthetic post-reader mutation\n")
                @test_throws ErrorException C.hashes(proof.input_hashes)
                write(file,old)
                @test C.hashes(proof.input_hashes)===true
                if Sys.islinux()
                    equivalent=joinpath(directory,"byte-identical-producer-copy.log")
                    write(equivalent,old);rm(file)
                    try
                        symlink(equivalent,file)
                        @test islink(file)
                        @test C.hashes(proof.input_hashes)===true
                        @test_throws ErrorException C.recheck_dependencies(proof)
                    finally
                        rm(file;force=true);write(file,old)
                    end
                    @test C.recheck_dependencies(proof)===true
                end
                auditfile=joinpath(C.folder(c,"native123",12),"audit.json")
                original=read(auditfile)
                for mutation in ("missing-audit-context","false-cold-evolution-claim","promoted-leaf-proof")
                    audit=JSON.parse(String(copy(original)))
                    if mutation=="missing-audit-context";delete!(audit,"source_context")
                    elseif mutation=="false-cold-evolution-claim";audit["native_evolution_reexecuted_by_reader"]=true
                    else;audit["leaf_arithmetic_review"]["job_resource_and_complete_provenance_ledger_verified"]=true
                    end
                    write(auditfile,JSON.json(audit)*"\n");reseal(c,jobs)
                    failed,result=fresh_reader(directory,mutation)
                    @test failed==2
                    @test result["status"]=="FAIL" && result["fit_performed"]===false
                    write(auditfile,original);reseal(c,jobs)
                end
                # Keep all local pair/category algebra valid and rehash every
                # protocol file. The old global common-order choice must still
                # be rejected by the independent saved-operand reducer.
                path=joinpath(C.folder(c,"native123",36),"report.json");original=read(path)
                changed=JSON.parse(String(copy(original)))
                leaf=changed["leaf_report"]
                points=leaf["native_values"][3]["momentum_levels"][3]["points"]
                sources=c.source.request.data["points"][leaf["definition"]["source_point_indices"]]
                for (point,source) in zip(points,sources)
                    for pair in Base.values(point["pairs"]), field in ("UU","Collins")
                        pair[field]*=1.001
                    end
                    pairs=point["pairs"]
                    u=Dict(field=>pairs["pm"][field]+pairs["mp"][field] for field in ("UU","Collins"))
                    l=Dict(field=>pairs["pp"][field]+pairs["mm"][field] for field in ("UU","Collins"))
                    categories=Dict("U"=>u,"L"=>l,"C"=>Dict(field=>u[field]+l[field] for field in ("UU","Collins")))
                    point["categories"]=categories
                    right=categories[source["ratio"]=="UL" ? "L" : "C"]
                    amplitude=source["source_angular"]["value"]*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"])
                    point["algebra_reference"]=amplitude;point["amplitude"]=amplitude
                end
                # Prove that this mutation reaches the global selector rather
                # than failing first on local sums, ratios or saved leaf proof.
                checked=C.N.audit(c.source.request.data,c.source.parent,leaf["definition"],c.binding,leaf)
                stored=JSON.parsefile(joinpath(C.folder(c,"native123",36),"audit.json"))
                @test C.S.same_json(checked,stored["leaf_arithmetic_review"])
                write(path,JSON.json(changed)*"\n");reseal(c,jobs)
                failed,result=fresh_reader(directory,"rehashed-global-momentum-substitution")
                @test failed==2 && result["status"]=="FAIL"
                @test occursin("changed first common all72 momentum decision",result["error"])
                write(path,original);reseal(c,jobs)
                restored,result=fresh_reader(directory,"fresh-restored")
                @test restored==0 && result["full72_saved_fixture_arithmetic_verified"]===true
            end
        end
    end
end
