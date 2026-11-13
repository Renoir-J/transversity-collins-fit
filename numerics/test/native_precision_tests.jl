function precision_capture(f)
    mktemp() do _,stream
        value=redirect_stderr(stream) do
            try f() catch exception;exception end
        end
        flush(stream);seekstart(stream)
        (value,read(stream,String))
    end
end

@testset "native diagnostic rejection and stderr preservation" begin
    sentinel=[-0.,1.,2.]
    @test CollinsFit.hoppet_native_checked(()->sentinel)===sentinel
    @test signbit(CollinsFit.hoppet_native_checked(()->-0.))
    value,diagnostic=precision_capture() do
        CollinsFit.hoppet_native_checked("ordinary warning") do
            println(stderr,"ordinary native warning")
            19
        end
    end
    @test value==19
    @test diagnostic=="ordinary native warning\n"
    value,diagnostic=precision_capture() do
        CollinsFit.hoppet_native_checked("zero-status diagnostic") do
            if Sys.islinux()
                native_stderr=unsafe_load(cglobal(:stderr,Ptr{Cvoid}))
                ccall(:fputs,Cint,(Cstring,Ptr{Cvoid}),"cgauss64D113.1TOO HIGH ACCURACY REQUIRED\n",native_stderr)
            else
                println(stderr,"cgauss64D113.1TOO HIGH ACCURACY REQUIRED")
            end
            0
        end
    end
    @test value isa ErrorException
    @test occursin("native precision error",sprint(showerror,value))
    @test count("TOO HIGH ACCURACY REQUIRED",diagnostic)==1
    value,diagnostic=precision_capture() do
        CollinsFit.hoppet_native_checked("wrapped exception") do
            println(stderr,"preserve before exception")
            throw(ArgumentError("original failure"))
        end
    end
    @test value isa ArgumentError
    @test diagnostic=="preserve before exception\n"
    value,diagnostic=precision_capture() do
        CollinsFit.hoppet_native_checked("outer") do
            CollinsFit.hoppet_native_checked("inner") do
                println(stderr,"one nested warning")
                23
            end
        end
    end
    @test value==23
    @test diagnostic=="one nested warning\n"
    @test CollinsFit.hoppet_native_checked(()->29)==29
end

if Sys.islinux()
    @testset "synthetic zero-status native constructor closes rejected context" begin
        mktempdir() do directory
            compiler=Sys.which("cc")
            compiler===nothing && error("C compiler required for the native error fixture")
            source=normpath(joinpath(@__DIR__,"../../fixtures/native_precision/zero_status_error.c"))
            library=joinpath(directory,"libzero_status_error.so")
            run(`$compiler -std=c11 -Wall -Wextra -Werror -fPIC -shared $source -o $library`)
            handle=CollinsFit.Libdl.dlopen(library)
            try
                active=CollinsFit.Libdl.dlsym(handle,:cf_test_active)
                closes=CollinsFit.Libdl.dlsym(handle,:cf_test_closes)
                @test ccall(active,Cint,())==0
                # Deliberately internal token: this is not a fabricated physics build receipt.
                token=CollinsFit.HoppetBuildVerification(Dict{String,Any}("library"=>library))
                for construction in (:independent,:batched_probes)
                    before=ccall(closes,Cint,())
                    value,diagnostic=precision_capture() do
                        plan=HoppetPlan(token,[10.],2.4,AsymptoticTwoLoopCoupling();xmin=.1,layout=:single,construction=construction)
                        close(plan)
                        :incorrectly_accepted
                    end
                    @test value isa ErrorException
                    @test occursin("native precision error",sprint(showerror,value))
                    @test count("TOO HIGH ACCURACY REQUIRED",diagnostic)==1
                    @test ccall(active,Cint,())==0
                    @test ccall(closes,Cint,())==before+1
                end
            finally
                CollinsFit.Libdl.dlclose(handle)
            end
        end
    end
end
