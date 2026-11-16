using SHA

# Explicitly synthetic fixed rows. These fixtures do not test native evolution.
function storage_synthetic_operator(process,kind,queries,radial,model)
    grid=sort(unique(vcat(.05,.1,queries,.75,1.)))
    q2=css_mu_b2.(Ref(model),radial.b)
    tau=evolution_time.(Ref(model.coupling),q2,model.Q02)
    alpha=alpha_s.(Ref(model.coupling),q2)
    blocks=CollinsFit.HoppetPointBlock[];saved=Dict{Float64,CollinsFit.FixedRowMatrix}()
    for query in queries
        rows=get!(saved,query) do
            values=zeros(length(q2),length(grid))
            values[:,findfirst(==(query),grid)].=1 .+ .015 .* radial.b
            CollinsFit.fixed_rows(values)
        end
        push!(blocks,CollinsFit.HoppetPointBlock(collect(eachindex(grid)),rows))
    end
    data=(grid,collect(queries),q2,model.Q02,tau,kind,5,32,:synthetic_storage_fixture,
        blocks,sum(sizeof(r) for r in values(saved)),0.)
    HoppetMatchingPlan(data...,process,alpha,CollinsFit.hoppet_point_snapshot(data...;process=process,alpha=alpha))
end

function storage_synthetic_context(process;shared=false,target=NucleonTarget(:proton),metadata_extra=Dict{String,Any}())
    model=KPSYCSS(order=16);radial=RadialRule([0.,.25,.8,2.];order=4)
    z2=process==:SIA ? (shared ? .3 : .5) : nothing
    plan=leptonic_plan(process,model,radial,process==:SIDIS ? 4. : 110.,.2,.3;z2=z2)
    first=[(1+.03i)*(1+.02b) for i in 1:10,b in radial.b]
    second=[(1+.05i)*(1+.04b) for i in 1:10,b in radial.b]
    kind=process==:SIDIS ? :transversity : :KPSY_Hhat3_homogeneous
    op=storage_synthetic_operator(process,kind,process==:SIDIS ? [.18] : [.3,z2],radial,model)
    bases=[(1+.1i)*x*(1-x)^3 for i in 1:(process==:SIDIS ? 2 : 4),x in op.grid]
    samples=PolarizedInitialSamples(kind,op.grid,bases;input_kind=:analytic)
    op2=process==:SIDIS ? storage_synthetic_operator(process,:KPSY_Hhat3_homogeneous,[.3],radial,model) : nothing
    s2=op2===nothing ? nothing : PolarizedInitialSamples(:KPSY_Hhat3_homogeneous,op2.grid,
        [(1+.1i)*x*(1-x)^3 for i in 1:4,x in op2.grid];input_kind=:analytic)
    metadata=Dict{String,Any}("scope"=>"reconstructed_theory","synthetic_storage_fixture"=>true,
        "literal_metadata"=>(-0.,UInt64(17),Int64(-3),false,nothing,:synthetic_storage),
        "bool_array"=>Bool[false,true],"uint_array"=>UInt64[0,typemax(UInt64)],
        "typed_metadata"=>[Dict("count"=>Int64(2))],"tuple_metadata"=>[(Int64(3),"synthetic")])
    merge!(metadata,metadata_extra)
    polarized_leptonic_plan(plan,first,second,op,samples;second_operator=op2,second_samples=s2,
        target=process==:SIDIS ? target : nothing,charges=process==:SIDIS ? (-1,0) : (1,-1),provenance=metadata)
end

function storage_synthetic_bindings()
    c=fixed_context_storage_contract()
    values=Dict(k=>bytes2hex(sha256("explicitly-synthetic-storage-binding:"*k)) for k in c["required_binding_fields"])
    values["storage_contract_sha256"]=CollinsFit.FIXED_STORAGE_CONTRACT_SHA
    values
end

storage_test_digest(path)=open(sha256,path) |> bytes2hex
storage_test_manifest(dir)=joinpath(dir,"fixed-context-manifest.json")
storage_test_payload(dir)=joinpath(dir,"fixed-context-payload.bin")
function storage_test_rewrite(dir,m)
    open(storage_test_manifest(dir),"w") do io;write(io,JSON.json(m));end
    storage_test_digest(storage_test_manifest(dir))
end
storage_test_load(dir,receipt,b)=load_fixed_leptonic_context(dir;expected_manifest_sha256=receipt.manifest_sha256,bindings=b,input_kind=:analytic)
storage_test_bits(x::Float64)=reinterpret(UInt64,x)
storage_test_bits(x::AbstractArray{Float64})=reinterpret(UInt64,vec(x))

@testset "synthetic fixed-context storage preserves full values,AD and aliases" begin
    b=storage_synthetic_bindings();theta=collect(kpsy_published_parameters().values)
    direction=[.02,.03,.01,.02,.01,.01,.02,-.02,.01,.02,.03,.03,.001]
    targets=[NucleonTarget(:proton),NucleonTarget(:neutron),
        NucleonTarget("synthetic free isoscalar",(.5,.5),(.5,.5);source="synthetic storage assembly fixture")]
    contexts=vcat([storage_synthetic_context(:SIDIS;target=t) for t in targets],
        [storage_synthetic_context(:SIA),storage_synthetic_context(:SIA;shared=true)])
    mktempdir() do temporary
        for (i,cold) in enumerate(contexts)
            dir=joinpath(temporary,"context-$i")
            receipt=save_fixed_leptonic_context(cold,dir;bindings=b,input_kind=:analytic)
            loaded=storage_test_load(dir,receipt,b)
            @test fixed_leptonic_storage_signature(loaded)==fixed_leptonic_storage_signature(cold)
            @test loaded!==cold
            @test loaded.first_operator.blocks[1].rows!==cold.first_operator.blocks[1].rows
            @test loaded.first_samples!==cold.first_samples
            @test loaded.provenance==cold.provenance
            @test signbit(loaded.provenance["literal_metadata"][1])
            @test typeof(loaded.provenance["typed_metadata"])==typeof(cold.provenance["typed_metadata"])
            @test receipt.payload_bytes==filesize(storage_test_payload(dir))
            @test receipt.input_kind===:analytic
            @test_throws ErrorException save_fixed_leptonic_context(cold,dir;bindings=b,input_kind=:analytic)
            for context in (cold,loaded)
                @test context.provenance["production_ready"]===false
                @test context.first_samples.input_kind===:analytic
            end
            if i==5
                @test cold.first_operator.blocks[1].rows===cold.first_operator.blocks[2].rows
                @test loaded.first_operator.blocks[1].rows===loaded.first_operator.blocks[2].rows
            end
            formulas=if cold.leptonic.process==:SIDIS
                [c->t->polarized_structure(c,t).Collins/polarized_structure(c,t).UU]
            else
                [c->t->predict_observable(sia_observable(c;ratio=r,
                    angular=SIAMeanOfRatio(.685;source="synthetic storage angle",averaging_measure="synthetic fixed weights"),
                    source="synthetic storage"),t) for r in (:UL,:UC)]
            end
            for build in formulas
                before,after=build(cold),build(loaded)
                for t in (theta,theta.+.01direction)
                    @test storage_test_bits(before(t))==storage_test_bits(after(t))
                    g1,g2=ForwardDiff.gradient(before,t),ForwardDiff.gradient(after,t)
                    h1=vec(ForwardDiff.jacobian(s->ForwardDiff.gradient(before,t+s[1]*direction),[0.]))
                    h2=vec(ForwardDiff.jacobian(s->ForwardDiff.gradient(after,t+s[1]*direction),[0.]))
                    @test storage_test_bits(g1)==storage_test_bits(g2)
                    @test storage_test_bits(h1)==storage_test_bits(h2)
                    cold.leptonic.process==:SIA && (@test all(iszero,g2[1:6]);@test all(iszero,h2[1:6]))
                end
            end
            rows=loaded.first_operator.blocks[1].rows
            @test_throws ErrorException (rows[1,1]=1.)
            @test_throws ErrorException (view(rows,:,1).=0.)
            copy_of_rows=copy(rows);copy_of_rows[1]+=1
            @test rows!=copy_of_rows
            # Once loaded, the owned matrices do not follow later disk edits.
            original=read(storage_test_payload(dir));changed=copy(original);changed[end]⊻=0x01
            write(storage_test_payload(dir),changed)
            @test fixed_leptonic_storage_signature(loaded)==fixed_leptonic_storage_signature(cold)
            @test_throws ErrorException storage_test_load(dir,receipt,b)
            write(storage_test_payload(dir),original)
            @test fixed_leptonic_storage_signature(storage_test_load(dir,receipt,b))==fixed_leptonic_storage_signature(cold)
            loaded.first_samples.bases[1,1]+=.01
            @test_throws ErrorException polarized_structure(loaded,theta)
        end
    end
end

@testset "storage rejects altered files,scope,bindings and typed structure" begin
    b=storage_synthetic_bindings();cold=storage_synthetic_context(:SIA;shared=true)
    mktempdir() do temporary
        dir=joinpath(temporary,"context");receipt=save_fixed_leptonic_context(cold,dir;bindings=b,input_kind=:analytic)
        original_manifest_bytes=read(storage_test_manifest(dir))
        original=JSON.parsefile(storage_test_manifest(dir));payload=read(storage_test_payload(dir))
        @test_throws ErrorException load_fixed_leptonic_context(dir;expected_manifest_sha256=repeat("0",64),bindings=b,input_kind=:analytic)
        @test_throws ErrorException load_fixed_leptonic_context(dir;expected_manifest_sha256=receipt.manifest_sha256,bindings=b,input_kind=:pinned_inputs)
        @test_throws ErrorException load_fixed_leptonic_context(dir;expected_manifest_sha256=receipt.manifest_sha256,bindings=b,input_kind=:production)
        for key in keys(b)
            altered=copy(b);altered[key]=repeat("1",64)
            @test_throws ErrorException load_fixed_leptonic_context(dir;expected_manifest_sha256=receipt.manifest_sha256,bindings=altered,input_kind=:analytic)
        end
        counterfactuals=[
            m->(m["production_ready"]=true),m->(m["fit_performed"]=true),
            m->(m["scope"]="production_fit"),m->(m["producer_julia_version"]="unavailable-version"),
            m->(m["payload"]["file"]="../unrelated.bin"),
            m->(m["native_endian"]=m["native_endian"]+1),m->(m["word_size"]=32),m->(m["word_size"]=64.),
            m->(m["root"]=Dict("ref"=>length(m["objects"])+1)),
            m->(m["objects"][1]["name"]="CT10Input"),
            m->(pop!(m["objects"][1]["fields"])),
            m->(first(o for o in m["objects"] if o["kind"]=="array")["offset"]+=1),
            m->(first(o for o in m["objects"] if o["kind"]=="array")["shape"]=[typemax(Int),typemax(Int)]),
            m->(first(o for o in m["objects"] if o["kind"]=="array")["type"]="Float32"),
            m->(push!(m["objects"],Dict("kind"=>"component","name"=>"unreachable"))),
            m->(first(o for o in m["objects"] if o["kind"]=="dict")["values"][1]=Dict("ref"=>1)),
        ]
        for change in counterfactuals
            m=deepcopy(original);change(m);new_digest=storage_test_rewrite(dir,m)
            @test_throws ErrorException storage_test_load(dir,receipt,b)
            @test_throws Exception load_fixed_leptonic_context(dir;expected_manifest_sha256=new_digest,bindings=b,input_kind=:analytic)
        end
        # A rehashed payload and manifest still fail the retained logical identity.
        changed=copy(payload);changed[length(changed)÷2]⊻=0x01;write(storage_test_payload(dir),changed)
        m=deepcopy(original);m["payload"]["sha256"]=storage_test_digest(storage_test_payload(dir));digest=storage_test_rewrite(dir,m)
        @test_throws Exception load_fixed_leptonic_context(dir;expected_manifest_sha256=digest,bindings=b,input_kind=:analytic)
        write(storage_test_payload(dir),payload)
        # Split a shared row into two equal owned objects: only alias identity changes.
        m=deepcopy(original);blocks=[o for o in m["objects"] if get(o,"name",nothing)=="HoppetPointBlock"]
        row_id=blocks[1]["fields"][2]["ref"]
        @test blocks[2]["fields"][2]["ref"]==row_id
        push!(m["objects"],deepcopy(m["objects"][row_id]));blocks[2]["fields"][2]=Dict("ref"=>length(m["objects"]))
        digest=storage_test_rewrite(dir,m)
        @test_throws ErrorException load_fixed_leptonic_context(dir;expected_manifest_sha256=digest,bindings=b,input_kind=:analytic)
        # Missing bytes are rejected before allocating reconstructed matrices.
        write(storage_test_manifest(dir),original_manifest_bytes);write(storage_test_payload(dir),payload[1:end-1])
        @test storage_test_digest(storage_test_manifest(dir))==receipt.manifest_sha256
        @test_throws ErrorException storage_test_load(dir,receipt,b)
        write(storage_test_payload(dir),payload);write(storage_test_manifest(dir),original_manifest_bytes)
        @test fixed_leptonic_storage_signature(storage_test_load(dir,receipt,b))==fixed_leptonic_storage_signature(cold)
    end
end

@testset "storage refuses unsupported or parameter-dependent metadata" begin
    b=storage_synthetic_bindings()
    mktempdir() do temporary
        for (i,key) in enumerate(("theta","gradient","HVP","replicas","optimization_state"))
            context=storage_synthetic_context(:SIDIS);context.provenance[key]=[.1,.2]
            path=joinpath(temporary,"dependent-$i")
            @test_throws ErrorException save_fixed_leptonic_context(context,path;bindings=b,input_kind=:analytic)
            @test !ispath(storage_test_manifest(path))
        end
        context=storage_synthetic_context(:SIDIS);context.provenance["cycle"]=context.provenance
        @test_throws ErrorException save_fixed_leptonic_context(context,joinpath(temporary,"cyclic");bindings=b,input_kind=:analytic)
        context=storage_synthetic_context(:SIDIS);context.first_samples.bases[1,1]=NaN
        @test_throws ErrorException save_fixed_leptonic_context(context,joinpath(temporary,"nonfinite");bindings=b,input_kind=:analytic)
        context=storage_synthetic_context(:SIDIS;metadata_extra=Dict("unsupported_float"=>Float32(.125)))
        path=joinpath(temporary,"unsupported")
        @test_throws ErrorException save_fixed_leptonic_context(context,path;bindings=b,input_kind=:analytic)
        @test !ispath(storage_test_manifest(path))
        context=storage_synthetic_context(:SIDIS)
        @test_throws ErrorException save_fixed_leptonic_context(context,joinpath(temporary,"wrong-kind");bindings=b,input_kind=:pinned_inputs)
        @test_throws ErrorException save_fixed_leptonic_context(context,joinpath(temporary,"scope");bindings=b,input_kind=:production)
    end
end
