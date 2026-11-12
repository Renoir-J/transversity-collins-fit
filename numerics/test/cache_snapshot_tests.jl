struct SyntheticSnapshotFloat <: AbstractFloat
    value::Float64
end
const synthetic_snapshot_type_label=Ref("synthetic-cache-float-v1")
Base.show(io::IO,::Type{SyntheticSnapshotFloat})=print(io,synthetic_snapshot_type_label[])

@testset "primitive type-name optimization preserves parent snapshot bytes" begin
    fixture=JSON.parsefile(joinpath(@__DIR__,"../../fixtures/closure/snapshot-digests.json"))
    @test fixture["scope"]=="synthetic_cache_integrity_fixture"
    @test fixture["contract"]=="KPSY-NUMERICAL-CACHE-001"
    @test fixture["production_ready"]===false
    @test fixture["native_endian"]==Int(Base.ENDIAN_BOM)
    primitives=(Bool,Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128,Float16,Float32,Float64)
    types=Dict(string(T)=>T for T in primitives)
    @test Set(r["type"] for r in fixture["records"])==Set(keys(types))
    for record in fixture["records"]
        T=types[record["type"]]
        @test CollinsFit.numerical_snapshot(one(T))==record["scalar_one"]
        @test CollinsFit.numerical_snapshot(reshape(T[one(T),zero(T)],1,2))==record["matrix_one_zero"]
    end
    composite=Dict("a"=>[.1,.2],"b"=>false,"c"=>(:SIDIS,2,(-0.,"synthetic")))
    @test CollinsFit.numerical_snapshot(composite)==fixture["synthetic_composite"]
    value=SyntheticSnapshotFloat(.25);array=[value]
    old_scalar=CollinsFit.numerical_snapshot(value);old_array=CollinsFit.numerical_snapshot(array)
    try
        synthetic_snapshot_type_label[]="synthetic-cache-float-v2"
        @test CollinsFit.snapshot_type_name(SyntheticSnapshotFloat)=="synthetic-cache-float-v2"
        @test CollinsFit.numerical_snapshot(value)!=old_scalar
        @test CollinsFit.numerical_snapshot(array)!=old_array
    finally
        synthetic_snapshot_type_label[]="synthetic-cache-float-v1"
    end
    @test CollinsFit.numerical_snapshot(value)==old_scalar
    @test CollinsFit.numerical_snapshot(array)==old_array
end

@testset "owned-row block type name retains parent framing" begin
    fixture=JSON.parsefile(joinpath(@__DIR__,"../../fixtures/closure/snapshot-owned-block.json"))
    blocks=[(Int[1,2],CollinsFit.fixed_rows(reshape([.25,.5,.75,1.],2,2))),
            (Int[3],CollinsFit.fixed_rows(reshape([.125,.375],2,1)))]
    @test fixture["scope"]=="synthetic_cache_integrity_fixture"
    @test fixture["native_endian"]==Int(Base.ENDIAN_BOM)
    @test fixture["native_word_size"]==Sys.WORD_SIZE
    @test CollinsFit.snapshot_type_name(eltype(blocks))==fixture["type_name"]
    @test CollinsFit.numerical_snapshot(blocks)==fixture["snapshot"]
end

@testset "complete fixed-array snapshots and typed metadata" begin
    snapshot=CollinsFit.numerical_snapshot
    values=reshape(collect(1.:100000.),100,1000)
    original=snapshot(values)
    # These include interior elements omitted by the parent general array hash.
    for index in (1,2,3,6,11,19,47,97,101,113,311,1001,1777,25001,50000,99999,100000)
        old=values[index];values[index]=old+.25
        @test snapshot(values)!=original
        values[index]=old
    end
    @test snapshot(values)==original
    @test snapshot(values)!=snapshot(vec(values))
    @test snapshot(("a","bc"))!=snapshot(("ab","c"))
    @test snapshot(false)!=snapshot(0)
    @test snapshot(1)!=snapshot(1.)
    @test snapshot(0.)!=snapshot(-0.)
    first=Dict{String,Any}("a"=>1,"b"=>false)
    reordered=Dict{String,Any}("b"=>false,"a"=>1)
    @test snapshot(first)==snapshot(reordered)
    reordered["b"]=0
    @test snapshot(first)!=snapshot(reordered)
    @test_throws ErrorException snapshot(Dict(:unsupported=>1))
    @test_throws ErrorException snapshot(Ref(1.))
end

@testset "owned read-only rows preserve numerical multiplication" begin
    raw=reshape(collect(.1:.1:2.),4,5)
    rows=CollinsFit.fixed_rows(raw)
    original=copy(rows);signature=CollinsFit.numerical_snapshot(rows)
    raw[1]+=1
    @test rows==original
    @test CollinsFit.numerical_snapshot(rows)==signature
    @test sizeof(rows)==sizeof(original)
    @test CollinsFit.fixed_rows(rows)===rows
    @test_throws ErrorException (rows[1]=0.)
    @test_throws ErrorException (rows[3,2]=0.)
    @test_throws ErrorException (view(rows,1:2,2).=0.)
    @test_throws ErrorException fill!(rows,0.)
    @test_throws ErrorException copyto!(rows,zeros(size(rows)))
    @test_throws ErrorException rows._storage
    @test_throws MethodError CollinsFit.FixedRowMatrix(original,"forged snapshot")
    copied=copy(rows);copied[1]+=1
    @test rows==original
    sliced=rows[:,1];sliced[1]+=1
    @test rows==original
    initial=reshape(collect(.1:.1:1.),5,2)
    cotangent=reshape(collect(.1:.1:.8),4,2)
    @test rows*initial==original*initial
    @test transpose(rows)*cotangent==transpose(original)*cotangent
    @test rows*BigFloat.(initial)==original*BigFloat.(initial)
    prediction(t)=sum(abs2,rows*reshape(t,5,2))
    reference(t)=sum(abs2,original*reshape(t,5,2))
    @test ForwardDiff.gradient(prediction,vec(initial))==ForwardDiff.gradient(reference,vec(initial))
    @test ForwardDiff.hessian(prediction,vec(initial))==ForwardDiff.hessian(reference,vec(initial))
end

function test_cache_mutation_rejected(f,array,index)
    before=f();old=array[index]
    try
        array[index]=old+(old isa Integer ? 1 : 1e-8)
        @test_throws ErrorException f()
    finally
        array[index]=old
    end
    @test f()==before
end

@testset "prepared operators, samples, contexts and observables reject changed fixed data" begin
    theta=collect(kpsy_published_parameters().values)
    parameters=KPSYParameters(theta)
    for process in (:SIDIS,:SIA)
        context=observable_fixture(process)
        operator=context.first_operator
        initial=polarized_initial_columns(context.first_samples,parameters)
        forward()=apply_operator(operator,initial)
        cotangent=ones(length(operator.queries),size(initial,2),length(operator.tau))
        adjoint()=operator_adjoint(operator,cotangent)
        block=operator.blocks[1]
        @test_throws ErrorException (block.rows[1,1]+=1e-8)
        for data in (block.indices,operator.grid,operator.queries,operator.Q2,operator.tau,operator.alpha)
            index=max(1,length(data)÷2)
            test_cache_mutation_rejected(forward,data,index)
            test_cache_mutation_rejected(adjoint,data,index)
        end
        samples=context.first_samples
        columns()=polarized_initial_columns(samples,parameters)
        for data in (samples.grid,samples.logarithm,samples.complement_logarithm,samples.bases)
            test_cache_mutation_rejected(columns,data,max(1,length(data)÷2))
        end
        structure()=polarized_structure(context,theta)
        plan=context.leptonic
        for data in (context.unpolarized_first,context.unpolarized_second,
                     plan.radial.b,plan.radial.db,plan.radial.edges,plan.model.nodes,plan.model.weights,
                     plan.flavors,plan.antiparticle,plan.charge2,plan.uu_weights,plan.collins_weights,
                     plan.uu_damping,plan.fixed_exponent,plan.fragmentation_shape)
            test_cache_mutation_rejected(structure,data,max(1,length(data)÷2))
        end
        momentum=MomentumRule(.05,.35;order=16)
        observable=process==:SIDIS ?
            sidis_observable(context;definition=:structure_function_ratio,momentum=momentum,source="synthetic cache fixture") :
            sia_observable(context;ratio=:UL,momentum=momentum,
                angular=SIAMeanOfRatio(.6;source="synthetic angles",averaging_measure="fixed synthetic weights"),source="synthetic cache fixture")
        prediction()=predict_observable(observable,theta)
        for data in (momentum.nodes,momentum.weights,observable.uu_weights,observable.collins_weights,observable.denominators)
            test_cache_mutation_rejected(prediction,data,max(1,length(data)÷2))
        end
        old=observable.provenance["source"]
        observable.provenance["source"]="changed synthetic source"
        @test_throws ErrorException prediction()
        observable.provenance["source"]=old
        @test isfinite(prediction())
    end
end
