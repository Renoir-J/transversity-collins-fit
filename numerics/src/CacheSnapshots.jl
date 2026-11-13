"""A complete, typed snapshot of fixed cache data; never fitted parameters.

Every dense primitive array byte is streamed to SHA256. Julia's general array
hash may sample elements and is not used for these integrity checks. This is an
in-memory snapshot format, not an executable serialization or a production
admission certificate.
"""
function numerical_snapshot(value)
    context=SHA.SHA2_256_CTX()
    SHA.update!(context,codeunits("KPSY-NUMERICAL-CACHE-001/native-endian=$(Base.ENDIAN_BOM)/"))
    snapshot_update!(context,value)
    bytes2hex(SHA.digest!(context))
end

snapshot_tag!(context,tag)=SHA.update!(context,codeunits(tag))
snapshot_update!(context,::Nothing)=snapshot_tag!(context,"N")

snapshot_type_name(::Type{T}) where T=string(T)
for numeric_type in (Bool,Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128,
                  Float16,Float32,Float64)
    @eval snapshot_type_name(::Type{$numeric_type})=$(string(numeric_type))
end

function snapshot_update!(context,value::AbstractString)
    snapshot_tag!(context,"S$(ncodeunits(value)):")
    SHA.update!(context,codeunits(value))
end
function snapshot_update!(context,value::Symbol)
    snapshot_tag!(context,"Y")
    snapshot_update!(context,string(value))
end
function snapshot_update!(context,value::Union{Integer,AbstractFloat})
    isbitstype(typeof(value)) || error("fixed snapshot requires a primitive numeric scalar")
    snapshot_tag!(context,"P")
    snapshot_update!(context,snapshot_type_name(typeof(value)))
    SHA.update!(context,reinterpret(UInt8,[value]))
end
function snapshot_update!(context,value::Tuple)
    snapshot_tag!(context,"T$(length(value)):")
    for item in value
        snapshot_update!(context,item)
    end
end
function snapshot_update!(context,value::NamedTuple)
    snapshot_tag!(context,"K")
    snapshot_update!(context,keys(value))
    snapshot_update!(context,values(value))
end
function snapshot_update!(context,value::AbstractDict)
    all(key->key isa String,keys(value)) || error("fixed snapshot metadata requires string keys")
    snapshot_tag!(context,"D$(length(value)):")
    for key in sort!(collect(keys(value)))
        snapshot_update!(context,key)
        snapshot_update!(context,value[key])
    end
end
function snapshot_update!(context,value::AbstractArray{T}) where T
    snapshot_tag!(context,"A")
    snapshot_update!(context,snapshot_type_name(T))
    snapshot_update!(context,size(value))
    if value isa Array && isbitstype(T) && T<:Union{Integer,AbstractFloat}
        SHA.update!(context,reinterpret(UInt8,vec(value)))
    else
        for item in value
            snapshot_update!(context,item)
        end
    end
end
snapshot_update!(context,value)=error("unsupported fixed snapshot type $(typeof(value))")

"Owned, read-only retained rows; ordinary indexing never exposes mutable storage."
struct FixedRowMatrix <: AbstractMatrix{Float64}
    _storage::Matrix{Float64}
    snapshot::String
    function FixedRowMatrix(rows::AbstractMatrix{Float64})
        storage=rows isa Matrix{Float64} ? copy(rows) : Matrix{Float64}(rows)
        new(storage,numerical_snapshot(storage))
    end
end
@eval snapshot_type_name(::Type{Tuple{Vector{Int},FixedRowMatrix}})=
    $(string(Tuple{Vector{Int},FixedRowMatrix}))
fixed_rows(rows::FixedRowMatrix)=rows
fixed_rows(rows::AbstractMatrix{Float64})=FixedRowMatrix(rows)
Base.size(rows::FixedRowMatrix)=size(getfield(rows,:_storage))
Base.axes(rows::FixedRowMatrix)=axes(getfield(rows,:_storage))
Base.IndexStyle(::Type{FixedRowMatrix})=IndexLinear()
Base.getindex(rows::FixedRowMatrix,indices...)=getindex(getfield(rows,:_storage),indices...)
Base.setindex!(::FixedRowMatrix,value,indices...)=error("fixed numerical cache rows are read-only")
Base.copy(rows::FixedRowMatrix)=copy(getfield(rows,:_storage))
Base.sizeof(rows::FixedRowMatrix)=sizeof(getfield(rows,:_storage))
Base.similar(::FixedRowMatrix,::Type{T},dims::Dims) where T=Array{T}(undef,dims)
function Base.getproperty(rows::FixedRowMatrix,name::Symbol)
    name===:_storage && error("fixed row storage is private; use copy for a mutable independent array")
    getfield(rows,name)
end
Base.propertynames(::FixedRowMatrix,private::Bool=false)=(:snapshot,)
Base.:*(rows::FixedRowMatrix,initial::AbstractMatrix)=getfield(rows,:_storage)*initial
Base.:*(rows::Transpose{Float64,FixedRowMatrix},initial::AbstractMatrix)=
    transpose(getfield(parent(rows),:_storage))*initial
function snapshot_update!(context,rows::FixedRowMatrix)
    snapshot_tag!(context,"F")
    snapshot_update!(context,rows.snapshot)
end
