const FIXED_STORAGE_CONTRACT=normpath(joinpath(@__DIR__,"../../config/fixed-leptonic-storage-contract.json"))
const FIXED_STORAGE_CONTRACT_SHA="cf7de20446d7b80c57559de54d807ff8b73729240cfa0c3f6bbde464ed1c0996"
const FIXED_STORAGE_COMPONENTS=Dict(string(nameof(T))=>T for T in
    (AsymptoticTwoLoopCoupling,KPSYCSS,RadialRule,LeptonicPlan,NucleonTarget,
     PolarizedInitialSamples,HoppetMatchingPlan,HoppetPointBlock,FixedRowMatrix,PolarizedLeptonicPlan))
const FIXED_STORAGE_PRIMITIVES=Dict("Float64"=>Float64,"Int64"=>Int64,"UInt64"=>UInt64,"Bool"=>Bool)
const FIXED_STORAGE_FALSE_FLAGS=("production_ready","physical_extraction","historical_fit_ready",
    "experimental_data_admitted","fit_performed","accepted_prediction_diagnostic",
    "accepted_fitting_gate","G2PERF_pass","interval_coverage_verified")

fixed_storage_digest(path)=open(sha256,path) |> bytes2hex
fixed_storage_sha(value)=value isa AbstractString && occursin(r"^[0-9a-f]{64}$",value)

function fixed_context_storage_contract()
    fixed_storage_digest(FIXED_STORAGE_CONTRACT)==FIXED_STORAGE_CONTRACT_SHA || error("changed frozen storage contract")
    c=JSON.parsefile(FIXED_STORAGE_CONTRACT)
    c["id"]=="KPSY-FIXED-LEPTONIC-STORAGE-001" && c["scope"]=="fixed_leptonic_context_storage_integrity" || error("wrong storage contract")
    all(c[k]===false for k in FIXED_STORAGE_FALSE_FLAGS) && c["equation_change"]===false || error("storage contract cannot self-admit")
    Sys.WORD_SIZE==64 || error("storage supports the declared64-bit layout only")
    for (path,h) in c["parent_sha256"]
        fixed_storage_digest(normpath(joinpath(@__DIR__,"../..",path)))==h || error("changed storage parent contract")
    end
    c
end

function fixed_storage_bindings(bindings,c)
    bindings isa AbstractDict && Set(keys(bindings))==Set(c["required_binding_fields"]) || error("exact storage bindings required")
    all(k isa String && fixed_storage_sha(v) for (k,v) in bindings) || error("storage bindings require SHA256 values")
    bindings["storage_contract_sha256"]==FIXED_STORAGE_CONTRACT_SHA || error("stale storage binding")
    Dict{String,String}(bindings)
end

function fixed_storage_input_kind(kind)
    kind in (:analytic,:pinned_inputs) || error("explicit analytic or pinned_inputs storage kind required")
    kind
end

function fixed_storage_keys(value,keys)
    value isa AbstractDict && Set(Base.keys(value))==Set(keys) || error("invalid fixed storage object fields")
end

function fixed_storage_int(value;minimum=0,maximum=typemax(Int))
    value isa Integer && !(value isa Bool) && minimum<=value<=maximum || error("invalid bounded storage integer")
    Int(value)
end

function fixed_storage_type(T,depth=0)
    depth<=64 || error("storage type nesting exhausted")
    T isa Type || error("unsupported storage type schema")
    for (name,type) in FIXED_STORAGE_PRIMITIVES
        T===type && return name
    end
    T===Any && return "Any"
    T===String && return "String"
    T===Symbol && return "Symbol"
    T===Nothing && return "Nothing"
    for (name,type) in FIXED_STORAGE_COMPONENTS
        T===type && return "component:"*name
    end
    if T<:Dict && keytype(T)===String
        return Dict("kind"=>"dict_type","value"=>fixed_storage_type(valtype(T),depth+1))
    elseif T<:Tuple && isconcretetype(T)
        return Dict("kind"=>"tuple_type","items"=>[fixed_storage_type(t,depth+1) for t in T.parameters])
    elseif T<:Array && ndims(T) in (1,2)
        return Dict("kind"=>"array_type","element"=>fixed_storage_type(eltype(T),depth+1),"dimensions"=>ndims(T))
    end
    error("unsupported fixed storage element type $T")
end

function fixed_storage_type(schema::Union{AbstractString,AbstractDict},depth=0)
    depth<=64 || error("storage type nesting exhausted")
    if schema isa AbstractString
        haskey(FIXED_STORAGE_PRIMITIVES,schema) && return FIXED_STORAGE_PRIMITIVES[schema]
        schema=="Any" && return Any
        schema=="String" && return String
        schema=="Symbol" && return Symbol
        schema=="Nothing" && return Nothing
        if startswith(schema,"component:")
            name=schema[11:end]
            haskey(FIXED_STORAGE_COMPONENTS,name) || error("unsupported component type")
            return FIXED_STORAGE_COMPONENTS[name]
        end
        error("unsupported storage type label")
    end
    kind=get(schema,"kind",nothing)
    if kind=="dict_type"
        fixed_storage_keys(schema,("kind","value"))
        return Core.apply_type(Dict,String,fixed_storage_type(schema["value"],depth+1))
    elseif kind=="tuple_type"
        fixed_storage_keys(schema,("kind","items"))
        schema["items"] isa AbstractVector && length(schema["items"])<=100000 || error("invalid tuple type")
        return Core.apply_type(Tuple,(fixed_storage_type(s,depth+1) for s in schema["items"])...)
    elseif kind=="array_type"
        fixed_storage_keys(schema,("kind","element","dimensions"))
        n=fixed_storage_int(schema["dimensions"];minimum=1,maximum=2)
        return Core.apply_type(Array,fixed_storage_type(schema["element"],depth+1),n)
    end
    error("unsupported storage type schema")
end

mutable struct FixedStorageEncoder
    io::IOStream
    objects::Vector{Any}
    ids::IdDict{Any,Int}
    active::IdDict{Any,Bool}
    policy::Dict{String,Any}
end

function fixed_storage_node!(writer,value,depth=0)
    depth<=writer.policy["maximum_nesting_depth"] || error("fixed storage nesting exhausted")
    value===nothing && return Dict("kind"=>"nothing")
    if typeof(value) in values(FIXED_STORAGE_PRIMITIVES)
        value isa AbstractFloat && !isfinite(value) && error("nonfinite fixed storage scalar")
        return Dict("kind"=>"scalar","type"=>string(typeof(value)),"bits"=>bitstring(value))
    elseif value isa String || value isa Symbol
        return Dict("kind"=>value isa String ? "string" : "symbol","value"=>string(value))
    elseif value isa Tuple
        return Dict("kind"=>"tuple","items"=>[fixed_storage_node!(writer,x,depth+1) for x in value])
    end
    haskey(writer.active,value) && error("cyclic fixed storage metadata")
    haskey(writer.ids,value) && return Dict("ref"=>writer.ids[value])
    id=length(writer.objects)+1
    id<=writer.policy["maximum_objects"] || error("too many fixed storage objects")
    writer.ids[value]=id;writer.active[value]=true;push!(writer.objects,nothing)
    object=if value isa Array && eltype(value) in values(FIXED_STORAGE_PRIMITIVES)
        1<=ndims(value)<=writer.policy["maximum_array_dimensions"] || error("unsupported array dimensions")
        eltype(value)===Float64 && !all(isfinite,value) && error("nonfinite fixed storage array")
        offset=position(writer.io);nbytes=sizeof(value)
        offset+nbytes<=writer.policy["maximum_payload_bytes"] || error("fixed storage payload budget exhausted")
        write(writer.io,value)
        Dict("kind"=>"array","type"=>string(eltype(value)),"shape"=>collect(size(value)),"offset"=>offset,"bytes"=>nbytes)
    elseif value isa Vector
        Dict("kind"=>"vector","element_type"=>fixed_storage_type(eltype(value)),
            "items"=>[fixed_storage_node!(writer,x,depth+1) for x in value])
    elseif value isa Dict && keytype(typeof(value))===String
        names=sort(collect(keys(value)))
        Dict("kind"=>"dict","value_type"=>fixed_storage_type(valtype(typeof(value))),"keys"=>names,
            "values"=>[fixed_storage_node!(writer,value[k],depth+1) for k in names])
    elseif value isa FixedRowMatrix
        raw=getfield(value,:_storage)
        numerical_snapshot(raw)==value.snapshot || error("mutated owned rows cannot be stored")
        Dict("kind"=>"fixedrows","array"=>fixed_storage_node!(writer,raw,depth+1),"snapshot"=>value.snapshot)
    elseif haskey(FIXED_STORAGE_COMPONENTS,string(nameof(typeof(value)))) && FIXED_STORAGE_COMPONENTS[string(nameof(typeof(value)))]===typeof(value)
        names=collect(String.(fieldnames(typeof(value))))
        Dict("kind"=>"component","name"=>string(nameof(typeof(value))),"field_names"=>names,
            "fields"=>[fixed_storage_node!(writer,getfield(value,f),depth+1) for f in fieldnames(typeof(value))])
    else
        error("unsupported fixed storage value $(typeof(value))")
    end
    writer.objects[id]=object;delete!(writer.active,value)
    Dict("ref"=>id)
end

function fixed_storage_context_valid(context,kind)
    context isa PolarizedLeptonicPlan || error("prepared leptonic context required")
    fixed_storage_input_kind(kind)
    all(get(context.provenance,k,false)===false for k in FIXED_STORAGE_FALSE_FLAGS) || error("stored context cannot claim admission")
    function fixed_metadata(value,depth=0)
        depth<=64 || error("fixed metadata nesting exhausted")
        if value isa AbstractDict
            forbidden=("theta","fitted_parameters","parameter_vector","parameter_response","prediction","predictions",
                "gradient","gradients","hvp","hessian_vector_product","optimization_state","replica","replicas","posterior")
            any(lowercase(k) in forbidden for k in keys(value)) && error("parameter-dependent storage metadata is unsupported")
            for v in values(value);fixed_metadata(v,depth+1);end
        elseif value isa Union{AbstractVector,Tuple}
            for v in value;fixed_metadata(v,depth+1);end
        end
    end
    fixed_metadata(context.provenance)
    polarized_context_valid(context)
    for (operator,samples) in ((context.first_operator,context.first_samples),(context.second_operator,context.second_samples))
        operator===nothing && samples===nothing && continue
        operator isa HoppetMatchingPlan && samples isa PolarizedInitialSamples || error("incomplete fixed leg")
        samples.input_kind===kind || error("wrong stored initial input kind")
        hoppet_point_signature(operator)==operator.signature || error("mutated detached matching operator")
        numerical_snapshot((samples.kind,samples.grid,samples.logarithm,samples.complement_logarithm,samples.bases,samples.input_kind))==samples.signature || error("mutated fixed initial samples")
    end
    context
end

function fixed_storage_hash!(sha,value,aliases,events,depth=0)
    depth<=64 || error("fixed storage logical nesting exhausted")
    if value isa FixedRowMatrix
        raw=getfield(value,:_storage)
        numerical_snapshot(raw)==value.snapshot || error("mutated owned rows")
        id=get!(aliases,value,length(aliases)+1);push!(events,id)
        snapshot_tag!(sha,"FIXED-ROWS:");snapshot_update!(sha,value.snapshot)
        snapshot_update!(sha,raw)
    elseif value isa AbstractArray
        snapshot_tag!(sha,"ARRAY:");snapshot_update!(sha,string(eltype(value)));snapshot_update!(sha,size(value))
        if value isa Array && eltype(value) in values(FIXED_STORAGE_PRIMITIVES)
            eltype(value)===Float64 && !all(isfinite,value) && error("nonfinite fixed logical array")
            SHA.update!(sha,reinterpret(UInt8,vec(value)))
        else
            for x in value;fixed_storage_hash!(sha,x,aliases,events,depth+1);end
        end
    elseif value isa Tuple
        snapshot_tag!(sha,"TUPLE:");snapshot_update!(sha,length(value))
        for x in value;fixed_storage_hash!(sha,x,aliases,events,depth+1);end
    elseif value isa AbstractDict
        all(k isa String for k in keys(value)) || error("string metadata keys required")
        snapshot_tag!(sha,"DICT:");snapshot_update!(sha,length(value))
        for k in sort(collect(keys(value)))
            snapshot_update!(sha,k);fixed_storage_hash!(sha,value[k],aliases,events,depth+1)
        end
    elseif haskey(FIXED_STORAGE_COMPONENTS,string(nameof(typeof(value)))) && FIXED_STORAGE_COMPONENTS[string(nameof(typeof(value)))]===typeof(value)
        snapshot_tag!(sha,"COMPONENT:");snapshot_update!(sha,string(nameof(typeof(value))))
        for f in fieldnames(typeof(value))
            snapshot_update!(sha,string(f));fixed_storage_hash!(sha,getfield(value,f),aliases,events,depth+1)
        end
    else
        value isa AbstractFloat && !isfinite(value) && error("nonfinite fixed logical scalar")
        snapshot_update!(sha,value)
    end
end

"Complete fixed-data identity and the canonical repeated-row reference sequence."
function fixed_leptonic_storage_signature(context::PolarizedLeptonicPlan)
    sha=SHA.SHA2_256_CTX();snapshot_tag!(sha,"KPSY-FIXED-CONTEXT-LOGICAL-001:")
    events=Int[];fixed_storage_hash!(sha,context,IdDict{Any,Int}(),events)
    (;logical_sha256=bytes2hex(SHA.digest!(sha)),row_alias_sha256=numerical_snapshot(events))
end

function fixed_storage_array_layout(object,policy)
    fixed_storage_keys(object,("kind","type","shape","offset","bytes"))
    haskey(FIXED_STORAGE_PRIMITIVES,object["type"]) || error("unsupported binary element type")
    T=FIXED_STORAGE_PRIMITIVES[object["type"]]
    shape=object["shape"]
    shape isa AbstractVector && 1<=length(shape)<=policy["maximum_array_dimensions"] || error("invalid binary dimensions")
    dims=Tuple(fixed_storage_int(d;maximum=policy["maximum_payload_bytes"]) for d in shape)
    count=foldl(Base.checked_mul,dims;init=1);nbytes=Base.checked_mul(count,sizeof(T))
    nbytes<=policy["maximum_payload_bytes"] || error("binary array allocation budget exhausted")
    offset=fixed_storage_int(object["offset"];maximum=policy["maximum_payload_bytes"])
    declared=fixed_storage_int(object["bytes"];maximum=policy["maximum_payload_bytes"])
    nbytes==declared || error("binary shape/size mismatch")
    (;T,dims,offset,nbytes)
end

mutable struct FixedStorageDecoder
    io::IOStream
    objects::Vector{Any}
    cache::Dict{Int,Any}
    active::Set{Int}
    policy::Dict{String,Any}
end

function fixed_storage_value(reader,node,depth=0)
    depth<=reader.policy["maximum_nesting_depth"] || error("fixed storage nesting exhausted")
    node isa AbstractDict || error("typed storage node required")
    if haskey(node,"ref")
        fixed_storage_keys(node,("ref",))
        id=fixed_storage_int(node["ref"];minimum=1,maximum=length(reader.objects))
        haskey(reader.cache,id) && return reader.cache[id]
        id in reader.active && error("cyclic stored object reference")
        push!(reader.active,id);object=reader.objects[id]
        object isa AbstractDict || error("missing fixed storage object")
        kind=get(object,"kind",nothing)
        value=if kind=="array"
            layout=fixed_storage_array_layout(object,reader.policy)
            a=Array{layout.T}(undef,layout.dims);seek(reader.io,layout.offset);read!(reader.io,a)
            layout.T===Float64 && !all(isfinite,a) && error("nonfinite stored array")
            layout.T===Bool && !all(x->x in (0x00,0x01),reinterpret(UInt8,vec(a))) && error("invalid stored Bool bytes")
            a
        elseif kind=="vector"
            fixed_storage_keys(object,("kind","element_type","items"))
            T=fixed_storage_type(object["element_type"])
            object["items"] isa AbstractVector || error("stored vector items required")
            values=[fixed_storage_value(reader,n,depth+1) for n in object["items"]]
            all(v isa T for v in values) || error("stored vector element type mismatch")
            Vector{T}(values)
        elseif kind=="dict"
            fixed_storage_keys(object,("kind","value_type","keys","values"))
            T=fixed_storage_type(object["value_type"]);names=object["keys"];nodes=object["values"]
            names isa AbstractVector && nodes isa AbstractVector && length(names)==length(nodes) || error("invalid stored dictionary")
            all(k isa String for k in names) && issorted(names) && length(unique(names))==length(names) || error("invalid stored dictionary keys")
            d=Dict{String,T}()
            for (k,n) in zip(names,nodes)
                v=fixed_storage_value(reader,n,depth+1);v isa T || error("stored metadata value type mismatch");d[k]=v
            end
            d
        elseif kind=="fixedrows"
            fixed_storage_keys(object,("kind","array","snapshot"))
            raw=fixed_storage_value(reader,object["array"],depth+1)
            raw isa Matrix{Float64} || error("fixed rows require a Float64 matrix")
            rows=FixedRowMatrix(raw)
            rows.snapshot==object["snapshot"] || error("stored row snapshot mismatch")
            rows
        elseif kind=="component"
            fixed_storage_keys(object,("kind","name","field_names","fields"))
            name=object["name"];haskey(FIXED_STORAGE_COMPONENTS,name) && name!="FixedRowMatrix" || error("unsupported stored component")
            T=FIXED_STORAGE_COMPONENTS[name]
            object["field_names"]==collect(String.(fieldnames(T))) && length(object["fields"])==fieldcount(T) || error("stored component schema mismatch")
            values=[fixed_storage_value(reader,n,depth+1) for n in object["fields"]]
            all(values[i] isa fieldtype(T,i) for i in eachindex(values)) || error("stored component field type mismatch")
            T(values...)
        else
            error("unsupported stored object kind")
        end
        reader.cache[id]=value;delete!(reader.active,id)
        return value
    end
    kind=get(node,"kind",nothing)
    if kind=="nothing"
        fixed_storage_keys(node,("kind",));return nothing
    elseif kind in ("string","symbol")
        fixed_storage_keys(node,("kind","value"));node["value"] isa String || error("stored text required")
        return kind=="string" ? node["value"] : Symbol(node["value"])
    elseif kind=="tuple"
        fixed_storage_keys(node,("kind","items"));node["items"] isa AbstractVector || error("stored tuple items required")
        return Tuple(fixed_storage_value(reader,n,depth+1) for n in node["items"])
    elseif kind=="scalar"
        fixed_storage_keys(node,("kind","type","bits"))
        haskey(FIXED_STORAGE_PRIMITIVES,node["type"]) || error("unsupported stored scalar")
        T=FIXED_STORAGE_PRIMITIVES[node["type"]];bits=node["bits"]
        bits isa String && ncodeunits(bits)==8sizeof(T) && all(c->c in ('0','1'),bits) || error("invalid scalar bits")
        value=if T===Bool
            bits in ("00000000","00000001") || error("invalid Bool scalar bits");bits=="00000001"
        else
            u=parse(UInt64,bits;base=2);T===UInt64 ? u : reinterpret(T,u)
        end
        value isa AbstractFloat && !isfinite(value) && error("nonfinite stored scalar")
        return value
    end
    error("unsupported typed storage node")
end

"Persist only fixed prepared data; a returned file digest is not a physical gate."
function save_fixed_leptonic_context(context::PolarizedLeptonicPlan,destination;
        bindings,input_kind)
    c=fixed_context_storage_contract();kind=fixed_storage_input_kind(input_kind)
    bound=fixed_storage_bindings(bindings,c);policy=c["format"]
    fixed_storage_context_valid(context,kind)
    ispath(destination) && error("preserve prior storage attempts; use a fresh destination")
    mkpath(dirname(abspath(destination)));mkdir(destination)
    payload=joinpath(destination,policy["payload_filename"])
    writer=nothing;root=nothing
    open(payload,"w") do io
        writer=FixedStorageEncoder(io,Any[],IdDict{Any,Int}(),IdDict{Any,Bool}(),policy)
        root=fixed_storage_node!(writer,context)
        flush(io)
    end
    fixed_storage_context_valid(context,kind)
    signature=fixed_leptonic_storage_signature(context)
    manifest=Dict{String,Any}("schema_version"=>1,"format_id"=>policy["id"],"scope"=>c["scope"],
        "contract_sha256"=>FIXED_STORAGE_CONTRACT_SHA,"input_kind"=>string(kind),"bindings"=>bound,
        "producer_julia_version"=>string(VERSION),"native_endian"=>Int(Base.ENDIAN_BOM),"word_size"=>Sys.WORD_SIZE,
        "payload"=>Dict("file"=>policy["payload_filename"],"bytes"=>filesize(payload),"sha256"=>fixed_storage_digest(payload)),
        "root"=>root,"objects"=>writer.objects,"logical_context_sha256"=>signature.logical_sha256,
        "row_alias_sha256"=>signature.row_alias_sha256)
    for flag in FIXED_STORAGE_FALSE_FLAGS;manifest[flag]=false;end
    text=JSON.json(manifest)
    ncodeunits(text)<=policy["maximum_manifest_bytes"] || error("fixed manifest budget exhausted")
    # Decode through the same strict reader before publishing completion.
    fixed_storage_decode(manifest,payload,c,bound,kind)
    # The completion manifest is published only after full payload verification.
    open(joinpath(destination,policy["manifest_filename"]*".writing"),"w") do io;write(io,text);end
    manifest_path=joinpath(destination,policy["manifest_filename"])
    mv(manifest_path*".writing",manifest_path;force=false)
    (;manifest_sha256=fixed_storage_digest(manifest_path),payload_sha256=manifest["payload"]["sha256"],
      logical_context_sha256=signature.logical_sha256,row_alias_sha256=signature.row_alias_sha256,
      payload_bytes=filesize(payload),objects=length(writer.objects),input_kind=kind)
end

"Load a byte-bound owned context; never infer trust from an internal completion marker."
function load_fixed_leptonic_context(destination;expected_manifest_sha256,bindings,input_kind)
    c=fixed_context_storage_contract();kind=fixed_storage_input_kind(input_kind)
    bound=fixed_storage_bindings(bindings,c);policy=c["format"]
    fixed_storage_sha(expected_manifest_sha256) || error("external manifest SHA256 required")
    path=joinpath(destination,policy["manifest_filename"]);payload=joinpath(destination,policy["payload_filename"])
    isfile(path) && isfile(payload) && !islink(path) && !islink(payload) || error("complete owned storage files required")
    filesize(path)<=policy["maximum_manifest_bytes"] || error("fixed manifest budget exhausted")
    fixed_storage_digest(path)==expected_manifest_sha256 || error("external manifest digest mismatch")
    m=JSON.parsefile(path)
    context=fixed_storage_decode(m,payload,c,bound,kind)
    fixed_storage_digest(path)==expected_manifest_sha256 || error("manifest changed during loading")
    context
end

function fixed_storage_decode(m,payload,c,bound,kind)
    policy=c["format"]
    fields=("schema_version","format_id","scope","contract_sha256","input_kind","bindings","producer_julia_version",
        "native_endian","word_size","payload","root","objects","logical_context_sha256","row_alias_sha256",FIXED_STORAGE_FALSE_FLAGS...)
    fixed_storage_keys(m,fields)
    m["schema_version"]===1 && m["format_id"]==policy["id"] && m["scope"]==c["scope"] &&
        m["contract_sha256"]==FIXED_STORAGE_CONTRACT_SHA && all(m[f]===false for f in FIXED_STORAGE_FALSE_FLAGS) || error("wrong stored scope or admission claim")
    m["input_kind"]==string(kind) && fixed_storage_bindings(m["bindings"],c)==bound || error("wrong stored bindings or initial kind")
    m["producer_julia_version"]==string(VERSION) || error("unsupported stored Julia version")
    fixed_storage_int(m["native_endian"];minimum=Int(Base.ENDIAN_BOM),maximum=Int(Base.ENDIAN_BOM))
    fixed_storage_int(m["word_size"];minimum=Sys.WORD_SIZE,maximum=Sys.WORD_SIZE)
    fixed_storage_keys(m["payload"],("file","bytes","sha256"))
    m["payload"]["file"]==policy["payload_filename"] || error("wrong payload filename")
    bytes=fixed_storage_int(m["payload"]["bytes"];maximum=policy["maximum_payload_bytes"])
    filesize(payload)==bytes && fixed_storage_digest(payload)==m["payload"]["sha256"] || error("changed or truncated storage payload")
    objects=m["objects"]
    objects isa Vector && 1<=length(objects)<=policy["maximum_objects"] || error("invalid stored object count")
    spans=Tuple{Int,Int}[]
    for object in objects
        object isa AbstractDict || error("missing stored object")
        if get(object,"kind",nothing)=="array"
            layout=fixed_storage_array_layout(object,policy);push!(spans,(layout.offset,layout.nbytes))
        end
    end
    cursor=0
    for (offset,nbytes) in sort(spans)
        offset==cursor || error("overlapping or missing binary array spans")
        cursor=Base.checked_add(cursor,nbytes);cursor<=bytes || error("array span exceeds payload")
    end
    cursor==bytes || error("unclaimed payload bytes")
    context=open(payload,"r") do io
        reader=FixedStorageDecoder(io,objects,Dict{Int,Any}(),Set{Int}(),policy)
        value=fixed_storage_value(reader,m["root"])
        length(reader.cache)==length(objects) || error("unreachable stored objects")
        value
    end
    fixed_storage_context_valid(context,kind)
    signature=fixed_leptonic_storage_signature(context)
    signature.logical_sha256==m["logical_context_sha256"] && signature.row_alias_sha256==m["row_alias_sha256"] || error("stored fixed-data or row-alias identity mismatch")
    fixed_storage_digest(payload)==m["payload"]["sha256"] || error("payload changed during loading")
    context
end
