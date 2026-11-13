module JetBornReference
using JSON, SHA

# Independent angular-variable forms transcribed from PDG51.4-51.12.
# The return is (s/alpha_s^2)*(d sigma/dOmega); production uses Mandelstam H_U.
function pdg_angular(kind,y)
    r=1-y
    kind==:qqprime_t && return (1+r^2)/(9y^2)
    kind==:qqprime_u && return (1+y^2)/(9r^2)
    kind==:qq_same && return ((1+r^2)/y^2+(1+y^2)/r^2-2/(3y*r))/9
    kind==:qbarq_same_t && return ((y^2+r^2)+(1+r^2)/y^2+2r^2/(3y))/9
    kind==:qbarq_same_u && return pdg_angular(:qbarq_same_t,r)
    kind==:qbarq_other && return (y^2+r^2)/9
    kind==:qbarq_gluon && return (8//27)*(y^2+r^2)*(1/(y*r)-9//4)
    kind==:qg_quark && return (1+r^2)*(1/r+9/(4y^2))/9
    kind in (:qg_gluon,:gq_quark) && return pdg_angular(:qg_quark,r)
    kind==:gq_gluon && return pdg_angular(:qg_quark,y)
    kind==:gg_quark && return (y^2+r^2)*(1/(y*r)-9//4)/24
    kind==:gg_gluon && return (9//8)*(3-y*r+r/y^2+y/r^2)
    error("unknown PDG reference channel")
end

# CHL1994 Eqs.21-26, y=(1-cos(theta))/2. These ratios are not copied
# from the production transfer-hard expressions.
function chl_transfer(kind,y)
    r=1-y
    kind==:zero && return zero(y)
    kind in (:qg,:qqprime) && return 2r/(1+r^2)
    kind==:qq_same && return r^2*(3-4y)/((3-y+y^2)*(1-3y+3y^2))
    kind==:qbarq_same && return r*(3+y)/(3-2y+y^2-2y^3+3y^4)
    kind==:qbarq_opposite && return -r^2/(3-6y+13y^2-10y^3+3y^4)
    error("unknown CHL reference channel")
end

function reference_hard(uk,tk,y)
    U=4pdg_angular(uk,y)
    (UU=U,UT=chl_transfer(tk,y)*U)
end

# Evaluate only numeric arithmetic and explicit source-index lookups.
# Never eval/include external Python or Wolfram source.
function arithmetic_expression(text,environment)
    expression=text isa AbstractString ? Meta.parse(replace(text,"**"=>"^")) : text
    function visit(value)
        value isa Integer && return BigInt(value)//BigInt(1)
        value isa Rational && return Rational{BigInt}(value)
        value isa AbstractFloat && return rationalize(BigInt,value;tol=0)
        value isa Symbol && return environment[value]
        value isa Expr || error("unsupported reference expression")
        if value.head==:ref
            length(value.args)==2 || error("unsupported reference index")
            index=value.args[2]
            index isa Integer && !(index isa Bool) || error("literal source index required")
            return visit(value.args[1])[index]
        end
        value.head==:call || error("non-arithmetic reference expression")
        operator=value.args[1]
        operator in (:+,:-,:*,:/,:^) || error("forbidden reference operation")
        arguments=map(visit,value.args[2:end])
        if operator==:^
            length(arguments)==2 && arguments[2] isa Rational && denominator(arguments[2])==1 || error("integer reference power required")
            return arguments[1]^Int(arguments[2])
        end
        operator==:+ && return +(arguments...)
        operator==:- && return -(arguments...)
        operator==:* && return *(arguments...)
        operator==:/ && return /(arguments...)
        error("unsupported reference operation")
    end
    visit(expression)
end

parse_arithmetic(text)=Meta.parse(replace(text,"**"=>"^"))

function load_source_references(path)
    root=normpath(joinpath(@__DIR__,"../.."))
    data=JSON.parsefile(path)
    data["status"]=="PREPARED_NOT_VALIDATED" && data["scope"]=="pinned_source_expressions_not_native_physics_validation" || error("incorrect source-reference scope")
    all(data[k]===false for k in ("production_ready","physical_extraction","native_FeynFacet_calculation_executed","legacy_imports_or_pickle_deserialization_executed")) || error("incorrect source-reference admission flags")
    data["contract_sha256"]==bytes2hex(sha256(read(joinpath(root,"config/star-born-contract.json")))) || error("stale source-reference contract")
    data["preparer_sha256"]==bytes2hex(sha256(read(joinpath(root,"scripts/prepare_star_born_references.py")))) || error("stale source-reference preparer")
    length(data["sources"])==3 || error("three pinned reference assets required")
    for source in data["sources"]
        bytes2hex(sha256(read(source["raw_path"])))==source["git_blob_sha256"] || error("changed source-reference bytes")
    end
    expressions=data["legacy"]["coefficients"]
    U=Dict(parse(Int,k)=>parse_arithmetic(v) for (k,v) in expressions["get_Hupol"])
    T=Dict(parse(Int,k)=>parse_arithmetic(v) for (k,v) in expressions["get_HTffa"])
    contractions=Dict()
    for (name,value) in data["legacy"]["contractions"]
        contractions[name]=(aliases=Dict(Symbol(k)=>parse_arithmetic(v) for (k,v) in value["aliases"]),
                            terms=parse_arithmetic.(value["terms"]))
    end
    fU=Dict(k=>parse_arithmetic(v) for (k,v) in data["FeynFacet"]["UU_D_dimensional"])
    fT=Dict(k=>(process=v["process"],XX=parse_arithmetic(v["XX"]),NN=parse_arithmetic(v["NN"])) for (k,v) in data["FeynFacet"]["transfer_Born"])
    (data=data,UU=U,UT=T,contractions=contractions,FeynUU=fU,FeynUT=fT)
end

function legacy_hard(references,s,t,u)
    H=Dict{Int,Any}(0=>zero(s));HT=Dict{Int,Any}()
    environment=Dict{Symbol,Any}(:s=>s,:t=>t,:u=>u,:Hupol=>H,:HTffa=>HT)
    for index in sort(collect(keys(references.UU)))
        H[index]=arithmetic_expression(references.UU[index],environment)
    end
    for index in sort(collect(keys(references.UT)))
        HT[index]=arithmetic_expression(references.UT[index],environment)
    end
    (UU=H,UT=HT)
end

"Independent pinned seven-flavor density expression; no upstream module is imported."
function legacy_density(references,s,t,u,fA,fB,D,hA,C)
    all(length(v)==7 for v in (fA,fB,D,hA,C)) || error("legacy reference has seven explicit flavors")
    hard=legacy_hard(references,s,t,u)
    indexed(v)=Dict(i-1=>v[i] for i in eachindex(v))
    fragmentation=Dict(i=>Dict(j-1=>hard.UT[i]*C[j] for j in eachindex(C)) for i in keys(hard.UT))
    environment=Dict{Symbol,Any}(:f=>indexed(fA),:ft=>indexed(fB),:d=>indexed(D),:h=>indexed(hA),
                                :Hupol=>hard.UU,:Hxxpz=>fragmentation)
    values=map(("_get_dsig","_get_dsigST")) do name
        local_environment=copy(environment)
        expression=references.contractions[name]
        for (alias,rhs) in expression.aliases
            local_environment[alias]=arithmetic_expression(rhs,local_environment)
        end
        sum(arithmetic_expression(term,local_environment) for term in expression.terms)
    end
    (UU=values[1],UT=values[2])
end

function tagged_ids(name)
    left,right=split(name," -> ")
    first,second=split(left)
    labels=Dict("q"=>2,"qb"=>-2,"qp"=>1,"qpb"=>-1,"g"=>21)
    (labels[first],labels[second],labels[right])
end

"Relabel released reference channels by flavor permutations and charge conjugation."
function exported_channel_table(references,pids,s,t,u)
    U=Dict{NTuple{3,Int},Any}();T=Dict{NTuple{3,Int},Any}()
    environment=Dict{Symbol,Any}(:s=>s,:t=>t,:u=>u,:D=>4//1,:CA=>3//1,:CF=>4//3,
        :Nc=>3//1,:v=>1+t/s,:alphaS=>1//5)
    flavors=filter(>(0),filter(!=(21),collect(pids)))
    function add!(table,key,value)
        haskey(table,key) && table[key]!=value && error("inconsistent relabeled source channels")
        table[key]=value
    end
    for a in flavors,b in flavors
        a==b && continue
        for charge in (-1,1)
            labels=Dict("q"=>charge*a,"qb"=>-charge*a,"qp"=>charge*b,"qpb"=>-charge*b,"g"=>21)
            relabel(name)=begin
                left,right=split(name," -> ");first,second=split(left)
                (labels[first],labels[second],labels[right])
            end
            for (name,expression) in references.FeynUU
                add!(U,relabel(name),arithmetic_expression(expression,environment))
            end
            for row in values(references.FeynUT)
                xx=arithmetic_expression(row.XX,environment)*s*(1+t/s)/(environment[:alphaS]^2)
                nn=arithmetic_expression(row.NN,environment)*s*(1+t/s)/(environment[:alphaS]^2)
                xx==nn || error("released Born XX/NN disagree")
                add!(T,relabel(row.process),xx)
            end
        end
    end
    (UU=U,UT=T)
end

end
