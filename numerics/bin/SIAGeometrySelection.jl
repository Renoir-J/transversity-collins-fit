"Arithmetic reducers for a complete72 ledger; no process dispatch or fitting gate."
module SIAGeometrySelection
using JSON, SHA
export packed_point, collect_level, select_momentum, select_native, select_tails

const PAIRS=("pp","pm","mp","mm")
const FLAGS=("production_ready","physical_extraction","historical_fit_ready","experimental_data_admitted",
    "fit_performed","G2PERF_pass","interval_coverage_verified")
const VALUES_STATUS="COMPLETE_GEOMETRY_VALUES_NOT_DIAGNOSTIC_GATE"
const SCOPE="historical_reconstructed_sia_complete_geometry_jobs"
false_flags()=Dict{String,Any}(key=>false for key in (FLAGS...,"accepted_fitting_gate"))
finite_number(x)=x isa Real && !(x isa Bool) && isfinite(x)
require(value,message)=value || error(message)
scaled(a,b,floor)=maximum(abs.(a.-b)./max.(abs.(b),floor))
packed_point(point)=vcat([point["pairs"][pair][key] for pair in PAIRS for key in ("UU","Collins")],
    [point["categories"][category][key] for category in ("U","L","C") for key in ("UU","Collins")],
    [point["amplitude"]])
packed(points)=reduce(vcat,packed_point.(points))

function same_json(a,b)
    if b isa Bool
        return a===b
    elseif b isa Integer
        return a isa Integer && !(a isa Bool) && a==b
    elseif b isa AbstractFloat
        return a isa AbstractFloat && a==b
    elseif b isa AbstractDict
        return a isa AbstractDict && Set(keys(a))==Set(keys(b)) && all(same_json(a[key],value) for (key,value) in b)
    elseif b isa AbstractVector
        return a isa AbstractVector && length(a)==length(b) && all(same_json(x,y) for (x,y) in zip(a,b))
    end
    typeof(a)==typeof(b) && a==b
end

function validate_parent(parent)
    path=joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json")
    require(bytes2hex(sha256(read(path)))=="29830e78c4e48d2eb8e5a70711885fc961b2b94e0ec3303315befd59a65fdf97","frozen scientific contract changed")
    require(same_json(parent,JSON.parsefile(path)),"caller cannot change scientific thresholds or assumptions")
end

function declaration(request)
    if haskey(request,"fixture_scope")
        require(request["fixture_scope"]=="explicitly_synthetic_arithmetic_fixture_not_primary_data" &&
            get(request,"theory_scope",nothing)=="synthetic_closure","synthetic arithmetic must retain its explicit scope")
        return "synthetic_closure"
    end
    require(get(request,"theory_scope",nothing)=="reconstructed_theory","wrong theory declaration")
    "reconstructed_theory"
end

function geometry_groups(request)
    declaration(request)
    require(get(request,"scope",nothing)=="historical_reconstructed_sia_integral_prediction_diagnostic" &&
        get(request,"source_fingerprint",nothing)=="01be55cc1449679f6dc26bd43ea9788604938c1fecbb25d276acaa35bc4bae47","wrong frozen request scope or identity")
    require(all(get(request,key,nothing)===false for key in FLAGS),"source request cannot self-admit")
    points=request["points"]
    ids=["KPSY15_SIA_$(view):r$row" for view in (1002,1003) for row in 2:37]
    require(length(points)==72 && [point["id"] for point in points]==vec(permutedims(ids)),"full ordered72 request required")
    groups=Vector{Vector{Int}}();keys=Tuple[]
    for (index,point) in enumerate(points)
        require(point["ratio"]==(index<=36 ? "UC" : "UL"),"changed source ratio")
        key=Tuple(point["kinematics"][name] for name in ("z1","z2","Q2_GeV2"))
        require(all(finite_number,key) && 0<key[1]<1 && 0<key[2]<1 && key[3]==112.36,"changed literal source geometry")
        position=findfirst(==(key),keys)
        if position===nothing
            push!(keys,key);push!(groups,Int[]);position=length(keys)
        end
        push!(groups[position],index)
    end
    require(length(keys)==36 && all(length(group)==2 && group[1]<=36<group[2] for group in groups),"exact36 UC/UL geometries required")
    groups,keys
end

function validate_point(row,source,order,policy)
    require(row["id"]==source["id"] && row["momentum_order"]===order,"changed point identity or momentum order")
    require(Set(keys(row["pairs"]))==Set(PAIRS) && Set(keys(row["categories"]))==Set(("U","L","C")),"incomplete yields")
    for value in values(row["pairs"])
        require(Set(keys(value))==Set(("UU","Collins")) && all(finite_number,values(value)) && value["UU"]>0,"invalid absolute pair yields")
    end
    pairs=row["pairs"]
    u=Dict(key=>pairs["pm"][key]+pairs["mp"][key] for key in ("UU","Collins"))
    l=Dict(key=>pairs["pp"][key]+pairs["mm"][key] for key in ("UU","Collins"))
    categories=Dict("U"=>u,"L"=>l,"C"=>Dict(key=>u[key]+l[key] for key in ("UU","Collins")))
    require(row["categories"]==categories,"category sums changed")
    require(all(finite_number,packed_point(row)),"nonfinite packed value")
    angular=source["source_angular"]["value"]
    require(finite_number(angular) && 0<=angular<=1 && row["angular_factor"]===angular,"changed source angular scalar")
    right=categories[source["ratio"]=="UL" ? "L" : "C"]
    amplitude=angular*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"])
    require(finite_number(row["algebra_reference"]) && row["algebra_reference"]==amplitude,"changed saved ratio algebra")
    require(scaled([row["amplitude"]],[amplitude],policy["relative_floor"])<=policy["observable_algebra_max_scaled"],"invalid ratio of integrals")
    require(all(finite_number,row["qT_bounds_GeV"]) && row["qT_bounds_GeV"]==[0.,3.5] &&
        all(finite_number,row["mapped_P_bounds_GeV"]) && row["mapped_P_bounds_GeV"]==[0.,source["kinematics"]["z1"]*3.5],"changed finite window or Jacobian map")
    row
end

function collect_level(request,records,level,B,parent)
    validate_parent(parent)
    require(level isa Int && level in 1:4 && B isa Int && B in (12,24,32),"undeclared level or finite bound")
    groups,keys=geometry_groups(request);policy=parent["numerical_policy"]
    require(length(records)==36,"all36 geometry records are required")
    outputs=[Vector{Any}(undef,72) for _ in policy["momentum_orders"]]
    for (index,record) in enumerate(records)
        require(record["status"]==VALUES_STATUS && get(record,"scope",nothing)==SCOPE &&
            get(record,"theory_scope",nothing)==declaration(request) &&
            all(get(record,key,nothing)===false for key in (FLAGS...,"accepted_fitting_gate")),"wrong or self-admitting geometry values")
        require(record["geometry_index"]===index && record["native_level"]===level && record["bound_inverse_GeV"]===B &&
            record["setting"]==parent["numerical_policy"]["refinements"][level] && Tuple(record["geometry"])==keys[index],"changed geometry, native setting or bound")
        levels=record["momentum_levels"]
        require(length(levels)==length(policy["momentum_orders"]) && [row["order"] for row in levels]==policy["momentum_orders"],"retain all declared momentum orders")
        for (oi,order) in enumerate(policy["momentum_orders"])
            require(length(levels[oi]["points"])==2,"both source ratio outputs required")
            for (row,point_index) in zip(levels[oi]["points"],groups[index])
                outputs[oi][point_index]=validate_point(row,request["points"][point_index],order,policy)
            end
        end
    end
    [Dict("order"=>order,"points"=>outputs[index]) for (index,order) in enumerate(policy["momentum_orders"])]
end

function select_momentum(request,records,level,B,parent)
    policy=parent["numerical_policy"];levels=collect_level(request,records,level,B,parent)
    selected=nothing
    for (index,row) in enumerate(levels)
        error=index==1 ? nothing : scaled(packed(row["points"]),packed(levels[index-1]["points"]),policy["relative_floor"])
        met=error!==nothing && error<=policy["momentum_max_scaled"]
        accepted=row["order"]>=policy["minimum_momentum_order"] && met
        merge!(row,Dict("previous_order_max_scaled"=>error,"target_met"=>met,"accepted"=>accepted))
        selected===nothing && accepted && (selected=index)
    end
    result=merge(false_flags(),Dict{String,Any}("status"=>"MOMENTUM_SELECTION_OBSERVATION_NOT_GATE",
        "scope"=>SCOPE,"theory_scope"=>declaration(request),
        "native_level"=>level,"bound_inverse_GeV"=>B,"momentum_levels"=>levels,"selected_order"=>nothing,
        "selection_complete"=>selected!==nothing,"points"=>Any[]))
    if selected!==nothing
        result["selected_order"]=levels[selected]["order"];result["points"]=levels[selected]["points"]
    end
    result
end

function select_native(request,native_records,parent)
    require(length(native_records) in (3,4),"all native1/2/3 records required before global selection")
    policy=parent["numerical_policy"];levels=Any[];chosen=nothing
    for (level,records) in enumerate(native_records)
        result=select_momentum(request,records,level,24,parent)
        require(result["selection_complete"],"declared momentum refinement exhausted")
        error=level==1 ? nothing : scaled(packed(result["points"]),packed(last(levels)["points"]),policy["relative_floor"])
        met=error!==nothing && error<=policy["physical_value_max_scaled"]
        merge!(result,Dict("previous_native_level_max_scaled"=>error,"native_target_met"=>met))
        push!(levels,result)
        if level>=policy["minimum_refinement_level"] && met
            chosen=level
            require(level==length(native_records),"native4 may run only after the complete native3 selector requests it")
        end
    end
    merge(false_flags(),Dict{String,Any}("scope"=>SCOPE,"theory_scope"=>declaration(request),"status"=>chosen===nothing ?
        (length(levels)==3 ? "NATIVE4_REQUIRED_NOT_DIAGNOSTIC_GATE" : "NATIVE_REFINEMENT_EXHAUSTED_NOT_GATE") :
        "COMMON_NATIVE_SELECTION_OBSERVATION_NOT_GATE","selection_complete"=>chosen!==nothing,
        "selected_level"=>chosen,"levels"=>levels))
end

function select_tails(request,native_records,records12,records32,parent)
    native=select_native(request,native_records,parent)
    require(native["selection_complete"],"complete common native selection required before finite tails")
    level=native["selected_level"];b24=last(native["levels"])
    b12=select_momentum(request,records12,level,12,parent);b32=select_momentum(request,records32,level,32,parent)
    require(b12["selection_complete"] && b32["selection_complete"],"finite-bound momentum refinement exhausted")
    policy=parent["numerical_policy"]
    error=scaled(packed(b32["points"]),packed(b24["points"]),policy["relative_floor"])
    merge(false_flags(),Dict{String,Any}("status"=>"FINITE_TAIL_SELECTION_OBSERVATION_NOT_GATE",
        "scope"=>SCOPE,"theory_scope"=>declaration(request),
        "selected_level"=>level,"finite24to32_max_scaled"=>error,
        "finite_tail_target_met"=>error<=policy["physical_value_max_scaled"],
        "infinite_tail_proved"=>false,"bounds"=>[b12,b24,b32]))
end
end
