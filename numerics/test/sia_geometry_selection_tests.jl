using Test, JSON
include("../bin/SIAGeometrySelection.jl")
const S=SIAGeometrySelection
module FrozenPointPackingComparator
include("../bin/predict_historical_sia_points.jl")
end
const HSP=FrozenPointPackingComparator

function fixture_request()
    points=Any[]
    for view in (1002,1003),i in 1:36
        push!(points,Dict("id"=>"KPSY15_SIA_$view:r$(i+1)","ratio"=>view==1002 ? "UC" : "UL",
            "kinematics"=>Dict("z1"=>.2+.1*((i-1)÷6),"z2"=>.2+.1*((i-1)%6),"Q2_GeV2"=>112.36),
            "source_angular"=>Dict("value"=>view==1002 ? .5 : .51)))
    end
    merge(S.false_flags(),Dict{String,Any}("scope"=>"historical_reconstructed_sia_integral_prediction_diagnostic",
        "theory_scope"=>"synthetic_closure",
        "source_fingerprint"=>"01be55cc1449679f6dc26bd43ea9788604938c1fecbb25d276acaa35bc4bae47",
        "fixture_scope"=>"explicitly_synthetic_arithmetic_fixture_not_primary_data","points"=>points))
end

function fixture_point(point,order,scale)
    pairs=Dict(name=>Dict("UU"=>scale*uu,"Collins"=>scale*collins) for (name,uu,collins) in
        (("pp",2.,.01),("pm",3.,.03),("mp",1.,-.02),("mm",1.,-.01)))
    categories=HSP.category_values(pairs);u=categories["U"];right=categories[point["ratio"]=="UL" ? "L" : "C"]
    angular=point["source_angular"]["value"]
    amplitude=angular*(u["Collins"]/u["UU"]-right["Collins"]/right["UU"])
    Dict("id"=>point["id"],"pairs"=>pairs,"categories"=>categories,"amplitude"=>amplitude,
        "algebra_reference"=>amplitude,"angular_factor"=>angular,"momentum_order"=>order,
        "qT_bounds_GeV"=>[0.,3.5],"mapped_P_bounds_GeV"=>[0.,point["kinematics"]["z1"]*3.5])
end

function fixture_records(request,parent,level;B=24,scale=(geometry,order)->1.)
    groups,keys=S.geometry_groups(request)
    [merge(S.false_flags(),Dict{String,Any}("status"=>S.VALUES_STATUS,"geometry_index"=>index,
        "scope"=>S.SCOPE,"theory_scope"=>"synthetic_closure",
        "native_level"=>level,"bound_inverse_GeV"=>B,"setting"=>deepcopy(parent["numerical_policy"]["refinements"][level]),
        "geometry"=>collect(keys[index]),"momentum_levels"=>[Dict("order"=>order,
            "points"=>[fixture_point(request["points"][i],order,scale(index,order)) for i in group])
            for order in parent["numerical_policy"]["momentum_orders"]])) for (index,group) in enumerate(groups)]
end

@testset "Synthetic complete72 common-order arithmetic versus inherited packing" begin
    parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"));request=fixture_request()
    records=fixture_records(request,parent,3);before=deepcopy(records)
    result=S.select_momentum(request,records,3,24,parent)
    @test result["selected_order"]==64
    @test result["theory_scope"]=="synthetic_closure"
    @test length(result["points"])==72
    @test [p["id"] for p in result["points"]]==[p["id"] for p in request["points"]]
    @test all(S.packed_point(p)==HSP.packed_point(p) for p in result["points"])
    @test records==before
    for (field,value) in (("physical_value_max_scaled",1.),("momentum_max_scaled",1.),("relative_floor",1.),("minimum_refinement_level",true))
        badpolicy=deepcopy(parent);badpolicy["numerical_policy"][field]=value
        @test_throws ErrorException S.select_momentum(request,records,3,24,badpolicy)
    end
    scales=Dict(16=>1.01,32=>1.001,64=>1.000001,128=>1.000000001,256=>1.)
    mixed=fixture_records(request,parent,3;scale=(geometry,order)->geometry==36 ? scales[order] : 1.)
    global_result=S.select_momentum(request,mixed,3,24,parent)
    @test global_result["selected_order"]==256
    levels=global_result["momentum_levels"]
    errors=[nothing;[HSP.scaled(HSP.packed(levels[i]["points"]),HSP.packed(levels[i-1]["points"]),1e-8) for i in 2:5]]
    @test [row["previous_order_max_scaled"] for row in levels]==errors
    @test all(row["accepted_fitting_gate"]===false for row in (result,global_result))
    exhausted=fixture_records(request,parent,3;scale=(geometry,order)->1+order/1000)
    @test S.select_momentum(request,exhausted,3,24,parent)["selection_complete"]===false
    for flag in (S.FLAGS...,"accepted_fitting_gate")
        bad=deepcopy(records);bad[36][flag]=true
        @test_throws ErrorException S.select_momentum(request,bad,3,24,parent)
    end
    badscope=deepcopy(records);badscope[36]["theory_scope"]="reconstructed_theory"
    @test_throws ErrorException S.select_momentum(request,badscope,3,24,parent)
    badfixture=deepcopy(request);badfixture["theory_scope"]="reconstructed_theory"
    @test_throws ErrorException S.select_momentum(badfixture,records,3,24,parent)
    for mutation in ("omit","duplicate","geometry","ratio-point","angle","category","window","nan","order","negative")
        bad=deepcopy(records);row=bad[36]["momentum_levels"][3]["points"][2]
        if mutation=="omit";pop!(bad)
        elseif mutation=="duplicate";bad[36]=deepcopy(bad[1])
        elseif mutation=="geometry";bad[36]["geometry"][1]+=1e-14
        elseif mutation=="ratio-point";row["id"]=request["points"][36]["id"]
        elseif mutation=="angle";row["angular_factor"]=.5
        elseif mutation=="category";row["categories"]["U"]["UU"]+=.01
        elseif mutation=="window";row["qT_bounds_GeV"]=Any[false,3.5]
        elseif mutation=="nan";row["pairs"]["pp"]["Collins"]=NaN
        elseif mutation=="order";row["momentum_order"]=128
        elseif mutation=="negative";row["pairs"]["pp"]["UU"]=-1.
        end
        @test_throws ErrorException S.select_momentum(request,bad,3,24,parent)
    end
    native=[fixture_records(request,parent,level;scale=(g,o)->level==1 ? .9 : level==2 ? .99995 : 1.) for level in 1:3]
    selection=S.select_native(request,native,parent)
    @test selection["selected_level"]==3
    @test selection["selection_complete"]===true
    @test selection["levels"][2]["native_target_met"]===false
    @test_throws ErrorException S.select_native(request,native[1:2],parent)
    @test_throws ErrorException S.select_native(request,[native;[fixture_records(request,parent,4)]],parent)
    need4=[fixture_records(request,parent,level;scale=(g,o)->level==1 ? .9 : level==2 ? 1.001 : 1.) for level in 1:3]
    @test S.select_native(request,need4,parent)["status"]=="NATIVE4_REQUIRED_NOT_DIAGNOSTIC_GATE"
    @test S.select_native(request,[need4;[fixture_records(request,parent,4)]],parent)["selected_level"]==4
    failed4=fixture_records(request,parent,4;scale=(g,o)->1.01)
    @test S.select_native(request,[need4;[failed4]],parent)["status"]=="NATIVE_REFINEMENT_EXHAUSTED_NOT_GATE"
    tails=S.select_tails(request,native,fixture_records(request,parent,3;B=12),fixture_records(request,parent,3;B=32),parent)
    @test tails["finite_tail_target_met"]===true
    @test tails["infinite_tail_proved"]===false
    @test tails["bounds"][2]==last(selection["levels"])
    badtail=fixture_records(request,parent,3;B=32,scale=(g,o)->g==36 ? 1.01 : 1.)
    @test S.select_tails(request,native,fixture_records(request,parent,3;B=12),badtail,parent)["finite_tail_target_met"]===false
    for field in ("scope","source_fingerprint")
        badrequest=deepcopy(request);badrequest[field]="unreviewed"
        @test_throws ErrorException S.select_momentum(badrequest,records,3,24,parent)
    end
end
