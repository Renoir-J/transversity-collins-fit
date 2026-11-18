using Test, JSON, CollinsFit
include("../bin/SIAGeometryGlobalAudit.jl")
const G=SIAGeometryGlobalAudit
module ExplicitSyntheticOperands
# Reuse the existing independently specified arithmetic fixtures and their
# controls. No native construction or experimental intake is performed.
include("sia_geometry_selection_tests.jl")
end
const F=ExplicitSyntheticOperands

function complete_values(data,parent,level;B=24,scale=(geometry,order)->1.)
    values=F.fixture_records(data,parent,level;B=B,scale=scale)
    for record in values
        merge!(record,G.N.flags(),Dict("initial_source_kind"=>"analytic",
            "same_native_rows_across_ratios_producer_attestation"=>true,
            "resources"=>[Dict("geometry"=>record["geometry"],"bound_inverse_GeV"=>B,
                "radial_nodes"=>1,"row_coefficient_payload_bytes"=>8,"cold_seconds_informational"=>0.,
                "context_provenance"=>Dict("scope"=>"explicitly_synthetic_arithmetic_not_physical_preparation"))]))
    end
    values
end

function native_reports(data,parent;fine=1.,mixed=false)
    levels=[complete_values(data,parent,level;scale=(geometry,order)->level==3 ? fine : 1.) for level in 1:3]
    if mixed
        scales=Dict(16=>1.01,32=>1.001,64=>1.000001,128=>1.000000001,256=>1.)
        levels[3]=complete_values(data,parent,3;scale=(geometry,order)->geometry==36 ? scales[order] : 1.)
    end
    [[Dict("native_values"=>[level[index] for level in levels]) for index in 1:36]]
end

@testset "Independent complete72 saved-operand selection, explicitly synthetic" begin
    parent=JSON.parsefile(joinpath(@__DIR__,"../../config/historical-sia-integral-contract.json"));data=F.fixture_request()
    reports=native_reports(data,parent);original=deepcopy(reports)
    result=G.produce_native(data,parent,reports)
    expected=G.audit_native(data,parent,reports,JSON.parse(JSON.json(result)))
    @test expected.selected==3 && length(last(expected.levels).points)==72
    @test all(row.order==64 for row in expected.levels)
    @test reports==original
    for mixed in (false,true)
        current=native_reports(data,parent;mixed=mixed)
        selected=G.produce_native(data,parent,current)
        audited=G.audit_native(data,parent,current,selected)
        @test last(audited.levels).order==(mixed ? 256 : 64)
        @test [p["id"] for p in last(audited.levels).points]==[p["id"] for p in data["points"]]
        for field in ("selected_level","status","selection_complete")
            bad=deepcopy(selected);bad[field]=field=="selected_level" ? 4 : field=="status" ? "PASS" : false
            @test_throws ErrorException G.audit_native(data,parent,current,bad)
        end
        for mutation in ("omit","order","residual","point","target","coarse-miss")
            bad=deepcopy(selected);row=last(bad["levels"])
            if mutation=="omit";pop!(bad["levels"])
            elseif mutation=="order";row["selected_order"]=128
            elseif mutation=="residual";row["momentum_levels"][3]["previous_order_max_scaled"]+=.01
            elseif mutation=="point";row["points"][72]["amplitude"]+=.01
            elseif mutation=="target";row["native_target_met"]=false
            elseif mutation=="coarse-miss";row["momentum_levels"][1]["accepted"]=true
            end
            @test_throws ErrorException G.audit_native(data,parent,current,bad)
        end
    end
    missed=native_reports(data,parent;fine=1.01)
    request4=G.produce_native(data,parent,missed)
    @test G.audit_native(data,parent,missed,request4).selected===nothing
    fine4=complete_values(data,parent,4;scale=(geometry,order)->1.01)
    push!(missed,[Dict("native_values"=>[value]) for value in fine4])
    @test G.audit_native(data,parent,missed,G.produce_native(data,parent,missed)).selected==4
    premature=deepcopy(reports);push!(premature,[Dict("native_values"=>[value]) for value in fine4])
    @test_throws ErrorException G.native(data,G.native_records(premature),parent)
    for B in (12,32),level in (3,4)
        values=complete_values(data,parent,level;B=B)
        selected=G.momentum(data,values,level,B,parent)
        @test selected.order==64
        for index in (1,15,36)
            @test G.operands(selected.points[index])==F.HSP.packed_point(selected.points[index])
        end
    end
    tails=[Dict("tail_values"=>[complete_values(data,parent,3;B=B)[i] for B in (12,32)]) for i in 1:36]
    result=G.produce_tails(data,parent,reports,tails)
    selected=G.audit_tails(data,parent,reports,tails,result)
    @test selected.error==0. && selected.hi.order==64
    for mutation in ("bound","residual","decision","infinite","B24-order","omit")
        bad=deepcopy(result)
        if mutation=="bound";bad["selected_level"]=4
        elseif mutation=="residual";bad["finite24to32_max_scaled"]=.01
        elseif mutation=="decision";bad["finite_tail_target_met"]=false
        elseif mutation=="infinite";bad["infinite_tail_proved"]=true
        elseif mutation=="B24-order";bad["bounds"][2]["selected_order"]=128
        elseif mutation=="omit";pop!(bad["bounds"])
        end
        @test_throws ErrorException G.audit_tails(data,parent,reports,tails,bad)
    end
    for (i,report) in enumerate(tails)
        report["independent_references"]=[Dict("reference"=>deepcopy(selected.hi.points[j]),
            "order_comparisons"=>[Dict("order"=>order,"max_scaled"=>0.) for order in (16,32,64,128,256)]) for j in (i,i+36)]
        for row in report["independent_references"]
            delete!(row["reference"],"momentum_order")
        end
    end
    @test G.reference_errors(data,parent,tails,selected.hi)==zeros(72)
    tails[36]["independent_references"][2]["order_comparisons"][3]["max_scaled"]=1.
    @test_throws ErrorException G.reference_errors(data,parent,tails,selected.hi)
end
