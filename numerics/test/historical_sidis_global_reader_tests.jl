using Test, CollinsFit, JSON
isdefined(@__MODULE__,:HistoricalSIDISJobRecords) || include(joinpath(@__DIR__,"../bin/HistoricalSIDISJobRecords.jl"))
isdefined(@__MODULE__,:synthetic_sidis_job_reference_records) || include(joinpath(@__DIR__,"historical_sidis_job_fixtures.jl"))

@testset "Independent global reader requires all132 signed native comparisons" begin
    ids=["SYNTHETIC_READER:"*sidis_job_point_key(i) for i in 1:132]
    ctx=(contract=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT),request=(data=Dict("points"=>[Dict("id"=>id) for id in ids]),))
    function row(level)
        uu=1+level*1e-6;fut=.2+level*1e-6
        Dict("level"=>level,"bound_inverse_GeV"=>24,"values"=>Dict("UU"=>uu,"FUT"=>fut,"amplitude"=>fut/uu))
    end
    baseline=[Dict{String,Any}("phase"=>"native123","point"=>Dict("id"=>id),"levels"=>[row(l) for l in 1:3]) for id in ids]
    @test HistoricalSIDISJobRecords.native_reader(ctx,baseline).level==3
    @test_throws ErrorException HistoricalSIDISJobRecords.native_reader(ctx,baseline[1:131])
    @test_throws ErrorException HistoricalSIDISJobRecords.native_reader(ctx,reverse(baseline))
    baseline[end]["levels"][3]["values"]["FUT"]=-.25
    baseline[end]["levels"][3]["values"]["amplitude"]=-.25/baseline[end]["levels"][3]["values"]["UU"]
    @test HistoricalSIDISJobRecords.native_reader(ctx,baseline).status=="NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"
    fourth=[Dict("phase"=>"native4","point"=>Dict("id"=>id),"levels"=>[deepcopy(p["levels"][3])]) for (id,p) in zip(ids,baseline)]
    for packet in fourth;packet["levels"][1]["level"]=4;end
    fourth[end]["levels"][1]["values"]["FUT"]=-.250001
    fourth[end]["levels"][1]["values"]["amplitude"]=-.250001/fourth[end]["levels"][1]["values"]["UU"]
    @test HistoricalSIDISJobRecords.native_reader(ctx,baseline;fourth).level==4
    fourth[end]["levels"][1]["values"]["FUT"]=-.3
    fourth[end]["levels"][1]["values"]["amplitude"]=-.3/fourth[end]["levels"][1]["values"]["UU"]
    @test HistoricalSIDISJobRecords.native_reader(ctx,baseline;fourth).status=="FAIL"
end

@testset "Independent saved-value links reject inconsistent source amplitudes" begin
    ctx=(contract=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT),)
    values=Dict("UU"=>2.,"FUT"=>.2,"amplitude"=>.05)
    packet=Dict{String,Any}("values"=>deepcopy(values),"source_factor"=>.5,
        "tails"=>[Dict("values"=>deepcopy(values)) for _ in 1:3],
        "candidate_derivatives"=>Dict("value"=>.05),"reference_centre"=>Dict("value"=>"0.05"),
        "independent_scalar_reference"=>deepcopy(values))
    packet["scalar_checks"]=HistoricalSIDISJobRecords.scalar_checks(values,values,nothing,ctx)
    @test HistoricalSIDISJobRecords.audit_value_links(packet,ctx)
    bad=deepcopy(packet);bad["tails"][3]["values"]["amplitude"]*=2
    @test_throws ErrorException HistoricalSIDISJobRecords.audit_value_links(bad,ctx)
    bad=deepcopy(packet);bad["candidate_derivatives"]["value"]*=2
    @test_throws ErrorException HistoricalSIDISJobRecords.audit_value_links(bad,ctx)
    bad=deepcopy(packet);bad["reference_centre"]["value"]="0.1"
    @test_throws ErrorException HistoricalSIDISJobRecords.audit_value_links(bad,ctx)
end

@testset "Independent global reader preserves adjacent all13 refinement failures" begin
    ordinary=synthetic_sidis_job_reference_records(7)
    adjacent=synthetic_sidis_job_reference_records(3000)
    unresolved=synthetic_sidis_job_reference_records(30000)
    ids=["SYNTHETIC_READER:"*sidis_job_point_key(i) for i in 1:132]
    ctx=(contract=JSON.parsefile(CollinsFit.HISTORICAL_SIDIS_CONTRACT),request=(data=Dict("points"=>[Dict("id"=>id) for id in ids]),))
    packet(id,reference,levels)=Dict{String,Any}("point"=>Dict("id"=>id),"native_level"=>3,
        "reference_steps"=>deepcopy(reference.rows[levels]),"centre_errors_decimal"=>reference.errors)
    initial=[packet(id,ordinary,1:3) for id in ids]
    @test HistoricalSIDISJobRecords.derivative_reader(ctx,initial,3).level==3
    initial[end]=packet(ids[end],adjacent,1:3)
    @test HistoricalSIDISJobRecords.derivative_reader(ctx,initial,3).status=="DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"
    fourth=[packet(id,ordinary,4:4) for id in ids];fourth[end]=packet(ids[end],adjacent,4:4)
    @test HistoricalSIDISJobRecords.derivative_reader(ctx,initial,3;fourth).level==4
    initial[end]=packet(ids[end],unresolved,1:3);fourth[end]=packet(ids[end],unresolved,4:4)
    @test HistoricalSIDISJobRecords.derivative_reader(ctx,initial,3;fourth).status=="FAIL"
    @test_throws ErrorException HistoricalSIDISJobRecords.derivative_reader(ctx,initial[1:131],3)
    @test_throws ErrorException HistoricalSIDISJobRecords.derivative_reader(ctx,initial,3;fourth=fourth[1:131])
    initial[end]["native_level"]=4
    @test_throws ErrorException HistoricalSIDISJobRecords.derivative_reader(ctx,initial,3)
end
