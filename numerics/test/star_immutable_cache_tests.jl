using Test, CollinsFit
if !isdefined(@__MODULE__,:JointLikelihoodReference)
    include("../bin/JointLikelihoodReference.jl")
end

function immutable_born_tree(value)
    value isa Union{Int,Float64,Symbol} && return true
    value isa Tuple && return all(immutable_born_tree,value)
    value isa Union{JetBornPlan,JetBornKinematics,CollinsFit.JetBornTerm} || return false
    !ismutabletype(typeof(value)) && all(n->immutable_born_tree(getfield(value,n)),fieldnames(typeof(value)))
end

@testset "immutable Born binding cannot hide mutable point data" begin
    p=JointLikelihoodReference.synthetic_star(200.)
    @test immutable_born_tree(p.born)
    @test_throws ErrorException setfield!(first(p.born),:fragmentation,:quark_and_gluon)
    @test_throws ErrorException setfield!(first(first(p.born).terms),:unpolarized,99.)
    @test_throws ErrorException setfield!(p,:born_signature,"invented")
    fields=(p.weights,p.couplings,p.fA,p.fB,p.operator,p.samples,p.fragmentation,p.unpolarized_coefficients,
        p.transfer_coefficients,p.jet_radius,p.axis,p.scale_prescription,p.measure_source,p.provenance,p.signature)
    external=collect(p.born)
    owned=StarPointPlan(external,fields...)
    external[1]=external[2]
    @test owned.born[1]===p.born[1] && owned.born[1]!==external[1]
    @test CollinsFit.star_point_valid(owned)===nothing
    @test_throws MethodError StarPointPlan(p.born,fields[1:end-1]...,"invented",p.signature)
    for field in (:weights,:couplings,:fA,:fB,:unpolarized_coefficients,:transfer_coefficients)
        damaged=deepcopy(p);array=getfield(damaged,field);array[1]=nextfloat(array[1])
        @test_throws ErrorException CollinsFit.star_point_valid(damaged)
    end
    for field in (:matched,:uu_weights,:collins_weights,:quark_exponent,:gluon_exponent,:shape,:unpolarized)
        damaged=deepcopy(p);array=getfield(damaged.fragmentation,field);array[1]=nextfloat(array[1])
        @test_throws ErrorException CollinsFit.star_point_valid(damaged)
    end
    for location in (:point,:fragmentation)
        damaged=deepcopy(p)
        metadata=location==:point ? damaged.provenance : damaged.fragmentation.provenance
        metadata["changed provenance"]="rejected"
        @test_throws ErrorException CollinsFit.star_point_valid(damaged)
    end
end
