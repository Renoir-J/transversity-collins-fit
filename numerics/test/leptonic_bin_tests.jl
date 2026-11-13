# All nodes, weights and acceptances in this file are explicitly synthetic.
function test_bin(nodes,weights;source="synthetic bin fixture",measure="two fixed synthetic nodes",
                  weight_definition="synthetic fixed quadrature times prefactor and acceptance",acceptance="explicit synthetic correlated acceptance")
    leptonic_bin(nodes,weights;source=source,measure=measure,weight_definition=weight_definition,acceptance=acceptance)
end

"Reference through separate charge-resolved point structure functions."
function manual_point_bin(nodes,weights,theta)
    if first(nodes).context.leptonic.process==:SIDIS
        rows=[polarized_structure(node.context,theta) for node in nodes]
        sum(weights[i]*nodes[i].angular*rows[i].Collins for i in eachindex(nodes))/
            sum(weights[i]*rows[i].UU for i in eachindex(nodes))
    else
        records=[[polarized_structure(with_pion_charges(node.context,c...),theta)
                  for c in ((1,1),(1,-1),(-1,1),(-1,-1))] for node in nodes]
        U=sum(weights[i,1]*nodes[i].angular*(records[i][2].Collins+records[i][3].Collins) for i in eachindex(nodes))
        L=sum(weights[i,2]*nodes[i].angular*(records[i][1].Collins+records[i][4].Collins) for i in eachindex(nodes))
        DU=sum(weights[i,1]*(records[i][2].UU+records[i][3].UU) for i in eachindex(nodes))
        DL=sum(weights[i,2]*(records[i][1].UU+records[i][4].UU) for i in eachindex(nodes))
        U/DU-(first(nodes).kind==:UL ? L/DL : (U+L)/(DU+DL))
    end
end

@testset "bin composition integrates yields before normalization" begin
    theta=collect(kpsy_published_parameters().values)
    sidis=observable_fixture(:SIDIS);sia=observable_fixture(:SIA)
    for definition in (:structure_function_ratio,:sin_phi_sum_amplitude)
        nodes=[sidis_observable(with_momentum(sidis,P);definition=definition,y=y,source="synthetic bin node")
               for (P,y) in ((.1,.2),(.35,.8))]
        weights=[2.,7.];bin=test_bin(nodes,weights)
        f(t)=predict_bin(bin,t);reference(t)=manual_point_bin(nodes,weights,t)
        @test f(theta) ≈ reference(theta) rtol=1e-13 atol=1e-15
        @test ForwardDiff.gradient(f,theta) ≈ ForwardDiff.gradient(reference,theta) rtol=1e-12 atol=1e-14
        @test ForwardDiff.hessian(f,theta) ≈ ForwardDiff.hessian(reference,theta) rtol=1e-11 atol=1e-13
        @test abs(f(theta)-sum(weights.*predict_observable.(nodes,Ref(theta)))/sum(weights))>1e-8
        @test predict_bin(test_bin(nodes,37 .* weights),theta) ≈ f(theta) rtol=1e-13
        @test predict_bin(test_bin(reverse(nodes),reverse(weights)),theta) ≈ f(theta) rtol=1e-13
        @test predict_bin(test_bin([nodes[1]],[1.]),theta)==predict_observable(nodes[1],theta)
        @test bin.provenance["production_ready"]===false && bin.provenance["G2PERF_pass"]===false
        original=predict_bin(bin,theta)
        weights[1]=99.;nodes[1]=nodes[2]
        @test predict_bin(bin,theta)==original # The caller's arrays are not the owned arrays.
    end
    for ratio in (:UL,:UC)
        nodes=[sia_observable(with_momentum(sia,P);ratio=ratio,
            angular=SIAAngularFactor(a;source="synthetic conditional angle"),source="synthetic bin node")
            for (P,a) in ((.1,.15),(.35,.9))]
        weights=[1. .3;7. 5.];bin=test_bin(nodes,weights)
        f(t)=predict_bin(bin,t);reference(t)=manual_point_bin(nodes,weights,t)
        @test f(theta) ≈ reference(theta) rtol=1e-13 atol=1e-15
        @test ForwardDiff.gradient(f,theta) ≈ ForwardDiff.gradient(reference,theta) rtol=1e-12 atol=1e-14
        @test ForwardDiff.hessian(f,theta) ≈ ForwardDiff.hessian(reference,theta) rtol=1e-11 atol=1e-13
        @test abs(f(theta)-sum(weights[:,1].*predict_observable.(nodes,Ref(theta)))/sum(weights[:,1]))>1e-8
        @test predict_bin(test_bin(nodes,37 .* weights),theta) ≈ f(theta) rtol=1e-13
        @test predict_bin(test_bin(reverse(nodes),reverse(weights;dims=1)),theta) ≈ f(theta) rtol=1e-13
        @test predict_bin(test_bin([nodes[1]],ones(1,2)),theta) ≈ predict_observable(nodes[1],theta) rtol=1e-13
        @test predict_bin(test_bin([nodes[1]],reshape([2.,5.],1,2)),theta) ≈
            manual_point_bin([nodes[1]],reshape([2.,5.],1,2),theta) rtol=1e-13
    end
end

@testset "bin conventions, signed numerators and finite momentum integrals" begin
    theta=collect(kpsy_published_parameters().values);context=observable_fixture(:SIDIS)
    nodes(phase)=[sidis_observable(with_momentum(context,P),SIDISAsymmetryConvention(
        phase=phase,normalization=:fixed_y_DNN_included,y=y,source="synthetic bin phase"))
        for (P,y) in ((.1,.2),(.35,.8))]
    positive=nodes(:sin_phi_sum);negative=nodes(:sin_phi_sum_minus_pi)
    @test predict_bin(test_bin(positive,[2.,7.]),theta)==-predict_bin(test_bin(negative,[2.,7.]),theta)
    @test_throws ErrorException test_bin([positive[1],negative[2]],[1.,1.])
    for process in (:SIDIS,:SIA)
        context=observable_fixture(process)
        momentum=MomentumRule(.1,.3;order=16)
        node=process==:SIDIS ? sidis_observable(context;definition=:structure_function_ratio,momentum=momentum,source="synthetic integral") :
            sia_observable(context;ratio=:UL,angular=SIAAngularFactor(.6;source="synthetic"),momentum=momentum,source="synthetic integral")
        weights=process==:SIDIS ? [1.] : ones(1,2)
        @test predict_bin(test_bin([node],weights),theta) ≈ predict_observable(node,theta) rtol=1e-13
    end
end

@testset "bin preparation and complete mutation guards reject invalid evidence" begin
    theta=collect(kpsy_published_parameters().values);context=observable_fixture(:SIDIS)
    node=sidis_observable(context;definition=:structure_function_ratio,source="synthetic guard")
    sia=sia_observable(observable_fixture(:SIA);ratio=:UL,angular=SIAAngularFactor(.6;source="synthetic"),source="synthetic")
    for weights in ([1.,-1.],[0.,0.],[1.,NaN],[1.,Inf],[true,false],[1.],ones(2,1))
        @test_throws ErrorException test_bin([node,node],weights)
    end
    for weights in ([1.,1.],ones(1,2),[1. 0.;2. 0.],[1. -1.;2. 3.],[1. NaN;2. 3.])
        @test_throws ErrorException test_bin([sia,sia],weights)
    end
    @test_throws ErrorException test_bin(LeptonicObservablePlan[],Float64[])
    @test_throws ErrorException test_bin([node,sia],[1.,1.])
    for keyword in (:source,:measure,:weight_definition,:acceptance)
        @test_throws ErrorException test_bin([node],[1.];(;keyword=>" ")...)
    end
    @test_throws ErrorException ForwardDiff.derivative(t->predict_bin(test_bin([node],[1+t]),theta),0.)
    different=sidis_observable(with_pion_charges(context,-1);definition=:structure_function_ratio,source="different charge")
    @test_throws ErrorException test_bin([node,different],[1.,1.])
    integrated=sidis_observable(context;definition=:structure_function_ratio,momentum=MomentumRule(.1,.3),source="integrated")
    @test_throws ErrorException test_bin([node,integrated],[1.,1.])
    mean_node=sia_observable(sia.context;ratio=:UL,angular=SIAMeanOfRatio(.6;source="synthetic",averaging_measure="synthetic event weights"),source="synthetic")
    @test_throws ErrorException test_bin([sia,mean_node],ones(2,2))
    for position in (1,17,65,129)
        bin=test_bin(fill(node,129),ones(129));bin.weights[position,1]=nextfloat(1.)
        @test_throws ErrorException predict_bin(bin,theta)
    end
    bin=test_bin([node],[1.]);bin.denominators[1]+=1
    @test_throws ErrorException predict_bin(bin,theta)
    bin=test_bin([node],[1.]);bin.provenance["production_ready"]=true
    @test_throws ErrorException predict_bin(bin,theta)
    bin=test_bin([node],[1.]);bin.nodes[1]=different
    @test_throws ErrorException predict_bin(bin,theta)
    bin=test_bin([node],[1.]);node.collins_weights[1]+=1
    @test_throws ErrorException predict_bin(bin,theta)
end
