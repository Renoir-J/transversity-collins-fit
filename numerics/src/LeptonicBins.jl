"Fixed prepared quadrature nodes; numerator and denominator integrals stay separate."
struct LeptonicBinPlan
    nodes::Vector{LeptonicObservablePlan}
    weights::Matrix{Float64}
    denominators::Vector{Float64}
    provenance::Dict{String,Any}
    signature::String
end

function bin_model_settings(model)
    (model.bmax,model.Q02,model.g2,model.gq,model.gh,
     model.coupling.lambda_GeV,model.coupling.nf,model.nodes,model.weights)
end

function bin_node_identity(node)
    context=node.context;target=context.target
    target_data=target===nothing ? nothing : (target.name,target.unpolarized,target.polarized,target.source)
    (context.leptonic.process,node.kind,bin_model_settings(context.leptonic.model),
     context.leptonic.flavors,context.leptonic.antiparticle,context.leptonic.charge2,
     target_data,context.leptonic.process==:SIDIS ? context.charges : nothing,
     get(node.provenance,"angle_convention",nothing),get(node.provenance,"normalization",nothing),
     get(node.provenance,"depolarization_removed",nothing),get(node.provenance,"angular_definition",nothing),
     node.momentum===nothing)
end

bin_signature(nodes,weights,denominators,metadata)=numerical_snapshot(
    (Tuple(node.signature for node in nodes),weights,denominators,metadata))

"""
Compose already prepared nodes with explicitly supplied fixed quadrature,
acceptance and cross-section weights. SIDIS needs one weight per node. SIA
needs an N-by-2 matrix of unlike/like weights; C is their summed yield.
No acceptance, phase-space Jacobian, sigma0 or experimental admission is inferred.
"""
function leptonic_bin(nodes,weights;source,measure,weight_definition,acceptance)
    nodes isa AbstractVector && !isempty(nodes) && all(node->node isa LeptonicObservablePlan,nodes) ||
        error("nonempty vector of prepared leptonic observable nodes required")
    descriptions=(source,measure,weight_definition,acceptance)
    all(text->text isa AbstractString && !isempty(strip(text)),descriptions) ||
        error("nonblank source, measure, weight definition and acceptance description required")
    foreach(observable_valid,nodes)
    identity=bin_node_identity(first(nodes))
    all(node->bin_node_identity(node)==identity,nodes) ||
        error("bin nodes must share process, model, target/species, observable convention and point/P-integral mode")
    sidis=first(nodes).context.leptonic.process==:SIDIS
    count=sidis ? 1 : 2
    raw=if sidis
        weights isa AbstractVector && length(weights)==length(nodes) || error("one fixed SIDIS weight per node required")
        reshape([fixed_scalar(w,"SIDIS quadrature weight") for w in weights],:,1)
    else
        weights isa AbstractMatrix && size(weights)==(length(nodes),2) || error("N-by-2 U/L weight matrix required")
        [fixed_scalar(weights[i,j],"SIA category quadrature weight") for i in eachindex(nodes),j in 1:2]
    end
    all(>=(0),raw) && all(j->sum(raw[:,j])>0,1:count) ||
        error("nonnegative finite fixed weights with nonzero support in every category required")
    owned=LeptonicObservablePlan[node for node in nodes]
    denominators=[sum(raw[i,j]*owned[i].denominators[j] for i in eachindex(owned)) for j in 1:count]
    foreach(physical_denominator,denominators)
    metadata=Dict{String,Any}("contract_id"=>"KPSY-LEPTONIC-BIN-001","scope"=>"reconstructed_theory",
        "production_ready"=>false,"physical_extraction"=>false,"G2PERF_pass"=>false,
        "source"=>String(source),"measure"=>String(measure),"weight_definition"=>String(weight_definition),
        "acceptance"=>String(acceptance),"process"=>string(first(owned).context.leptonic.process),
        "definition"=>string(first(owned).kind),"node_count"=>length(owned),
        "weight_columns"=>(sidis ? ["SIDIS"] : ["U","L"]),
        "prescription"=>"integrate signed angular-weighted numerators and unpolarized yields before ratios",
        "node_provenance"=>[deepcopy(node.provenance) for node in owned])
    LeptonicBinPlan(owned,copy(raw),denominators,metadata,bin_signature(owned,raw,denominators,metadata))
end

function bin_valid(bin::LeptonicBinPlan)
    foreach(observable_valid,bin.nodes)
    bin_signature(bin.nodes,bin.weights,bin.denominators,bin.provenance)==bin.signature ||
        error("mutated bin nodes, weights, denominators or metadata")
end

function predict_bin(bin::LeptonicBinPlan,p::KPSYParameters)
    bin_valid(bin)
    totals=fill(zero(p.values[1]),length(bin.denominators))
    for (i,node) in enumerate(bin.nodes)
        integrals=observable_polarized_integrals(node,p)
        for j in eachindex(totals)
            totals[j]+=bin.weights[i,j]*node.angular*integrals[j]
        end
    end
    if first(bin.nodes).context.leptonic.process==:SIDIS
        totals[1]/bin.denominators[1]
    else
        U,L=totals;DU,DL=bin.denominators
        right=first(bin.nodes).kind==:UL ? L/DL : (U+L)/(DU+DL)
        U/DU-right
    end
end
predict_bin(bin::LeptonicBinPlan,theta::AbstractVector)=predict_bin(bin,
    KPSYParameters(theta;gh=first(bin.nodes).context.leptonic.model.gh))
