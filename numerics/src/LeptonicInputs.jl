# Prepare fixed unpolarized matching matrices outside the fitted-parameter loop.
function matched_pdf_inputs(input::CT10Input,model,radial,x,flavors,order)
    gridfile = verified_parton_asset(input.input_root,"CT10_central")
    knots = parse.(Float64,split(readlines(gridfile)[4]))
    q2 = [css_mu_b2(model,b) for b in radial.b]
    rules = [matching_rule(:SIDIS,:pdf,x,alpha_s(model.coupling,q);order=order,knots=knots) for q in q2]
    coordinates = first(rules).coordinates
    count,nf = length(coordinates),length(flavors)
    output = zeros(nf,length(radial.b))
    # Bound temporary provider arrays; the prepared output is only flavors x b.
    block = max(1,32768 ÷ count)
    for low in 1:block:length(q2)
        high = min(length(q2),low+block-1)
        xs = repeat(coordinates,high-low+1)
        qs = repeat(q2[low:high];inner=count)
        values = reshape(ct10_batch(input,vcat(flavors,21),xs,qs),nf+1,count,high-low+1)
        for (j,bindex) in enumerate(low:high)
            gluon = values[end,:,j]./coordinates
            for f in 1:nf
                output[f,bindex] = match_leg(rules[bindex],values[f,:,j]./coordinates;gluon=gluon)
            end
        end
    end
    output
end

function matched_ff_inputs(input::LegacyPartonInput,process,model,radial,z,charge,flavors,order,policy)
    input.kind == :DSS14 || error("DSS14 pion fragmentation input required")
    output = zeros(length(flavors),length(radial.b))
    for (j,b) in enumerate(radial.b)
        q2 = css_mu_b2(model,b)
        rule = matching_rule(process,:ff,z,alpha_s(model.coupling,q2);order=order,knots=input.x)
        values = hcat([pion_ff(input,vcat(flavors,21),charge,y,q2;policy=policy) for y in rule.coordinates]...)
        for f in eachindex(flavors)
            output[f,j] = match_leg(rule,view(values,f,:);gluon=view(values,size(values,1),:))
        end
    end
    output
end

function prepared_input_metadata(process,model,radial,order,policy,audits)
    Dict("scope"=>"reconstructed_theory","production_ready"=>false,"process"=>string(process),
        "matching_order_per_input_cell"=>order,"radial_order_per_panel"=>radial.order,
        "radial_bound_inverse_GeV"=>last(radial.edges),"range_policy"=>string(policy),
        "input_manifest_sha256"=>bytes2hex(sha256(read(PARTON_MANIFEST))),
        "scale_range_GeV2"=>extrema([css_mu_b2(model,b) for b in radial.b]),
        "source_range_audits"=>audits,"historical_fit"=>false)
end

"Pinned inputs, fixed unpolarized SIDIS legs; no data likelihood or polarized model is inferred."
function prepare_unpolarized_sidis(pdf::CT10Input,ff::LegacyPartonInput,model,radial,Q2,P,x,z,charge;
                                   matching_order=8,range_policy=:strict)
    pdf.policy==range_policy || error("CT10 and requested range policies must agree explicitly")
    plan = leptonic_plan(:SIDIS,model,radial,Q2,P,z)
    first = matched_pdf_inputs(pdf,model,radial,css_z(x),plan.flavors,matching_order)
    second = matched_ff_inputs(ff,:SIDIS,model,radial,z,charge,plan.flavors,matching_order,range_policy)
    q2 = [css_mu_b2(model,b) for b in radial.b]
    audit = Dict("CT10_outside_strict_nodes"=>count(q -> !parton_domain(pdf,x,q)["within_grid"],q2),
        "DSS14_outside_strict_nodes"=>count(q -> !parton_domain(ff,z,q;policy=range_policy)["within_grid"],q2),
        "radial_nodes"=>length(q2))
    (plan=plan,first=first,second=second,
     provenance=prepared_input_metadata(:SIDIS,model,radial,matching_order,range_policy,audit))
end

"Pinned pion FF legs for the hadron-2-frame SIA unpolarized structure function."
function prepare_unpolarized_sia(ff::LegacyPartonInput,model,radial,Q2,P,z1,z2,charge1,charge2;
                                matching_order=8,range_policy=:strict)
    range_policy in (:strict,:reference) || error("explicit known range policy required")
    plan = leptonic_plan(:SIA,model,radial,Q2,P,z1;z2=z2)
    first = matched_ff_inputs(ff,:SIA,model,radial,z1,charge1,plan.flavors,matching_order,range_policy)
    second = matched_ff_inputs(ff,:SIA,model,radial,z2,charge2,plan.flavors,matching_order,range_policy)
    q2 = [css_mu_b2(model,b) for b in radial.b]
    audit = Dict("DSS14_outside_strict_nodes"=>count(q -> !parton_domain(ff,z1,q;policy=range_policy)["within_grid"],q2),
        "radial_nodes"=>length(q2))
    (plan=plan,first=first,second=second,
     provenance=prepared_input_metadata(:SIA,model,radial,matching_order,range_policy,audit))
end
