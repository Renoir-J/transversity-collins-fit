using SpecialFunctions: besselj

"Explicit finite radial interval [0,B]; all evaluation nodes are strictly positive."
struct RadialRule
    b::Vector{Float64}
    db::Vector{Float64}
    edges::Vector{Float64}
    order::Int
end

function RadialRule(edges;order=16)
    grid = [fixed_scalar(v,"radial edge") for v in edges]
    length(grid)>=2 && first(grid)==0 && all(diff(grid).>0) ||
        error("radial edges must increase from zero to an explicit finite bound")
    nodes,weights = gauss_rule(order)
    b,db = Float64[],Float64[]
    for i in 1:length(grid)-1
        low,high = grid[i],grid[i+1]
        half = (high-low)/2
        append!(b,(low+high)/2 .+ half.*nodes)
        append!(db,half.*weights)
    end
    all(>(0),b) || error("radial nodes must be positive")
    RadialRule(b,db,grid,order)
end

function hankel_weights(rule::RadialRule,n,k)
    n isa Integer && !(n isa Bool) && n in 0:2 || error("only harmonics 0,1,2 supported")
    k = fixed_scalar(k,"P/z")
    k >= 0 || error("nonnegative transverse momentum required")
    rule.db .* rule.b.^(n+1) .* besselj.(n,k.*rule.b) ./ (2pi)
end

function finite_transform(weights,amplitude)
    length(weights)==length(amplitude) && all(isfinite,amplitude) || error("invalid radial amplitudes")
    value = dot(weights,amplitude)
    isfinite(value) || error("nonfinite radial transform")
    value
end

"KPSY Eq73/74, for amplitudes that already include their complete CSS factors."
function sidis_fourier(rule::RadialRule,P,z,uu,collins)
    z = css_z(z)
    k = fixed_scalar(P,"P")/z
    (UU=finite_transform(hankel_weights(rule,0,k),uu)/z^2,
     Collins=-finite_transform(hankel_weights(rule,1,k),collins)/(2z^3))
end

"KPSY Eq108/109 in the hadron-2 frame; P belongs to hadron 1."
function sia_fourier(rule::RadialRule,P,z1,z2,uu,collins)
    z1,z2 = css_z(z1),css_z(z2)
    k = fixed_scalar(P,"P")/z1
    (UU=finite_transform(hankel_weights(rule,0,k),uu)/z1^2,
     Collins=finite_transform(hankel_weights(rule,2,k),collins)/(4z1^3*z2))
end

struct LeptonicPlan
    process::Symbol
    model::KPSYCSS
    radial::RadialRule
    Q2::Float64
    P::Float64
    z1::Float64
    z2::Union{Nothing,Float64}
    flavors::Vector{Int}
    antiparticle::Vector{Int}
    charge2::Vector{Float64}
    uu_weights::Vector{Float64}
    collins_weights::Vector{Float64}
    uu_damping::Vector{Float64}
    fixed_exponent::Vector{Float64}
    fragmentation_shape::Vector{Float64}
end

"Cache fixed kinematics/CSS; supplied matching matrices are flavors by b nodes."
function leptonic_plan(process,model::KPSYCSS,radial::RadialRule,Q2,P,z1;
                       z2=nothing,flavors=[-5,-4,-3,-2,-1,1,2,3,4,5])
    process in (:SIDIS,:SIA) || error("unsupported leptonic process")
    q2,p,z = fixed_scalar(Q2,"Q2"),fixed_scalar(P,"P"),css_z(z1)
    p >= 0 || error("nonnegative P required")
    coupling_log(model.coupling,q2)
    second_z = if process == :SIA
        isnothing(z2) && error("SIA requires z2")
        css_z(z2)
    else
        isnothing(z2) || error("SIDIS has a single fragmentation fraction")
        nothing
    end
    !isempty(flavors) && all(f -> f isa Integer && !(f isa Bool) && 1<=abs(f)<=5,flavors) &&
        length(unique(flavors))==length(flavors) && all(f -> -f in flavors,flavors) ||
        error("unique charge-conjugate quark/antiquark pairs required")
    ids = Int.(flavors)
    anti = [findfirst(==(-f),ids) for f in ids]
    charges = [abs(f) in (2,4) ? 4/9 : 1/9 for f in ids]
    k = p/z
    uu = hankel_weights(radial,0,k)./z^2
    pol = process == :SIDIS ? -hankel_weights(radial,1,k)./(2z^3) :
                             hankel_weights(radial,2,k)./(4z^3*second_z)
    fixed = [sudakov_perturbative(model,b,q2)+css_np_common(model,b,q2)+
        (process == :SIDIS ? model.gq*b^2 : 0.) for b in radial.b]
    shape = radial.b.^2 .* (process == :SIDIS ? inv(z^2) : inv(z^2)+inv(second_z^2))
    damping = exp.(-fixed.-model.gh.*shape)
    all(isfinite,damping) || error("nonfinite CSS damping")
    LeptonicPlan(process,model,radial,q2,p,z,second_z,ids,anti,charges,uu,pol,damping,fixed,shape)
end

function matched_matrix(plan,values)
    size(values)==(length(plan.flavors),length(plan.radial.b)) && all(isfinite,values) ||
        error("matching matrix must be flavors by radial nodes with finite values")
    values
end

function flavor_product(plan,first,second)
    matched_matrix(plan,first); matched_matrix(plan,second)
    [sum(plan.charge2[f]*first[f,j]*second[plan.process==:SIA ? plan.antiparticle[f] : f,j]
        for f in eachindex(plan.flavors)) for j in eachindex(plan.radial.b)]
end

function leptonic_unpolarized(plan::LeptonicPlan,first,second)
    finite_transform(plan.uu_weights,flavor_product(plan,first,second).*plan.uu_damping)
end

function leptonic_collins(plan::LeptonicPlan,first,second,gc)
    width = collins_width(plan.model,gc)
    values = flavor_product(plan,first,second)
    damped = values .* exp.(-plan.fixed_exponent.-width.*plan.fragmentation_shape)
    finite_transform(plan.collins_weights,damped)
end

"Structure functions only: angular/depolarization factors, bins and double ratios are separate."
function leptonic_structure(plan::LeptonicPlan,first,second,polarized_first,polarized_second;gc)
    (UU=leptonic_unpolarized(plan,first,second),
     Collins=leptonic_collins(plan,polarized_first,polarized_second,gc))
end
