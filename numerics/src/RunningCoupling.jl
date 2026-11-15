"A fixed-nf asymptotic two-loop coupling, distinct from an exact two-loop RG solution."
struct AsymptoticTwoLoopCoupling
    lambda_GeV::Float64
    nf::Int
    function AsymptoticTwoLoopCoupling(lambda_GeV, nf)
        lambda = fixed_scalar(lambda_GeV, "Lambda")
        lambda > 0 || error("positive Lambda required")
        nf isa Integer && 3 <= nf <= 6 || error("fixed integer nf from three to six required")
        new(lambda, nf)
    end
end

AsymptoticTwoLoopCoupling(; lambda_GeV=0.225, nf=5) = AsymptoticTwoLoopCoupling(lambda_GeV, nf)

function fixed_scalar(value, name)
    value isa Union{AbstractFloat,Integer} && !(value isa Bool) ||
        error("fixed numeric $name required; setting derivatives are not implemented")
    converted = Float64(value)
    isfinite(converted) || error("finite $name required")
    converted
end

"""Owned fixed evolution times with exact physical-scale labels and provenance.

This is numerical input, not an alpha_s function or an admission certificate.
It carries no derivatives of the grid, scales, coupling or supplied times.
"""
struct PreparedEvolutionTimes
    Q2::Tuple{Vararg{Float64}}
    tau::Tuple{Vararg{Float64}}
    Q02::Float64
    nf::Int
    alpha0::Float64
    model_id::String
    source_id::String
    source_sha256::String
    function PreparedEvolutionTimes(Q2,tau,Q02,nf,alpha0,model_id,source_id,source_sha256)
        Q2 isa Union{AbstractVector,Tuple} && tau isa Union{AbstractVector,Tuple} ||
            error("one-dimensional fixed scales and times required")
        1<=length(Q2)==length(tau)<=4096 || error("one to 4096 matching scales/times required")
        qs=Tuple(fixed_scalar(q,"prepared Q2") for q in Q2)
        ts=Tuple(fixed_scalar(t,"prepared tau") for t in tau)
        q0=fixed_scalar(Q02,"prepared Q02");a0=fixed_scalar(alpha0,"prepared alpha0 metadata")
        nf isa Integer && !(nf isa Bool) && 3<=nf<=6 || error("fixed integer nf from three to six required")
        q0>0 && 0<a0<=1 && all(>(0),qs) && all(>(0),diff(collect(qs))) || error("ordered positive prepared scales required")
        all(t->abs(t)<=1,ts) && all(>(0),diff(collect(ts))) || error("ordered prepared times in [-1,1] required")
        anchor=findfirst(==(q0),qs)
        anchor!==nothing && iszero(ts[anchor]) || error("exact Q02 label and zero anchor required")
        all(sign(q-q0)==sign(t) for (q,t) in zip(qs,ts)) || error("prepared time/scale signs disagree")
        for id in (model_id,source_id)
            id isa AbstractString && occursin(r"^[A-Za-z0-9][A-Za-z0-9_.:/@+\-]{0,255}$",id) ||
                error("bounded nonempty model and source identifiers required")
        end
        source_sha256 isa AbstractString && occursin(r"^[0-9a-f]{64}$",source_sha256) ||
            error("lowercase SHA256 source identity required")
        new(qs,ts,q0,Int(nf),a0,String(model_id),String(source_id),String(source_sha256))
    end
end

function prepared_time_metadata(profile::PreparedEvolutionTimes)
    (Q2=profile.Q2,tau=profile.tau,Q02=profile.Q02,nf=profile.nf,alpha0=profile.alpha0,
     model_id=profile.model_id,source_id=profile.source_id,source_sha256=profile.source_sha256)
end

"Exact lookup only, preserving query order and repeats; no extrapolation or time AD."
function prepared_times(profile::PreparedEvolutionTimes,Q2,Q02;nf=profile.nf)
    q0=fixed_scalar(Q02,"prepared query Q02")
    nf isa Integer && !(nf isa Bool) && nf==profile.nf && q0==profile.Q02 ||
        error("prepared nf/Q02 mismatch")
    Q2 isa Union{AbstractVector,Tuple} && 1<=length(Q2)<=4096 || error("nonempty fixed query scales required")
    [begin
        value=fixed_scalar(q,"prepared query Q2")
        index=findfirst(==(value),profile.Q2)
        index===nothing && error("unavailable prepared physical scale; interpolation is not implemented")
        profile.tau[index]
    end for q in Q2]
end

evolution_time(profile::PreparedEvolutionTimes,Q2,Q02)=only(prepared_times(profile,[Q2],Q02))

function coupling_log(coupling::AsymptoticTwoLoopCoupling, Q2)
    q2 = fixed_scalar(Q2, "Q2")
    q2 > 0 || error("positive Q2 required")
    L = log(q2) - 2log(coupling.lambda_GeV)
    L > 0 || error("Q2 is outside the coupling Landau-pole domain")
    L
end

"The specified 1/L - beta1/beta0^2*log(L)/L^2 expansion; no flavor matching."
function alpha_s(coupling::AsymptoticTwoLoopCoupling, Q2)
    L = coupling_log(coupling, Q2)
    beta0, beta1 = 11 - 2coupling.nf/3, 102 - 38coupling.nf/3
    alpha = 4pi/(beta0*L) * (1 - beta1/beta0^2 * log(L)/L)
    isfinite(alpha) && alpha > 0 || error("invalid asymptotic two-loop coupling")
    alpha
end

"Exact integral of this explicitly truncated coupling, integral alpha_s/(2*pi) d ln Q2."
function evolution_time(coupling::AsymptoticTwoLoopCoupling, Q2, Q02)
    L0 = coupling_log(coupling, Q02)
    coupling_log(coupling, Q2)
    beta0, beta1 = 11 - 2coupling.nf/3, 102 - 38coupling.nf/3
    # log1p retains precision for nearby scales and avoids an overflowing ratio.
    delta = abs(Q2-Q02) < Q02/2 ? log1p((Q2-Q02)/Q02) : log(Q2)-log(Q02)
    ratio = delta/L0
    logratio = log1p(ratio)
    correction = (logratio - (log(L0)+1)*ratio)/(L0*(1+ratio))
    tau = (2/beta0) * (logratio + beta1/beta0^2 * correction)
    isfinite(tau) || error("nonfinite integrated coupling")
    tau
end
