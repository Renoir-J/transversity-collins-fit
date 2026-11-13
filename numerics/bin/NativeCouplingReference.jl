module NativeCouplingReference

using CollinsFit, LinearAlgebra, SHA
export native_request, parse_native_response, run_native_raw, native_batch,
       rg_reference, log_rule, scaled_error, comparison

function native_inputs(alpha0,Q02,Q2,nf,nloop,dt)
    a=CollinsFit.fixed_scalar(alpha0,"reference alpha0")
    q0=CollinsFit.fixed_scalar(Q02,"reference Q02")
    step=CollinsFit.fixed_scalar(dt,"reference dt")
    q=[CollinsFit.fixed_scalar(v,"reference Q2") for v in Q2]
    nf isa Integer && !(nf isa Bool) && 3<=nf<=6 || error("reference fixed nf 3..6 required")
    nloop isa Integer && !(nloop isa Bool) && nloop in (1,2) || error("one/two loops required")
    0<a<=.4 && 2.4<=q0<=1e8 && .0001<=step<=.2 || error("reference settings outside contracted domain")
    1<=length(q)<=20000 && all(v->q0<=v<=1e8,q) || error("forward scales in contracted support required")
    (alpha0=a,Q02=q0,Q2=q,nf=Int(nf),nloop=Int(nloop),dt=step)
end

function native_request(alpha0,Q02,Q2,nf,nloop,dt)
    s=native_inputs(alpha0,Q02,Q2,nf,nloop,dt)
    "CF-COUPLING-1\n$(length(s.Q2)) $(s.nf) $(s.nloop)\n$(s.alpha0) $(s.Q02) $(s.dt)\n"*
        join(string.(s.Q2),'\n')*"\n"
end

function parse_native_response(text,n,nf,nloop)
    lines=split(chomp(text),'\n')
    length(lines)==n+2 && lines[1]=="CF-COUPLING-1" || error("incomplete/wrong native response protocol")
    parse.(Int,split(lines[2]))==[n,nf,nloop] || error("native response header mismatch")
    values=parse.(Float64,strip.(lines[3:end]))
    all(v->isfinite(v) && v>0,values) || error("nonfinite/nonpositive native alpha")
    values
end

"One complete native comparator batch, with every protocol and diagnostic byte retained."
function run_native_raw(executable,request,directory,index)
    mkpath(directory)
    prefix=joinpath(directory,"batch-"*lpad(string(index),5,'0'))
    input,output,diagnostics=prefix*".request.txt",prefix*".response.txt",prefix*".stderr.txt"
    any(ispath,(input,output,diagnostics)) && error("refuse to replace existing native batch evidence")
    write(input,request)
    code=open(input,"r") do incoming
        open(output,"w") do outgoing
            open(diagnostics,"w") do errors
                process=run(pipeline(ignorestatus(Cmd([executable]));stdin=incoming,stdout=outgoing,stderr=errors))
                process.exitcode
            end
        end
    end
    evidence=Dict("request"=>input,"response"=>output,"stderr"=>diagnostics,"exit_code"=>code,
        "sha256"=>Dict(p=>bytes2hex(sha256(read(p))) for p in (input,output,diagnostics)),
        "precision_errors"=>count("TOO HIGH ACCURACY REQUIRED",read(diagnostics,String)))
    (text=read(output,String),stderr=read(diagnostics,String),evidence=evidence)
end

function native_batch(executable,alpha0,Q02,Q2,nf,nloop,dt,directory,index)
    request=native_request(alpha0,Q02,Q2,nf,nloop,dt)
    raw=run_native_raw(executable,request,directory,index)
    raw.evidence["exit_code"]==0 || error("native coupling comparator failed; evidence retained")
    raw.evidence["precision_errors"]==0 || error("native coupling precision error; evidence retained")
    values=parse_native_response(raw.text,length(Q2),nf,nloop)
    (values=values,evidence=raw.evidence)
end

"Independent high-precision solution of the fixed-nf one/two-loop beta equation."
function rg_reference(alpha0,Q2,Q02,nf,nloop;bits=256)
    s=native_inputs(alpha0,Q02,[Q2],nf,nloop,.05)
    setprecision(bits) do
        a0,q,q0=BigFloat(s.alpha0),BigFloat(s.Q2[1]),BigFloat(s.Q02)
        beta0=BigFloat(11)-BigFloat(2)*nf/3
        beta1=BigFloat(102)-BigFloat(38)*nf/3
        B0=beta0/(4big(pi));c=nloop==1 ? BigFloat(0) : beta1/(4big(pi)*beta0)
        F(a)=inv(a)+c*log(a/(1+c*a))
        target=F(a0)+B0*log(q/q0)
        alpha=if q==q0
            a0
        elseif nloop==1
            a0/(1+B0*a0*log(q/q0))
        else
            low,high=BigFloat(0),a0
            for _ in 1:bits+16
                middle=(low+high)/2
                (middle==low || middle==high) && break
                if F(middle)>target;low=middle;else;high=middle;end
            end
            (low+high)/2
        end
        tau=q==q0 ? BigFloat(0) : (2/beta0)*log(a0*(1+c*alpha)/(alpha*(1+c*a0)))
        (alpha=Float64(alpha),tau=Float64(tau),implicit_residual=Float64(abs(F(alpha)-target)))
    end
end

"Gauss integration in ln(Q2), with alpha/(2*pi) weights; no missing factor of two."
function log_rule(Q02,Q2,n)
    q0,q=CollinsFit.fixed_scalar(Q02,"Q02"),CollinsFit.fixed_scalar(Q2,"Q2")
    0<q0<=q && n isa Integer && !(n isa Bool) && n>=2 || error("forward scales and at least two Gauss nodes required")
    q==q0 && return (nodes=Float64[],weights=Float64[])
    spectral=eigen(SymTridiagonal(zeros(n),[j/sqrt(4j^2-1) for j in 1:n-1]))
    span=log(q/q0)
    nodes=q0.*exp.((spectral.values.+1).*(span/2))
    weights=span .* spectral.vectors[1,:].^2 ./ (2pi)
    (nodes=nodes,weights=weights)
end

scaled_error(a,b,floor)=maximum(abs.(a.-b)./max.(abs.(b),floor))

function comparison(name,actual,reference,tolerance;floor=1e-12,absolute=false)
    a,b=Float64.(vec(collect(actual))),Float64.(vec(collect(reference)))
    !isempty(a) && length(a)==length(b) || error("matching nonempty comparison operands required")
    isfinite(tolerance) && tolerance>=0 || error("finite nonnegative tolerance required")
    difference=absolute ? maximum(abs.(a-b)) : scaled_error(a,b,floor)
    Dict("name"=>name,"actual"=>a,"reference"=>b,"error"=>difference,"tolerance"=>tolerance,
        "metric"=>absolute ? "absolute" : "scaled","relative_floor"=>floor,
        "passed"=>all(isfinite,a) && all(isfinite,b) && isfinite(difference) && difference<=tolerance)
end

end
