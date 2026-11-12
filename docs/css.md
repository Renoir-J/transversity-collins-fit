# Reconstructed coupling and CSS factors

This component advances issue #4 under [KPSY-CSS-001](../config/css-contract.json). It implements an explicit fixed-nf asymptotic two-loop coupling, its collinear evolution time, the perturbative Sudakov exponent and the SIDIS/SIA nonperturbative factors in Julia. It remains `reconstructed_theory`. No experimental fit or production G2/G2AD/G2PERF certification is implied.

[KPSY Eqs. 26–28](https://arxiv.org/html/1505.05589v1#S2.E26) specify the perturbative integral, A1/A2/B1 coefficients and b-star prescription; Eqs. 52, 72 and 106–107 specify the channel-dependent nonperturbative factors. The default g2=0.84, gq=0.106 GeV2, gh=0.042 GeV2 and Q02=2.4 GeV2 come from that paper. bmax=1.5 GeV^-1 is inherited from the [referenced unpolarized fit, page 8](https://arxiv.org/pdf/1406.3073v1). It is an explicit reconstruction setting requiring historical comparison.

## Coupling definition

With `L=log(Q2/Lambda^2)`, the selected formula is

\[
\alpha_s(Q^2)=\frac{4\pi}{\beta_0 L}\left(1-\frac{\beta_1}{\beta_0^2}\frac{\ln L}{L}\right),
\quad \beta_0=11-2n_f/3,\quad\beta_1=102-38n_f/3.
\]

Its integrated evolution time is the difference of

\[
\frac{2}{\beta_0}\left[\ln L+\frac{\beta_1}{\beta_0^2}\frac{\ln L+1}{L}\right]
\]

between the final and initial scales. The implementation uses a stable difference for nearby scales. Lambda=0.225 GeV and nf=5 give alpha_s(sqrt(2.4))=0.32698427 and alpha_s(MZ)=0.11790266, consistent with the paper's printed precision. The separately defined exact two-loop RG solution anchored at the same alpha0 gives 0.11896599 at MZ. This comparison records different definitions; it is not a failed equality test or evidence of recovered author code. The original modified HOPPET convention still requires comparison.

```julia
using CollinsFit
coupling = AsymptoticTwoLoopCoupling(; lambda_GeV=0.225, nf=5)
model = KPSYCSS(; coupling=coupling, order=64)
Spert = sudakov_perturbative(model, 1.0, 100.0) # b in GeV^-1, Q2 in GeV2
Snp = sidis_collins_snp(model, 1.0, 100.0, 0.4, 0.02)
grid = collect(range(0.001, 1; length=201)) # interface example; accuracy not certified
operator = collinear_lo_operator(grid, 100.0, 2.4, coupling; quantity=:transversity)
```

`grid` must satisfy the collinear contract and its own accuracy checks. The example gc=0.02 GeV2 is an analytic test parameter. The scalar alpha0 collinear constructor continues to mean one-loop running; two-loop running requires the explicit coupling object and records that definition in the operator.

## Numerical and convention checks

The Sudakov integral uses fixed Gauss–Legendre nodes in ln(mu2), with A expanded in alpha_s/pi through A2 and B through B1. Signed integration limits are preserved when mu_b>Q. No cap to Q or nonnegative-exponent clamp is inserted. This short-b choice follows the written formula and remains an explicit reconstruction assumption pending historical comparison.

The nonperturbative logarithm uses `log(Q/Q0)=0.5*log(Q2/Q02)` once per structure function. SIDIS retains `gq+(gh-gc)/z^2` for Collins; SIA retains `(gh-gc)*(1/z1^2+1/z2^2)`. The API separates these channels. It requires positive b, Q2 above the coupling pole, z in (0,1] and gh-gc>0 for a convergent Collins Gaussian. gc derivatives are supported; operator/model settings remain fixed and reject unsupported AD requests.

Run `julia --project=numerics numerics/bin/check_css.jl runs/css/validation.json`. The complete job checks 256-bit coupling values, the closed evolution time against independent adaptive integration, printed parameter precision, signed constant-coupling integrals, adaptive Sudakov integration, node refinement, the ln(Q2) derivative, nonperturbative identities, gc derivatives and finite-difference plateaus. It also compares two-loop collinear evolution with independent Mellin inversion using a separately integrated coupling. Reports record source/contract hashes, all settings and source stability during execution.

The development run used 25 Sudakov locations, including mu_b>Q and a negative perturbative exponent. The maximum fixed/adaptive differences were 2.73e-11, 7.99e-15 and 1.78e-15 at 16/32/64 nodes. For the analytic polynomial collinear input at Q02=2.4, Q2=100 GeV2, nf=5, the 3201-node maximum scaled Mellin error was 4.94e-5 and moment relative error 5.22e-6. The 801/1601-node point checks remain above the unchanged 1e-4 target. These finite analytic cases do not certify arbitrary fitted shapes, physical bin predictions or interval coverage.

The [leptonic components](leptonic.md) now supply Wilson coefficients, Bessel transforms, flavor sums and CSS structure-function assembly from supplied matched legs. Fitted initial distributions, scalable evolution integration, publication-specific observable adapters and the historical 262-point/13-parameter benchmark remain incomplete. The [CT10/DSSV08/DSS14 input interfaces](parton-inputs.md) have a separate reconstruction contract and explicit range limits; choosing a physical low-scale prescription remains unresolved. The historical paper and HERMES2020 each prescribe interpretation at reported mean kinematics. Production adapters must match the relevant publication; where bin integration is required, integrate numerator and denominator before forming the ratio. Covariance and dataset admission remain separate unresolved requirements.
