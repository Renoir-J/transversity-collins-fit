# Numerical jobs and data contracts

Python passes one versioned JSON job to `julia --project=numerics numerics/bin/run_job.jl job.json`. A job owns a complete validation, fit, replica collection or holdout analysis. Julia never calls Python inside an objective or quadrature node.

The finite-kernel executor currently implements **synthetic closure only**, under [`NUMERICAL-JOB-SCOPE-001`](../config/numerical-job-contract.json). Job, data and kernel must each be JSON objects with the explicit exact string `scope="synthetic_closure"`. Missing, non-string, production, candidate, reconstructed-theory and other unsupported scopes fail. The job's scope is checked before input IO or action dispatch. The loader checks the data scope before opening the kernel, then the kernel scope before constructing the context; all three scopes must match. A production request never runs a synthetic fallback, and relabeling fixture JSON cannot enable a production backend.

```json
{
  "scope": "synthetic_closure",
  "action": "benchmark",
  "data": "fixtures/closure/data.json",
  "kernel": "fixtures/closure/kernel.json",
  "output": "runs/example/benchmark.json"
}
```

Separate reconstructed-physics validation scripts retain their explicit component interfaces and scopes. Static production-data normalization/review does not admit those files to this executor. Production execution still requires a separately implemented and validated scientific admission path; a metadata boolean is insufficient. These JSON checks do not sandbox arbitrary trusted Julia code or turn low-level numerical primitives into an authorization system.

The initial numerical implementation accepts finite preintegrated linear/bilinear kernels. Each bin identifies the exact point ID, an integrated positive unpolarized denominator, an explicit convention sign, optional constant and linear terms, and bilinear coefficient terms. Parameter indices are 1-based. The numerator is summed first and divided by the denominator once. Parameter transformations (`identity`, `exp`, `tanh`) are differentiated in Julia. These transformations define the declared coefficient domain; they do not automatically impose a physical bound on a renormalized Collins twist-3 moment.

This interface represents fixed-basis kernel contractions. It does not construct a CSS or finite-R jet kernel. If basis shapes, alpha_s, Q0, scales or other kernel dependencies are fitted, their effect must enter the differentiable implementation or invalidate/rebuild kernels. The delivered synthetic fixture has four algebraic coefficients and no experimental interpretation.

Experimental source intake uses the separate `experimental_candidate` scope and is rejected by both the fit-data validator and Julia loader. Its errors are retained as reported components, without silently assembling a diagonal systematic covariance. Per-publication prediction prescriptions are part of the scientific contract: HERMES2020 uses structure-function asymmetries and an acceptance-systematic interpretation at reported mean kinematics. Belle and BESIII retain different angular-average definitions; BESIII exchanged z bins are a union of orientations. STAR supplement means do not supply missing bin edges. See [experimental-data.md](experimental-data.md).

Data JSON includes a version, `scope`, ordered unique point IDs, dataset provenance, process, hadron, bins, fraction-valued observations and an explicit convention. Covariance must be finite, symmetric and positive definite; validation never adds diagonal jitter. Shared additive Gaussian systematics use a separate `nuisance_design` B with unit-normal priors. Julia minimizes the profiled quadratic including nuisance penalties; it agrees with fixed effective covariance C+BBᵀ. This is a frequentist objective. Parameter-dependent covariance, non-Gaussian errors, multiplicative/t0 systematics and Bayesian normalization require separate implementations and validation.

`collinsflow normalize raw.xlsx reviewed-mapping.json output.json` (or CSV) requires a SHA256, explicit columns/selection, fraction/percent conversion, sign, hadron and convention. It does not infer missing depolarization factors, covariance, bin semantics or hadron composition. Production provenance includes publication/table/version/source URL/hash, covariance status and selection reason. Full covariance can be supplied in a reviewed JSON file. Nuisance entries are supplied in already normalized fraction units. Combining datasets requires an overlap review and an assembled block/joint covariance.

`collinsflow plan-update current.json incoming.json` removes named replacements and requires refitting frequentist replicas. Only unchanged-contract actual posterior samples are candidates for importance reweighting; this command plans that route but does not execute it. Unknown replacement targets, duplicate IDs and unresolved overlapping projections fail.

The fixed-operator Fortran ABI uses column-major Float64 arrays, explicit dimensions and an integer error status. `cf_apply_operator` batches applications of an already constructed operator; it does not construct HOPPET transversity evolution. Julia applies operators in generic arithmetic to retain validated AD. Upstream HOPPET is unchanged, and no speedup over it has been claimed.

Each stage emits an atomic receipt with source/configuration, input/output hashes, environment digest, seed, scope, code commit, architecture, policy, precision and thread counts. Consumers validate ancestors and output integrity. Synthetic and production scopes never mix. These content hashes detect stale or changed artifacts; they are not cryptographic signatures or an authorization system.

