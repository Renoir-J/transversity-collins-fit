# Fixed numerical cache integrity

`KPSY-NUMERICAL-CACHE-001` covers the fixed data used by retained HOPPET plans,
detached evolution/matching rows, initial samples, prepared polarized contexts,
and derived leptonic observables. It changes no physical equation or numerical
acceptance tolerance. Production and physical-extraction flags remain false.

Julia's ordinary array `hash` can sample elements, so it cannot establish that
every coefficient or source cell is unchanged. Fixed-cache snapshots use SHA256
over all dense primitive-array bytes, numeric types and array dimensions, framed
sequences and sorted string-keyed metadata. An explicit native-endian marker is
included. Accessible mutable arrays are checked completely before application.

Names for the 14 built-in Julia primitive numeric types are stored as immutable
strings at module load. Other types keep runtime `string(T)` formatting. This
reduces repeated type-display work without changing a snapshot byte or sampling
array content. Parent-generated synthetic scalar/matrix digests cover all 14
types and a composite metadata example; a custom-type display regression checks
the fallback. These digest fixtures explicitly declare their native-endian format.

Retained row matrices are large and fixed. Preparation copies each distinct
matrix into owned storage, computes its complete snapshot once and exposes a
read-only `AbstractMatrix` interface. Ordinary writes, writes through views and
mutating array operations fail. `copy` and slicing produce independent mutable
arrays. Multiplication uses the same dense Float64 BLAS or generic high-precision
arithmetic as the parent implementation, including nested ForwardDiff columns.
Repeated queries can share the same owned row object; payload accounting retains
the parent meaning. Block indices, scale/grid metadata and all accessible source,
CSS, quadrature, flavor and observable arrays remain part of validation.

This is a trusted-program in-memory integrity contract. Private `getfield`/raw
storage access, unsafe memory operations and forged internal objects are outside
the supported API. It is not an authenticated serialization format or a
production-cache admission certificate. File/source/build receipts remain
separate and mandatory. Changing a constructor argument after ownership transfer
does not mutate the cache; modifying an independent copy requires constructing a
new plan.

Validation must include interior coefficient/source mutations, write-through
views, changed block mappings and source metadata, restoration checks, physical
input comparisons and parameter derivatives. Preserve failed parent probes.
Measure construction-time copying and snapshots separately from warm numerical
application; matrix timing alone cannot certify production G2PERF, fit throughput
or interval coverage.

The locked-source producer exercises SIDIS and SIA caches, ordinary row writes and
views, block mappings, scale/grid metadata, initial bases, provider/CSS/flavor
arrays and projected momentum/observable weights. Every changed case must be
rejected with the expected integrity error and the original observable must be
restored exactly. The producer uses explicit reference-scale policy and a
diagnostic angular factor; its software receipt does not certify experimental
acceptance or numerical convergence.

```bash
julia --project=numerics numerics/bin/check_numerical_cache.jl \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json \
  runs/numerical-cache/integrity.json
```
