# Numerical oracle — `main` @ `69dee29`

Captured for the outstanding half of Phase 0 of [branch-integration-plan.md](branch-integration-plan.md).
The pass/fail baseline in [baseline-2026-09-10.md](baseline-2026-09-10.md) records **what runs**;
this records **what the numbers are**, so that Phase 5 can tell a merge that works from a merge
that quietly changed the physics. It is the evidence R7 gates on.

**Method.** Detached worktree of `origin/main`, run under `matlab -batch` (R2025b) with figures
suppressed. The one-line `plotSupport` fix (`for k=1:size(irots)` → `1:max(size(irots))`) was
applied **to the worktree only and never committed** — `main` is read-only under the write
policy in §3 of the plan. The benchmark scripts were copied to untracked `oracle_*.m` variants
that call an `oracleDump` helper after each `solve()`; no tracked file on `main` was modified
apart from that one line.

Each benchmark runs **two** optimisers in sequence against the same analysis:

1. `StressIntensityTopologyOptimizationVol` — ESO-style, terminates on its own criterion
2. `SIMP_MMA_TopologyOptimizationElasticCompliance` — MMA, terminates on `change < 0.001`
   with no iteration cap

`objF` is each optimiser's own objective, not a common quantity — the two columns are not
comparable to each other, only to a later re-run of the same optimiser.

---

## The headline: only one of the three benchmarks is a valid oracle

The plan assumed all three of `Cantilever.m`, `Beam99.m` and `Lshape.m` would yield reference
values. Diffing them against their `develop` counterparts first shows that two do not:

| Benchmark | `develop` counterpart | Comparable? |
|---|---|---|
| `Beam99.m` | `examples/topologyOpt/benchmarks/Beam99.m` | **Yes** — the only difference is the Mesh API migration (`mesh.elems` → the return value of `addRectMesh2D`). Same physics, same resolution. |
| `Cantilever.m` | same path | **Only at matched resolution** — `main` uses `res = 100`, `develop` uses `res = 40`. Different meshes are different problems. Captured at both. |
| `Lshape.m` | same path | **No** — `main`'s version clamps almost the whole domain. See below. |

### `Lshape` on `main` solves a degenerate problem

`main` and `develop` disagree about `Lshape`'s boundary conditions, and the disagreement is not
cosmetic:

```matlab
% main
fixedEdgeSelector  = Selector( @(x)( x(:,2) - 2*l ) );                       % signed distance
loadedEdgeSelector = Selector( @(x)( not( (abs(x(:,1)-2*l) < 1.0E-4) & ... ) ) );

% develop
fixedEdgeSelector  = Selector( @(x)( abs(x(:,2) - 2*l) < 0.001 ) );          % predicate
loadedEdgeSelector = Selector( @(x)( (abs(x(:,1)-2*l) < 1.0E-4) & ... ) );
```

The two branches also disagree about what `Selector` *means*:

| | `Selector.select` |
|---|---|
| `main` | `s = abs(fn(points)) > tolerance;` |
| `develop` | predicate used as-is if logical; otherwise `abs(v) < tolerance` |

`main`'s single rule is right for the **predicate** convention (`abs(true)=1 > tol`,
`abs(false)=0` is not) and inverts the **signed-distance** convention, selecting everything
*far from* the surface. `Lshape` is the one benchmark that passes a signed distance.

Measured on `main`'s own L-shape mesh, 6683 nodes:

```
main   Selector(@(x) x(:,2)-2*l)       selects 6642 / 6683 nodes  (99.4%)
develop Selector(@(x) abs(..)<0.001)   selects   41 / 6683 nodes  ( 0.6%)
true top edge (y == 2l)                            41 nodes
```

So `main`'s `Lshape` fixes `ux` and `uy` on 99.4% of the structure and applies its load
everywhere *except* the intended point. Its optimiser output is not a different answer to
`develop`'s — it is the answer to a different and degenerate question.

**Consequence for R7:** `Lshape` cannot gate the merge. `develop` already corrected it (both
selectors moved to the predicate form, which reads identically under either `Selector`), so
there is nothing on `main` to reproduce. Its numbers are recorded below for completeness and
explicitly marked as not a reference. Parity rests on `Beam99` and `Cantilever`.

> This is also worth noting as a hazard for Phase 5 generally: `develop`'s `Selector` accepts
> both conventions and `main`'s accepts one, so a call site merged from either feature branch
> can change meaning silently depending on which `Selector` it lands next to. `develop`'s
> version documents the trap in its own comments; keep that version.

---

## Reference values

<!-- ORACLE_TABLE_START -->
### `main` @ `69dee29`

| Benchmark | Optimiser | `objF` | Volume fraction | Elements | Iterations |
|---|---|---:|---:|---:|---:|
| `Beam99` (res 40) | `StressIntensityTopologyOptimizationVol` | 1918.96595194 | 0.399784573321 | 4800 | 96 |
| `Beam99` (res 40) | `SIMP_MMA_…ElasticCompliance` | 230.191384232 | 0.399999999988 | 4800 | 1017 |
| `Cantilever` (res 40) | `StressIntensityTopologyOptimizationVol` | 1275.94055647 | 0.398731423896 | 3200 | 70 |
| `Cantilever` (res 40) | `SIMP_MMA_…ElasticCompliance` | 76.7366325352 | 0.399999999808 | 3200 | 113 |
| `Lshape` (res 20) | both | *did not converge* | — | — | aborted at 2234 |

`Cantilever` at `res = 40` is the resolution-matched variant; `main`'s script ships `res = 100`.

**`Lshape` does not converge on `main`.** Its SIMP run was still reporting `Vrel=100.0` — no
material removed at all — after 2234 iterations, and was aborted. That is the degeneracy
described above showing up at runtime: with 99.4% of the nodes clamped there is nothing for the
optimiser to remove.

### The comparison against `develop` @ `ff92833`

Both benchmarks, both optimisers, compared field by field:

| Benchmark | Optimiser | Quantity | `main` | `develop` | |
|---|---|---|---:|---:|:--|
| `Cantilever` (res 40) | SIMP | `objF` | 76.7366325352 | 76.7366325352 | **identical** |
| | | volume fraction | 0.399999999808 | 0.399999999808 | **identical** |
| | | iterations | 113 | 113 | **identical** |
| `Beam99` (res 40) | SIMP | `objF` | 230.191384232 | 230.191384232 | **identical** |
| | | volume fraction | 0.399999999988 | 0.399999999988 | **identical** |
| | | iterations | 1017 | 1017 | **identical** |
| `Cantilever` (res 40) | ESO | `objF` | 1275.94055647 | 1270.06320313 | differs |
| | | volume fraction | 0.398731423896 | 0.396894750978 | differs |
| | | iterations | 70 | 311 | differs |
| `Beam99` (res 40) | ESO | `objF` | 1918.96595194 | 1915.94611678 | differs |
| | | volume fraction | 0.399784573321 | 0.399155440995 | differs |
| | | iterations | 96 | 422 | differs |

**Both SIMP compliance runs reproduce `main` exactly — every digit captured, and the same
iteration count, on two independent benchmarks.** That is the result worth having. The SIMP path
exercises assembly, the linear solver, the sensitivity computation and the DOF bookkeeping, so
it says that develop's DOF refactor, the Phase 1 repairs and the two `FEAnalysis` fixes have not
moved the physics by one ulp.

Two independent benchmarks agreeing to twelve figures also rules out coincidence: an assembly
error that happened to leave one problem's compliance unchanged would not leave a second one
unchanged as well.

The ESO runs differ **in the same direction on both benchmarks** — develop reaches a *lower*
objective at a slightly *lower* volume fraction, in four to five times as many iterations. That
is the expected signature of a smaller removal step finding a better optimum, not of a
perturbed stiffness matrix, which would have no reason to improve the objective consistently.

### The ESO difference is the optimiser, not the physics

`StressIntensityTopologyOptimizationVol.m` is byte-identical between the branches, so the
divergence is in its base class. `StressIntensityTopologyOptimization.updateDesign` computes the
element-removal threshold differently:

```matlab
% main — a quadratic in the volume ratio
obj.maxais = -0.139*v^2 + 0.2694*v - 0.081;

% develop — a sigmoid schedule, plus a per-iteration removal cap
max_elem_removal_factor = 0.025;
```

A cap of 2.5% of elements per iteration is what turns 70 iterations into 311: develop removes
material in many more, much smaller steps, and lands on a slightly different final design.
`TopologyOptimization` also changed base class (`ConstrainedOptimization` → `handle`) on
develop's line.

So `main`'s ESO numbers are **not** a parity target either. The ESO algorithm itself differs
between the branches, deliberately and on develop's side.

> **A hypothesis that was tested and rejected.** The obvious explanation was commit `2dae36f`,
> which fixed the per-element density reshape into `results` slot 18 — a field ESO reads and
> SIMP does not. Applying *only* that fix to a `main` worktree and re-running gives
> `objF = 1275.94055647`, volume fraction `0.398731423896`, 70 iterations — **bit-identical to
> `main` as-is**. The reshape is not the cause. Recorded because it is the answer most likely to
> be assumed by the next reader.

### What actually gates Phase 5

Of the six optimiser runs the plan nominated, exactly **one** is a valid numerical parity target:

| Target | Status |
|---|---|
| `Cantilever` res 40, SIMP | **The gate.** Already reproduces exactly. |
| `Beam99`, SIMP | Valid — comparable, `main` value recorded above |
| `Cantilever`/`Beam99`, ESO | Not comparable — the removal schedule differs by design |
| `Lshape`, either | Not comparable — `main`'s problem is degenerate and does not converge |

R7 should be read against that table rather than against "the three benchmarks".
<!-- ORACLE_TABLE_END -->

---

## How to compare in Phase 5

Re-run the same three benchmarks on the merged tree and diff against the table above.

- `Beam99` compares directly.
- `Cantilever` must be compared against the **`res = 40`** row, which matches `develop`'s copy
  of the script. The `res = 100` row records `main`'s published configuration.
- `Lshape` is **not** a parity target; assert only that it still runs.

**Tolerance.** These are iterative optimisations with a floating-point convergence test, not
closed-form results, so exact equality is the wrong bar. Compare:

| Quantity | Bar |
|---|---|
| Volume fraction | tight — the constraint is explicit, so `0.4` should be hit to ~1e-3 |
| Objective `objF` | ~1e-6 relative if the assembly is unchanged; a shift in the third significant figure means the stiffness matrix changed |
| Iteration count | indicative only — a one- or two-iteration difference is convergence noise, a large change is not |

The iteration count is the least reliable of the three and the most sensitive: MMA stops on
`max|x - x_old| < 0.001`, so a tiny perturbation near the threshold can move it. Treat a
changed objective with an unchanged volume fraction as the signal that matters.
