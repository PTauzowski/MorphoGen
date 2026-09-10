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
_(filled in as runs complete — see the tables below)_
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
