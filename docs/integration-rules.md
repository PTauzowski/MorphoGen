# MorphoGen — Integration Rules

**Status:** working rules for the branch integration described in
[branch-integration-plan.md](branch-integration-plan.md) · **Last updated:** 2026-09-10

These rules govern the merge of `Vibrations` and `CAS_Arm` into `develop`. They exist because an
integration of this size fails in a predictable way: the merge becomes an excuse to improve
things, scope grows without bound, and nobody can tell afterwards whether a changed result came
from the merge or from an improvement made along the way.

---

## Rule 1 — Merging is the only purpose

> No development of new features, analyses, finite elements, or similar is allowed. Only what
> already exists may be moved. **One exception: tests.**

### Permitted

- Moving a file from one branch to `develop`.
- Moving a file from the library into a study folder (`examples/<Study>/HelperFunctions/`).
- Editing a call site so existing code keeps working under the API chosen in Phase 2.
- Deleting code that is dead — see the ruling below.
- **Fixing errors** — code that cannot run, or that computes the wrong answer.
  One condition: **in its own commit, never bundled with a merge resolution.** See the ruling below.
- Writing tests. This is the stated exception; see Rule 3.

### Forbidden

- New analysis types, element types, shape functions, optimisers, or solvers.
- New capability added to an existing class because it "would be easy while we're in here".
- Performance work, style refactors, renaming that the merge does not force.
- Reorganising directories beyond moving studies into their own folders.

### Rulings on edge cases

**Deleting dead code is permitted.** Rule 1 says "only already existing should be moved", which
does not obviously cover deletion. The intent of the rule is to prevent *growth*, not to force
the carrying of corpses. Six files are referenced by nothing on the branches that carry them
(`Frame2D`, `StressIntensityTopologyOptimizationVol2`, `FreeVibrationsTopologyOpt`,
`StressIntensityTopologyOptimizationMultiLoad`, `SolidElasticElem-old`). Merging three
divergent copies of a file nobody calls is work with no product.

> **`Frame2D` is the exception, corrected 2026-09-11.** It is dead on both *feature* branches,
> which is where the original census was taken, and live on `develop` — `FrameOnElasticGround.m`
> and `FrameOnElasticGroundModel.m` both construct it, and develop repaired its `Ke` bug in
> Phase 1. **It is kept.** Deleting it on the two prep branches was still right: it turned a
> five-hunk add/add conflict into develop's copy surviving untouched, which is the outcome
> wanted. The lesson generalises — deadness is a property of a *tree*, and the merge target is
> the tree that decides. Deletion is permitted, subject to two conditions:

1. Deadness is proven by a whole-tree reference check, not by inspection.
2. The branch tip is tagged first, so the file stays reachable forever.

**Fixing errors is permitted — in its own commit.** An error is code that cannot run, or that
computes the wrong answer. Merging something known to be broken helps nobody, and a merge is
when these surface, so they get fixed.

The single condition exists because of how this codebase fails. There is no numerical baseline
yet, results are plausible-looking floating-point fields, and a fix bundled into a merge
resolution is **indistinguishable from a merge error** when the numbers later disagree — you
cannot bisect a commit that did two things. So:

> Fix it, commit it alone, say in the message what was wrong. Never inside a conflict resolution.

Where the fix is needed to make a merge step verifiable at all, do it *before* that step rather
than during it. Four such errors are known:

| Error | Fix | When |
|---|---|---|
| `Frame2D.computeStifnessMatrix` multiplies by `Ke`, never assigned | `Ke` → `Kl(:,:,k)` | Phase 1 — **done** *(moot if `Frame2D` is deleted as dead)* |
| `FEAnalysis.plotSupport`: `1:size(irots)` is not a scalar colon operand | `1:max(size(irots))` | Before capturing main's baseline; already fixed on develop |
| `FEAnalysis:17`: `for k=max(size(obj.felems))` is missing its `1:`, so the `eDofs` union runs on the last element group only | `for k=1:max(size(obj.felems))` | **Done** — `0674e86` |
| `FEAnalysis` stored the `felems` cell array in the caller's orientation and indexed it as a column, so a row-oriented multi-group model assembled only its first group | Normalise with `felems(:)` in the constructor | **Done** — `2969edc`, covered by `TestMultiGroupAssembly` |
| `plotSolid` was renamed to `plot` on both `PlaneElem` and `SolidElasticElem`, leaving 13 example scripts calling a method that no longer exists | `.plotSolid(` → `.plot(` | **Done** — found while capturing the numerical oracle |

The last two are invisible with a single element group and wrong with several, which is why
neither had been noticed — every example in the suite has one group. Both were found while
deciding D6 and are recorded in full there.

The second was not latent. Measured on a two-group model before the fix:

```
getTotalElemsNumber   4  (true total 8)
getElemIndices        1 range   (should be 2)
weighted assembly     numel(I)=512  numel(K)=256   -> 1 of 2 groups assembled
```

`CorbelModelMultiMat` builds `{fe1 fe2 fe3}` for its three materials and calls `solveWeighted`,
so it could not have run on `develop` at all.

The `plotSolid` case is the same shape as the `ndofs`→`eDofs` rename Phase 1 finished: a method
renamed on the element classes with its callers left behind. It went unnoticed for the same
reason — all 13 stale callers are topology-optimisation benchmarks and modular-structure
examples, and the smoke set contains neither. It surfaced only because the numerical oracle
needed `Beam99.m` and `Lshape.m` to actually run. It belongs on the "still broken" list in
[baseline-2026-09-10.md](baseline-2026-09-10.md) beside the five multi-block fixtures — and
unlike those, it is now fixed.

#### Not an error — deferred under Rule 2

**Gauss-point result layout differs between element families.** `PlaneStressElem` has its
`permute` call commented out and indexes `results.gp.stress` as `(component, elem, ip)`;
`SolidElasticElem` permutes to `(elem, ip, component)`. So `gp.stress(1,:,:)` means
*sxx everywhere* for plane elements and *element 1, all points, all components* for solids.

Nothing here computes a wrong answer — each family is self-consistent and its own callers index
it correctly — so this is not an error under the rule above. Unifying it would change a public
result layout and break every caller of one family or the other, which makes it an architectural
change the merge does not force: **out of scope under Rule 2**, and worth an issue afterwards.

The code is identical on `develop`, `Vibrations` and `CAS_Arm`, so it is pre-existing rather
than a merge artefact. Found by the constant-stress patch test on its first run, which is a fair
advertisement for writing that test early.

**Transitional wrappers — decided, see D1.** The plan proposes thin deprecated wrappers so the
11 legacy `globalMatrixAggregation` call sites keep working during Phase 5, removed in Phase 6.
Strictly, a wrapper is new code. Two readings, both defensible:

- *Permitted as scaffolding* — it is not a feature, it is a migration aid, and it must be gone
  by the end of Phase 6.
- *Forbidden* — migrate all 11 call sites in one commit and never introduce the wrapper.

**Decided: permitted as scaffolding** — see D1 below. The obligation that comes with it is
removal in Phase 6; a wrapper that outlives the merge has become a feature by accident.

---

## Rule 2 — Preserve the existing architecture

> Only architectural changes explained by the merging process are allowed.

The test for any architectural change is a single question: **would this change be unnecessary
if the branches had never diverged?** If yes, it is out of scope. If the change exists only to
reconcile two things that both already exist, it is in scope.

### In scope by this test

| Change | Why the merge forces it |
|---|---|
| Choosing one assembly API from the three that exist | Three dialects cannot coexist. Choosing between existing implementations invents nothing. |
| Finishing develop's `ndofs`→`eDofs` / `sf`→`shapeFn` rename | develop does not run. A merge onto a base that cannot execute is unverifiable — you cannot attribute a failure to your resolution or to the pre-existing breakage. |
| Folding `LinearElasticityWeighted` into `LinearElasticity` | Vibrations already did this; CAS_Arm did not. One of the two must win. |
| Choosing one test convention | Both branches have `tests/`, in incompatible styles. See Rule 3. |
| Moving studies into per-study folders | Shared example paths are the direct cause of 24 of the 73 conflict hunks. |

### Out of scope

Redesigning `DOFManager` beyond what is needed to make it work with the merged element classes;
introducing an interface or base class that neither branch has; changing the `Mesh` API because
the merged version is awkward; adopting a package/namespace layout.

### The architectural question — decided (D2)

develop's DOF refactor is **half-finished**: `FiniteElement` declares `eDofs`/`shapeFn`, but 11
files still read `obj.ndofs` and 8 read `obj.sf`. There are two ways to satisfy Rule 2, and the
rules as written do not choose between them:

- **(a) Finish it.** Complete the rename on develop, keep `DOFManager` and the `FEModel` that
  owns DOF numbering. This is the plan's Phase 1. It preserves the newest architecture and is
  the only path that supports frames and solids on one mesh — which is what the arm study needs.
- **(b) Revert it.** Drop the refactor, take the working element classes from `CAS_Arm`, and
  land the studies on the older architecture that demonstrably runs.

**(a) was chosen** — see [D2](#d2--develops-dof-refactor-is-finished-not-reverted-2026-09-10).
The capability is real, `CAS_Arm` has already partly adopted `eDofs` (11 files), and reverting
would discard 39 commits of deliberate work. (b) remained a legitimate reading of "preserve the
existing architecture", which is why the choice is recorded rather than assumed. The known risk
stands: because develop does not run, `DOFManager` has never been exercised against real models,
so Phase 5 may surface genuine design gaps in `initDOFs` rather than mere merge conflicts.

---

## Rule 3 — A `tests/` folder, with real tests

> `tests/` should be created, and tests of finite elements, shape functions, solvers etc.
> should be implemented.

### What already exists

Neither `develop` nor `main` has a `tests/` folder. Both feature branches do — three files
each, entirely disjoint, so all six move to develop without conflict:

| Branch | Files | Style |
|---|---|---|
| `CAS_Arm` | `test_Frame3DSectionProps.m`, `test_frameBasedSolver.m`, `test_rotationAwareElementMapping.m` | Plain script; manual `PASS`/`FAIL` printing; `addpath` inline |
| `Vibrations` | `TestPillarFEAPExportChemistry.m`, `TestPillarVerticalPlanes.m`, `TestWoodburyReanalysis.m` | `matlab.unittest.TestCase` classdef; `PathFixture`; deterministic `rng` |

### Convention: adopt the `matlab.unittest` style

Vibrations' style wins, and this is a Rule 2 decision the merge forces. It is the MATLAB
standard framework: it runs under `runtests('tests')`, reports properly, supports fixtures and
setup/teardown, and returns a non-zero status suitable for CI later. The three CAS_Arm scripts
are **converted**, not rewritten — same assertions, same tolerances, wrapped in a `TestCase`.
Naming: `Test<Subject>.m`, one class per subject.

### Coverage to implement

Ordered by value. The first two groups are cheap and catch the errors a merge actually causes;
the third is the real validation.

**Shape functions** — `math/ShapeFunction*.m` (L2, L3, L4, L8, L9, L16, L27, T3, T4, T6, Frame2D)

- Partition of unity: `sum(N(p)) == 1` at arbitrary interior points.
- Kronecker delta: `N_i(node_j) == δ_ij`.
- Gradient sum: `sum(dN(p)) == 0`.
- Integration: measure of the reference element via each integrator matches the known value.

**Finite elements** — `elements/*.m`

- Symmetry: `K == K'`.
- Rigid-body modes: `K * u == 0` for translations and rotations; rank deficiency equals the
  number of rigid-body modes (3 in 2D, 6 in 3D). *This is the single most valuable element test
  — it catches almost every transformation-matrix error, which is exactly what a merge of
  divergent `Frame3D` implementations risks.*
- `Frame3DSectionProps`: EA, EIy, EIz, GJ, GAy, GAz each enter the local stiffness matrix
  independently *(already covered by the CAS_Arm test — convert it)*.
- **Constant-stress patch test** — specified separately below; it is the one test that
  validates inter-element load application.

#### The constant-stress patch test

The most valuable test in the suite, because it is the only one that exercises the edge/face
load integral — the path where a wrong Jacobian or a lumped-instead-of-consistent nodal force
produces results that look plausible and are wrong.

**Setup.** A rectangular beam, uniform pressure applied to one edge (2D) or face (3D), shear
(roller) support on the opposite edge or face — restrained normal to it, free to slide — plus
one node pinned in the remaining directions to remove rigid-body motion.

**Assertion.** At **every Gauss point of every element**, the stress component in the loading
direction equals the applied pressure, and the other components are zero:

```
σxx == p    at all (element, integration point)
σyy == 0,  σxy == 0                    (2D)
σyy == σzz == σxy == σyz == σzx == 0   (3D)
```

**Tolerance is machine precision, not engineering tolerance.** Any element able to represent a
constant strain state must reproduce this exactly, on any mesh, at any resolution — that is the
defining property of a patch test. Use a relative tolerance around `1e-10`; a result that is
merely close is a failure, not a pass.

Assert on Gauss-point values, not nodal ones: `fe.results.gp.stress(:, elem, ip)` is already
populated by `computeResults`. Nodal results are extrapolated and averaged between elements,
which is exactly the smoothing that would mask a bad load integral.

**Parameterise over shape functions.** The existing example already sweeps `ShapeFunctionQ4`,
`Q9`, `Q16` and `T3`; all must pass identically.

**The fixtures already exist** — this is a conversion under Rule 1, not new development:

| Asset | Role |
|---|---|
| `examples/elasticity/planeProblems/ConstStressTest.m` | 2D driver — currently plots and prints, asserts nothing |
| `examples/elasticity/solidProblems/ConstStressSolidTest.m` | 3D driver |
| `examples/models/ConstPlaneStressModel.m` | 2D fixture: `elementLoadLineIntegral` on `x = 2l`, `ux` fixed on `x = 0`, `uy` pinned at the origin |
| `examples/models/ConstPlaneStressModelTriangular.m` | triangular-element variant |
| `examples/models/ConstStressSolidModel.m` | 3D fixture: `elementLoadSurfaceIntegral` on `x = l`, `ux` fixed on `x = 0`, `uy`/`uz` pinned at the origin |

The work is to turn the driver scripts into `TestCase` classes that assert the condition above,
keeping the models as fixtures.

**Useful side effect.** Both fixtures construct `LinearElasticityWeighted` and call
`solveWeighted`, so this test is a direct canary for the D3 API change: when that class folds
into `LinearElasticity`, the patch test is the first thing that should be made to pass again.

**Solvers** — `math/LinearEquationsSystem.m`, `LinearEquationsSystemTr2D.m`, `subsolv`, `mmasub`

- A small system with a known exact solution.
- Constrained-DOF handling: supports are actually enforced.
- Rotated supports: `LinearEquationsSystemTr2D` against a hand-computed 2-DOF case.
- Woodbury reanalysis *(already covered by the Vibrations test — keep)*.

**Analyses — closed-form validation.** These are the tests worth the most, because they check
physics rather than regression, and they exercise precisely the three capabilities being merged:

| Analysis | Check against |
|---|---|
| `LinearElasticity` | Cantilever tip deflection `PL³/3EI` |
| `LinearStability` | Euler buckling load of a pinned column, `π²EI/L²` |
| `LinearNaturalVibration` | First natural frequency of a simple beam |

### Scope bound

Rule 1 exempts tests from the freeze, which means test-writing is the one activity that can
grow without limit and stall the merge indefinitely. The bound:

**Tests must cover the 14 library files the merge actually touches. Coverage beyond that is
follow-up work, tracked separately, and does not gate the merge.**

Study code is explicitly **not** in scope. `examples/<Study>/` is the record of one experiment
and its correctness is evidenced by the paper it produced; testing it now would be auditing
published results, which is a different project. The asymmetry is deliberate and permanent:

| Code | Tests required |
|---|---|
| Study code in `examples/<Study>/` | No |
| Library code in `analysis/`, `elements/`, `design/`, `mesh/`, `math/` | Yes |

The reason is cost, not principle. These branches exist because a conference deadline made
forking the library cheaper than extending it properly. A test rule that ignores that pressure
will be ignored in turn at the next deadline. Requiring tests only for library code puts the
cost where the benefit is — many studies depend on the library, so a silent break there costs
more than one paper — and makes the library's test suite the deliberate price of admission,
paid when a capability is promoted, never under deadline.

### Relationship to the smoke harness

The plan's Phase 0 smoke harness and these unit tests answer different questions and both are
needed:

| | Question | Gate for |
|---|---|---|
| Smoke harness | Does the example still *run*? | Every phase |
| Unit tests (Rule 3) | Is the element matrix *right*? | Phases 5 and 6 |
| Numerical parity vs `main` | Do we get the *same answers* as the published release? | Phase 5 |

A smoke harness passing tells you nothing about correctness — a silently wrong stiffness matrix
still converges to something plausible. That gap is why Rule 3 matters.

---

## Applying the rules — quick reference

| Activity | Rule | Verdict |
|---|---|---|
| Move a study into `examples/<Study>/` | 1 | Allowed |
| Demote study-only helpers out of `design/` | 1 | Allowed |
| Delete the 7 dead files | 1 | Allowed, with proof + tag |
| Finish develop's `eDofs` rename | 2 | Allowed — merge-blocking |
| Fix the `Ke` bug | 1 | Allowed — merge-blocking |
| Fix any *other* error found en route | 1 | Allowed — in its own commit |
| Unify the `gp.stress` index order | 2 | Forbidden — not an error; breaks callers |
| Choose one assembly API | 2 | Allowed (D6) |
| Relocate the SIMP weighting from the element into the assembler | 2 | Forbidden (D6) — the merge forces a name, not a layering |
| Hoist `globalMatrixAggregationWeighted` to `FEAnalysis` | 1 | Allowed (D6) — moving an existing method up one class |
| Fix the `felems` orientation and missing-`1:` errors | 1 | Allowed — own commit, before Phase 5 |
| Adopt CAS_Arm's `matrixIndexCache` | 2 | Deferred (D6) — performance the merge does not force |
| Add transitional wrappers | 1 | Allowed (D1) — must be removed in Phase 6 |
| Restore `alphas` into the `Multi*` base | 1 | Allowed (D3) — moving existing code |
| Restore `prepareRHSVectors()` that Vibrations commented out | 1 | Allowed (D7) — an omission, not an API decision |
| Port main's 11 orphans (4 are design classes) | 1 | Allowed — moving existing code |
| Write element/solver tests | 3 | Required |
| Improve `Mesh` API ergonomics | 2 | Forbidden |
| Add a new analysis type | 1 | Forbidden |

---

## Proposed additional rules

Not yet adopted — offered for a decision.

**R4 — No new dependencies.** No new toolboxes, no external packages. A merge that changes what
the code needs to run is not only a merge.

**R5 — Every study keeps a runnable example.** A study is only successfully merged when at
least one of its examples runs on `develop`. This is what stops a study from being "merged" as
a directory of files nobody can execute.

**R6 — One resolution per commit.** Each conflict resolution, deletion, or study move is its own
commit with a message saying which rule permits it. With 43 library hunks and no CI, `git bisect`
is the only practical way to find which resolution broke a number.

**R7 — Numerical parity gates the merge.** Phase 5 does not complete until `Cantilever.m`,
`Beam99.m` and `Lshape.m` reproduce the reference values captured from `main` in Phase 0, within
a stated tolerance. Without this, the merge can ship a quiet numerical bug.

**R8 — No renaming the merge does not force.** Renames are nearly invisible in a diff but break
every caller. Beyond `ndofs`→`eDofs` and `sf`→`shapeFn`, which the merge forces, nothing is
renamed.

---

## Decisions taken

The three questions the rules left open have been decided.

### D1 — Transitional wrappers are permitted *(2026-09-10)*

Thin deprecated wrappers may be introduced so legacy call sites keep working across the API
change in Phase 2. They are scaffolding, not features. **They must be removed in Phase 6**, and
their removal is part of that phase's gate — a wrapper that survives the merge has become a
feature by accident.

### D2 — develop's DOF refactor is finished, not reverted *(2026-09-10)*

Complete the `ndofs`→`eDofs` and `sf`→`shapeFn` rename on develop, keep `DOFManager` and the
`FEModel` that owns DOF numbering. This is Phase 1 of the plan.

Consequence to plan for: `DOFManager` has never been exercised against real models, because
develop does not currently run. Phase 5 should be expected to surface genuine design gaps in
`initDOFs` — particularly around mixed frame/solid meshes, which is the case the arm study
needs and the reason the refactor exists. Budget for that as engineering, not as merge conflict.

### D3 — `StressIntensityMultiAv/MultiMax` stay in the library *(2026-09-10)*

**They were misclassified.** The consumer heuristic that flagged them as study-only mistook
"only one study has used it so far" for "specific to that study". Reading them settles it: both
extend the library base `StressIntensityTopologyOptimizationVol`, both are named for a generic
formulation (multi-load-case stress intensity, averaged vs worst-case), and neither contains
arm- or vibration-specific terminology. Duplicating them into two study folders would fork a
shared optimiser — precisely what this integration exists to undo.

**They are also not duplicates of each other.** CAS_Arm *refactored* them: it introduced a
common base `StressIntensityMultiTopologyOptimization` parameterised by an aggregation strategy,
reducing both subclasses to five-line constructors. Vibrations still carries the original pair
with duplicated bodies.

**Resolution — take CAS_Arm's structure, restore Vibrations' weighting:**

1. Adopt CAS_Arm's three files: the strategy-parameterised base plus two thin subclasses. It is
   the better organisation and carries newer features (`useParallel`, `maxdisplacement`,
   `plMaxStress`).
2. **Restore per-load-case weights (`alphas`) into that base**, defaulting to uniform.
   CAS_Arm's refactor dropped them; Vibrations depends on them.

The weighting is not a detail. The vibration study's call site is:

```matlab
StressIntensityMultiMaxTopologyOptimization(Rfilter, [analysisHarmonic1 analysisHarmonic2], ...
    [ (1-alphas(k)) alphas(k) ], cutTreshold, penal, volumeFractions, true);
```

`alphas(k)` is swept in a loop — the trade-off between two harmonic load cases **is** the
experiment. CAS_Arm's version takes no weights and aggregates uniformly. Its call sites pass six
arguments, Vibrations' pass seven, so a naive adoption fails loudly on arity; the danger is the
"fix" of dropping the argument, which leaves the sweep running with every weight equal and every
result quietly meaningless.

Defaulting `alphas` to uniform keeps CAS_Arm's six-argument call sites working untouched.

**Rule 1 compliance:** restoring `alphas` is *moving existing code* — the weighted aggregation
exists on `Vibrations` today and is in active use by a published study. Preserving it is the
merge's job, not new development. The optional parameter plumbing falls under D1.

**Rule 3 consequence:** this is a required test. A multi-load optimiser whose weighting silently
degrades to uniform produces plausible output and passes any smoke test. Assert that non-uniform
`alphas` yield a different intensity field than uniform ones.

### D4 — Article material is excluded from the merge *(2026-09-10)*

The article-writing prompt system has moved to its own repository and is developed there as a
distinct project. The articles themselves are mastered on Overleaf. Neither belongs in a
computational codebase, and neither is merged.

**Excluded — 35 files, dropped rather than resolved:**

| Material | Branch | Files |
|---|---|---|
| `ai/` — prompt system, workflows, agents, schemas | Vibrations | 27 |
| `ai/out/` — generated literature and style analyses | CAS_Arm | 3 |
| `docs/frameBasedSolver_method.tex` | CAS_Arm | 1 |
| `docs/frameBasedSolver_method.{aux,log,out,toc}` — LaTeX build artefacts | CAS_Arm | 4 |

**Verified safe:** no `.m` file on either branch references `ai/` or
`frameBasedSolver_method`. The exclusion breaks no code path.

**Consequences:**

1. Add `ai/`, `*.tex`, `*.bib` and the LaTeX build artefacts (`*.aux`, `*.log`, `*.out`,
   `*.toc`, `*.bbl`, `*.blg`, `*.synctex.gz`) to `.gitignore` in Phase 6, so the exclusion is
   enforced rather than remembered.
2. `docs/` on `develop` becomes unambiguously *documentation about the code* — this plan, these
   rules — and never paper text. The name collision with CAS_Arm's `docs/` resolves itself once
   the `.tex` material is dropped.
3. Where a study folder benefits from naming its publication, that goes in a short
   `examples/<Study>/README.md` — a citation and a link, not the manuscript. Only one study
   currently has a README; the rest can gain one during Phase 3 at negligible cost.

This does not apply to genuine code documentation: method notes explaining *what the
implementation does* remain welcome in `docs/`. The line is between documenting the software and
drafting a paper.

### D5 — main's orphaned design classes are preserved, not superseded *(2026-09-10)*

Four topology-optimisation classes exist only on `main`, added in `58535be` for the SoftwareX
paper after develop had diverged. All four are **ported, not dropped**:

| File | Disposition |
|---|---|
| `design/StressConstrainedTopologyOptimization.m` | **Preserve.** A distinct formulation, to be completed in future work. It is *not* superseded by the `StressIntensity*` family despite the similar name — stress-constrained optimisation and stress-intensity ESO are different methods. |
| `design/FatigueConstrainedTopologyOptimization.m` | Preserve — no successor anywhere downstream. |
| `design/ReliabilityConstrainedTopologyOptimization.m` | Preserve — no successor anywhere downstream. |
| `design/TopologyOptimization99.m` | Preserve — the classic 99-line reference implementation. |

`analysis/ElastoPlasticAnalysis.m` is the one exception: deliberately deleted on develop's line
in `a6b112a` (2023-11-13) and superseded by `ElastoPlasticity.m`. Safe to drop.

This matters because an incomplete formulation looks exactly like dead code to the reference
check that Rule 1 requires before deletion. `StressConstrainedTopologyOptimization` has no
callers today and would be deleted by a mechanical sweep. It is exempt: **absence of callers is
not evidence of deadness for a formulation still under development.**

### D6 — Assembly entry point: keep develop's two-method shape, hoisted *(2026-09-10)*

**This reverses the plan's Phase 2 recommendation.** §8 proposed adopting Vibrations'
`assemblyGlobalMatrix(fname, x, is_const)` on the grounds that `is_const` subsumes the
unweighted case. Reading the three implementations rather than their signatures overturns that.

**Adopted:** CAS_Arm's shape — `globalMatrixAggregation(fname)` and
`globalMatrixAggregationWeighted(fname, x)`, both on `FEAnalysis`.

This *is* develop's design. CAS_Arm changed two things about it, both of which are moves rather
than inventions: it hoisted `globalMatrixAggregationWeighted` from `LinearElasticityWeighted`
up to the base class, and it replaced grow-in-loop concatenation with a preallocated cell plus
`vertcat`. Moving an existing method up one class is Rule 1 "moving what already exists".

#### Why not Vibrations' entry point

**1. It is a layering change, not a rename.** develop's elements apply the density themselves:

```matlab
% elements/PlaneElem.m — computeStifnessMatrix(obj, nodes, varargin)
if ( nargin == 3 ), x = varargin{1}; else, x = ones(nelems,1); end
...
K(:,:,k) = x(k)*Ke;
```

Vibrations moved that multiply up into the assembler (`Ke = reshape(x,1,1,ne) .* Ke`) and left
the elements unweighted. On Vibrations, `PlaneElem` and `SolidElasticElem` **no longer define
`computeStifnessMatrix` at all** — the base `FiniteElement.computeStifnessMatrix(nodes, el_idx)`
does, and its third argument means an element index, not a density. Adopting the signature means
adopting the layering, which rewrites every element class in `elements/`.

Rule 2's test settles it: *would this change be unnecessary if the branches had never diverged?*
Choosing a name between three dialects — forced by the merge. Relocating the SIMP weighting
between two layers — not forced. The first is in scope, the second is not.

**2. Its weighting silently no-ops on multi-group models.** The multiply is guarded:

```matlab
if (numel(x)==ne)          % ne = this element group's count
    Ke = reshape(x,1,1,ne) .* Ke;
end
```

When `x` is the global density vector and `ne` is one group's element count, the guard fails,
the weighting is **skipped with no error**, and the solve proceeds against an unweighted
stiffness matrix. develop and CAS_Arm slice per group — `x(ei{k})` via `getElemIndices()` —
which is correct by construction.

The configuration that triggers this is the multi-element-group model, which is precisely what
develop's DOF refactor exists to enable and what the five outstanding multi-block fixtures
exercise. Importing a silent-wrong-answer path into the one capability the merge is being done
for is not a trade worth making.

#### Migrating Vibrations' call sites

Vibrations' 8 `assemblyGlobalMatrix` call sites get a D1 wrapper, removed in Phase 6:

```matlab
function K = assemblyGlobalMatrix(obj, fname, x, is_const)   % deprecated, Phase 6
    if is_const, fname = 'computeStifnessMatrixConst'; end
    K = obj.globalMatrixAggregationWeighted(fname, x);
end
```

#### Coupled change — the concatenation axis moves with it

`globalMatrixIndices` and the aggregator must use the **same** concatenation axis, because
`sparse(I,J,K)` consumes all three as `I(:), J(:), V(:)` and the linear order must agree:

| Branch | `globalMatrixIndices` | Aggregator |
|---|---|---|
| `develop` | horizontal — `[ I reshape(Ie',[],1) ]` | horizontal — `[ K elemK ]` |
| `Vibrations` | vertical — `[ I; reshape(...) ]` | vertical — `[ K; Ke(:) ]` |
| `CAS_Arm` | vertical | vertical — `vertcat(parts{:})` |

Each branch is internally consistent, which is why all three run. A resolution that takes the
aggregator from one branch and `globalMatrixIndices` from another produces a global matrix whose
entries land at the wrong `(I,J)` — and **a single element group masks it completely**, so every
smoke example would still pass.

> **Adopt CAS_Arm's vertical form for both methods, in one commit. Never split them across two.**

CAS_Arm also caches the index arrays (`matrixIndexCache`). That is a performance change the
merge does not force — **defer it**, and take it only if a profile later justifies it.

#### Two Rule 1 errors found while deciding this

Both are in the multi-element-group path, both are invisible with a single group, and both get
their own commit before the Phase 5 merge touches these methods:

| Error | Location | Effect |
|---|---|---|
| `for k=max(size(obj.felems))` — missing `1:` | `analysis/FEAnalysis.m:17` | The `eDofs` union runs on the **last group only**. With one group `max(size)==1`, so it is accidentally correct today. |
| `felems` orientation is assumed three different ways | `size(...,1)` at lines 31, 36, 45; `size(...,2)` at 266; `max(size(...))` elsewhere | `getElemIndices` returns one entry instead of *n* when the cell array's orientation disagrees, so `x(ei{k})` weights the wrong elements. |

The second matters directly to this decision: the per-group slicing that makes CAS_Arm's
aggregator correct depends on `getElemIndices` being right, and it is not yet.

### D7 — Density weighting folds into `LinearElasticity`, and the port runs three ways *(2026-09-10)*

**Adopted:** Vibrations' consolidation — one `LinearElasticity` carrying an `isConst` flag,
`LinearElasticityWeighted` deleted. The plan's §8 recommendation stands.

What does **not** stand is its framing of the cost — "port CAS_Arm's 2026 stress-constrained
additions first". The consolidated class has to absorb work from all three branches, and one of
the three contributions is a *removal* that must not be carried over.

**From CAS_Arm's `LinearElasticityWeighted` (2026-05)** — Vibrations' consolidated class has no
equivalent of any of it:

- the assembled-stiffness cache — `cachedWeightedX`, `cachedWeightedKvals`,
  `cachedWeightedFunction`, the `retainStiffness` flag, `hasCachedWeightedStiffness`,
  `clearWeightedStiffnessCache`
- `solveAdjointWithLoad(xPenal, P_adj_fem)` — adjoint sensitivity analysis reusing that cache
- `weightedStiffnessFunction()`, `saveMatrices(filename)`

The cache and the adjoint solve are one feature: `solveAdjointWithLoad` exists to avoid
reassembling `K` for the adjoint system. Porting either alone loses the point.

**From Vibrations:** `selfLoadFactor` and the `loadElementsSelfWeight(x, 0.1)` path in `solve`.

**From develop and CAS_Arm — restore `prepareRHSVectors()`.** Vibrations' `LinearElasticity.solve`
has it commented out:

```matlab
% bj.prepareRHSVectors();          % Vibrations — typo'd out, not removed deliberately
```

develop and CAS_Arm both call it, Vibrations still *defines* it, and six other Vibrations
analyses still call it — so this is local to `LinearElasticity`, not a considered API change.
Taking Vibrations' method body wholesale imports the omission silently.

It is not cosmetic. `prepareRHSVectors` calls `setCurrentLoadToRightHandSideVectors` and then
zeroes `Pnodal`, so it is the step that moves accumulated nodal load into the RHS. The plan's
§6 claim that the multi-RHS load API is "present and compatible on all three" is true of the
*methods* and false of this *call path*: on Vibrations' `LinearElasticity` it is never invoked.

#### Signature reconciliation

Three arities and two different return values meet here:

| Branch | Signature | Returns |
|---|---|---|
| `develop` | `LinearElasticity.solve()` | `[qn, K]` — nodal displacements **and** the global matrix |
| `develop` | `LinearElasticityWeighted.solveWeighted(x)` | `qfem` |
| `CAS_Arm` | `LinearElasticityWeighted.solveWeighted(x, retainStiffness)` | `qfem` |
| `Vibrations` | `LinearElasticity.solve(x)` | `qfem` |

Consolidated signature: **`solve(x, retainStiffness)`, both optional.** Absent `x` means
unweighted, preserving develop's `solve()` call sites; `retainStiffness` defaults false,
preserving CAS_Arm's six-argument sites and Vibrations' `solve(x)` alike.

The return value is the sharp edge: develop's `solve()` yields `[qn, K]` and the others yield
`qfem`. Callers taking two outputs must keep working, so the consolidated `solve` returns
`[qfem, K]` with `K` computed only when a second output is requested (`nargout > 1`).

#### Scope note

**44 files on `develop` reference `LinearElasticityWeighted`**, independently of the 51 stale
callers the plan attributes to Vibrations. Those 44 are not broken today — the class exists on
develop — but D7 makes them all Phase 6 migration work. Phase 6's caller migration is therefore
substantially larger than §4.3's "51 files" implies, and should be planned as the bulk of that
phase rather than one bullet in it.

**Rule 3 consequence.** The patch test already constructs `LinearElasticityWeighted` and calls
`solveWeighted`, so it is the canary for this change — it is the first thing that must be made
to pass again once the fold lands, exactly as §"The constant-stress patch test" anticipated.

---

### D8 — The shape-function family takes Vibrations' geometry-letter names *(2026-09-11)*

**Adopted:** `Q`/`H`. develop and CAS_Arm migrate to Vibrations' names before Phase 4.

Found during the Phase 3 triage, 2026-09-11; not identified by the survey, and it would have
blocked Phase 4 silently.

`Vibrations` renamed five shape-function classes, adopting a scheme in which the letter names
the element geometry:

| develop · CAS_Arm | Vibrations | Geometry |
|---|---|---|
| `ShapeFunctionL4` | `ShapeFunctionQ4` | quadrilateral |
| `ShapeFunctionL9` | `ShapeFunctionQ9` | quadrilateral |
| `ShapeFunctionL16` | `ShapeFunctionQ16` | quadrilateral |
| `ShapeFunctionL8` | `ShapeFunctionH8` | hexahedron |
| `ShapeFunctionL27` | `ShapeFunctionH27` | hexahedron |

`ShapeFunctionL2`, `L3`, `L4l`, `T3`, `T4`, `T6` are untouched on every branch, so after the
rename `L` means *line*, `Q` *quad*, `H` *hex*, `T` *triangle/tet* — internally consistent,
which the old scheme was not: `L4` was a quad and `L8` a hex.

**The rename is pure.** `ShapeFunctionL4.m` and `ShapeFunctionQ4.m` are byte-identical apart
from the classdef line and the constructor name; `L8`/`H8` differ additionally by one internal
`facesf` reference. No shape function, derivative or integration rule changes. **R7 numerical
parity is therefore unaffected by the choice** — this is a naming decision with no physics in
it.

**Why it blocks.** Call-site census over every `.m` file:

| | L-family | Q/H-family |
|---|---:|---:|
| `develop` | **61 files** | 0 |
| `prep/arm` | **73 files** | 0 |
| `prep/vibrations` | 0 | **84 files** |

The danger is that this conflicts *nowhere*. The five old files are deletes on the Vibrations
side and the five new ones are adds; git merges both cleanly. The Phase 4 merge
(`prep/vibrations ← prep/arm`) therefore produces, without a single conflict marker, a tree in
which 73 arm files construct classes that do not exist — and Phase 5 adds develop's 61. This is
§10's "clean textual merge, semantically wrong" risk realised in a file set §10 did not name.

**The options.**

| | Files to migrate | Notes |
|---|---:|---|
| Adopt `Q`/`H` (Vibrations' scheme) | 134 (61 develop + 73 arm) | Keeps the better scheme. Larger sweep, but it runs on the two branches we are already editing. |
| Keep `L` (develop's scheme) | 84 (vibrations) | Smaller sweep, but reverts a deliberate improvement and leaves `L4`-the-quad next to `L2`-the-line permanently. |

Either way the migration is a mechanical, verifiable substitution — five whole-word names, no
semantic review — and it must land **before** Phase 4, on whichever branches lose the vote, so
that the merge sees one vocabulary.

Note the interaction with **R8**: R8 forbids renames *the merge does not force*. This one is
forced — the branches already disagree, so there is no option that renames nothing.

**Carried out on develop**, 67 files and 87 substitutions in one commit: the five class files
renamed, every `.m` call site swept, the five `resources/project/*.xml` records repointed, and
the one `README.md` example updated. The sweep is a single word-boundary-guarded substitution
rule with no per-file judgement, which is what makes it reviewable at that size — and the guard
is load-bearing: `ShapeFunctionL4l` is a *different, live* class that a naive `ShapeFunctionL4`
substitution would have silently renamed to `ShapeFunctionQ4l`, and it is referenced by
`ShapeFunctionQ16`.

Local variable names (`sfL4 = ShapeFunctionQ4()` and friends) were deliberately **not** touched.
They are R8 territory — the merge does not force them — and changing them would have turned a
mechanical sweep into 53 files of cosmetic edits that have to be read. The mismatch is visible
and harmless; it can be tidied after the merge, or never.

---

### D9 — The `analysis/DOFManager*.m` family is abandoned scaffolding *(open — needs a decision)*

**Found 2026-09-11, running the only example that exercises mixed element classes.**

The plan calls develop "the branch that holds the DOF-manager rewrite" and treats that as the
reason develop is the merge target. The measurement is more specific than that, and splits in
two.

**The capability is real and it works.** DOF numbering lives in `FEModel.initDOFs()` — 80
documented lines that union each node's DOFs across every element class touching it. Run against
`FrameOnElasticGround`, the one mixed frame/plane model in the repository, it is correct:

```
mesh nodes           : 4336
DOFs per node (uniq) : [2;3]      ux,uy on plane-only nodes; ux,uy,fiz where the frame attaches
model DOFs total     : 8688
distinct DOF types   : ux,uy,fiz
```

This is the capability §10 of the plan lists as "unproven at scale". It is now proven for the
nonuniform case, and the example runs in 0.4 s, so it has been added to the smoke set.

**The four `DOFManager` classes are something else: an abandoned attempt to extract that logic
into a strategy hierarchy.** All four were last touched on 2025-01-10, nine months before
develop's tip. No library file references any of them. Worse, they cannot work as written:

| File | Defect |
|---|---|
| `DOFManagerNonuniform.m` | Constructor computes `nodalDOFS`, `DOFsInds`, `globalDOFs`, `elemDOFs` as locals and assigns **none** of its four properties, so `obj.nodesToDofs` is `[]` when `getIndices` indexes it |
| `DOFManagerUniform.m` | `getIndices` reads `obj.eDofs` and `obj.elems`; neither it nor its base `DOFManager` declares either |
| `DOFManagerNodalUniform.m` | `getIndices` reads `obj.elems`; the class declares only `nDOFs` |
| `DOFManager.m` | Abstract base. Sound, and unused. |

Three of the four would throw on first call. The only caller anywhere was a two-line probe in
`FrameOnElasticGround.m`, which is why that example failed — removed in `5b98521` under Rule 1,
leaving a working model.

**The question.** These are dead by exactly the test Rule 1 requires, but D5 cautions that
absence of callers is not evidence of deadness for work still under development, and these sit
on the merge target rather than on a branch being retired.

| | |
|---|---|
| **Delete** | Honest: the logic they were extracting works where it already is, and a reader who finds four DOF managers next to a working `initDOFs` will reasonably assume the managers are the live path. Git keeps them. |
| **Keep with a note** | If the extraction is meant to resume, a header comment saying so costs nothing and D5 is precedent for it. |

Either way, **the Phase 5 gate must be reworded**: "DOF-manager path exercised by at least one
mixed frame/solid model" names the wrong thing. The path that needs exercising is
`FEModel.initDOFs`, and the model that exercises it is `FrameOnElasticGround`.

---

## When the rules do not answer

For anything the rules do not cover: the default is **no**. The merge is finished when the
branches are merged, not when the codebase is better.
