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
- Writing tests. This is the stated exception; see Rule 3.

### Forbidden

- New analysis types, element types, shape functions, optimisers, or solvers.
- New capability added to an existing class because it "would be easy while we're in here".
- Performance work, style refactors, renaming that the merge does not force.
- Reorganising directories beyond moving studies into their own folders.

### Rulings on edge cases

**Deleting dead code is permitted.** Rule 1 says "only already existing should be moved", which
does not obviously cover deletion. The intent of the rule is to prevent *growth*, not to force
the carrying of corpses. Seven files are referenced by nothing (`Frame2D`,
`StressIntensityTopologyOptimizationVol2`, `FreeVibrationsTopologyOpt`,
`StressIntensityTopologyOptimizationMultiLoad`, `SolidElasticElem-old`), and `Frame2D` is a
five-hunk conflict in every merge combination. Merging three divergent copies of a file nobody
calls is work with no product. Deletion is permitted, subject to two conditions:

1. Deadness is proven by a whole-tree reference check, not by inspection.
2. The branch tip is tagged first, so the file stays reachable forever.

**Bug fixes are permitted only where they block the merge.** Two are known:

| Fix | Why it is permitted |
|---|---|
| `Frame2D.computeStifnessMatrix`: `Ke` → `Kl(:,:,k)` | `Ke` is never assigned; the function cannot run. *(Moot if `Frame2D` is deleted as dead.)* |
| `FEAnalysis.plotSupport`: `1:size(irots)` → `1:max(size(irots))` | Blocks every example on `main`, which is the numerical oracle. Already fixed on develop. |

Any other bug found during the merge is **recorded, not fixed** — open an issue and move on.
A bug fix silently bundled into a merge commit is indistinguishable from a merge error when the
numbers later disagree.

**Transitional wrappers — open question.** The plan proposes thin deprecated wrappers so the
11 legacy `globalMatrixAggregation` call sites keep working during Phase 5, removed in Phase 6.
Strictly, a wrapper is new code. Two readings, both defensible:

- *Permitted as scaffolding* — it is not a feature, it is a migration aid, and it must be gone
  by the end of Phase 6.
- *Forbidden* — migrate all 11 call sites in one commit and never introduce the wrapper.

**This needs a decision before Phase 2.** The second is cleaner and the call-site count is
small enough to make it practical.

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

### The unsettled architectural question

develop's DOF refactor is **half-finished**: `FiniteElement` declares `eDofs`/`shapeFn`, but 11
files still read `obj.ndofs` and 8 read `obj.sf`. There are two ways to satisfy Rule 2, and the
rules as written do not choose between them:

- **(a) Finish it.** Complete the rename on develop, keep `DOFManager` and the `FEModel` that
  owns DOF numbering. This is the plan's Phase 1. It preserves the newest architecture and is
  the only path that supports frames and solids on one mesh — which is what the arm study needs.
- **(b) Revert it.** Drop the refactor, take the working element classes from `CAS_Arm`, and
  land the studies on the older architecture that demonstrably runs.

**(a) is recommended** — the capability is real, `CAS_Arm` has already partly adopted `eDofs`
(11 files), and reverting would discard 39 commits of deliberate work. But (b) is a legitimate
reading of "preserve the existing architecture", and choosing (a) should be a recorded decision
rather than an assumption. Note that (a) carries a known risk: because develop does not run,
`DOFManager` has never been exercised against real models, so Phase 5 may surface genuine design
gaps in `initDOFs` rather than mere merge conflicts.

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
- Patch test: a constant-strain field is reproduced exactly.
- `Frame3DSectionProps`: EA, EIy, EIz, GJ, GAy, GAz each enter the local stiffness matrix
  independently *(already covered by the CAS_Arm test — convert it)*.

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
| Fix any *other* bug found en route | 1 | **Record, do not fix** |
| Choose one assembly API | 2 | Allowed |
| Add transitional wrappers | 1 | **Undecided — settle before Phase 2** |
| Port main's 6 orphan examples | 1 | Allowed — moving existing code |
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

---

## When the rules do not answer

For anything the rules do not cover: the default is **no**. The merge is finished when the
branches are merged, not when the codebase is better.
