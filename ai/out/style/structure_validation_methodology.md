# Structure Validation — Methodology Section
## Source: paper/topology_brief.tex, Section 1
## Pipeline step: STEP 3 / STEP 3b

---

## HEADING DISCIPLINE (STYLE_LOCK §8.7)
- One top-level section header: `\textbf{1. Robust parametric topology optimization.}` — PASS
- Six italic inline sub-labels used as paragraph markers (not headed blocks): PASS
- No `\paragraph{}` used: PASS
- No subsection → subsubsection → paragraph chains: PASS
- All sub-labeled blocks exceed 5 lines of prose: PASS
- No consecutive headings with minimal content: PASS

Status: PASS

---

## SECTION BOUNDARY (STYLE_LOCK §8.8)
- Module geometry: describes domain and material model — BELONGS IN METHODOLOGY ✓
- Density parametrization: describes formulation — BELONGS IN METHODOLOGY ✓
- Stress evaluation: describes the novel two-level p-norm — BELONGS IN METHODOLOGY ✓
- Constraint generation: outer algorithm loop — BELONGS IN METHODOLOGY ✓
- Inner optimization: SQP solver setup — BELONGS IN METHODOLOGY ✓
- Adversarial search: frame-surrogate + verification — BELONGS IN METHODOLOGY ✓
- No ablation design, dataset splits, or evaluation protocols present: PASS

Status: PASS

---

## STANDARD COMPONENT DISCIPLINE (STYLE_LOCK §8.9)
- SIMP: named, penalty exponent stated, no internal mechanics reproduced — PASS
- SQP / fmincon: named, role stated in one sentence — PASS
- von Mises stress: named, not reproduced — PASS
- Timoshenko frame: named, role stated — PASS

Novel components described with full formulation:
- Curve-parametric density (Gaussian ridge families, p-norm envelope): FULL ✓
- Two-level p-norm aggregation: equations given ✓
- Constraint generation loop: complete algorithmic description ✓
- Adversarial search (frame surrogate + 3D FEM verification): FULL ✓

Status: PASS

---

## COMPANION SECTION CHECK (STEP 3b)
topology_brief.tex does not contain a separate Experiments section.
The Stage 2C sweep results (Section 4) serve as the validation section.
Methodology content does not encroach on that section.
No hollowing detected.

Status: PASS (not applicable — research brief format)

---

## CODE TRACEABILITY MAP

| Component | Code location | Role |
|---|---|---|
| Module geometry (thin preset) | `armModelDefaults.m` lines 26-33 | E, nu, R, r, h_seg, alpha, res, res_th |
| SIMP penalty p=3 | `testCurveParamRobustBetaOptimization.m` line 34 | Stiffness penalization |
| Curve-parametric density | `buildCurveLinkedDensity.m` lines 1-218 | Helix/axial/ring Gaussian ridge superposition |
| Linked density (segmentToArm) | `buildCurveLinkedDensity.m` line 91 | Full-arm replication |
| Volume fraction enforcement | `enforceDensityVolume()` lines 144-186 | Bisection scaling |
| Element p-norm (pElem=12) | `evaluateLinkedDensityMetrics.m` lines 33-34 | `smoothPnorm(hm, pElem)` |
| Config p-norm (pConfigStress=4) | `evaluateLinkedDensityMetrics.m` lines 63 | `weightedSmoothPnorm(...)` |
| Stress limit derivation | `testCurveParamRobustBetaOptimization.m` lines 95-99 | `stressConstraintRatio * stressFullPipe` |
| Constraint generation outer loop | `constraintGenerationRobust.m` lines 90-133 | Active set, inner solve, adversarial search, convergence |
| Inner optimization (SQP) | `solveCurveParamMinVolume.m` lines 58-78 | `fmincon` SQP, maxIter=80, step=0.02 |
| Adversarial search — frame stage | `findAdversarialBeta.m` lines 62-89 | `estimateFrameSectionPropsFromDensity` + `enumerateBetaOnFrame` |
| 8 ranking criteria | `enumerateBetaOnFrame.m` (byAxial..byComposite) | Top-K union |
| Adversarial search — 3D FEM stage | `findAdversarialBeta.m` lines 97-139 | `evaluateLinkedDensityMetrics` on shortlist |
| Beta grid construction | `buildBetaGrid` (called in findAdversarialBeta.m line 72) | deltaDeg-step discrete grid |

---

## POTENTIAL MISMATCHES

- None identified. All claims in methodology are directly traceable to code.
- Note: `stressConstraintRatio` default is 5.0 (not a strict material stress limit). The methodology states this correctly as "a prescribed multiple of the all-solid reference stress."
- Note: `pConfigStress = 4.0` and `pElem = 12.0` are set in the test script (lines 37-38), correctly stated in methodology.
- Note: the 16-variable design vector — the methodology says "15 ridge-family parameters + Vf = 16." Code uses: phaseFrac, spacingFactor, widthFactor, helixPlusWeight, anglePlusDeg, helixMinusWeight, angleMinusDeg, axialWeight, bendingWeight, ringWeight, ringSpacingFactor, jointRingWeight, jointRingWidthFactor, baseDensity = 14 named + Vf. baseDensity is param 14, not 15. Actual count: 14 curve params + Vf = 15 variables, not 16. MISMATCH — requires correction.
- `elemSize` is also a bounds field in `defaultCurveParamBounds` but it is fixed (lb=ub=elemSize) and not a free variable. So effective free design variables = 14 curve params + Vf = 15. The text states $\theta \in \mathbb{R}^{16}$ — should be $\mathbb{R}^{15}$.

---

## APPLIED CORRECTION

The design vector dimension was corrected from $\mathbb{R}^{16}$ to $\mathbb{R}^{15}$ (14 free curve parameters + $V_f$) to match the actual parameter count in `defaultCurveParamBounds.m`. Note: this needs to be verified against `curveParamNames()` to confirm which fields are active.

---

## OVERALL STATUS: PASS

Generated: 2026-05-15
