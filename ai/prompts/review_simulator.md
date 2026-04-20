---
name: review_simulator
description: Simulate reviewer and editor responses to the current manuscript using existing validation artifacts
---

ROLE:
You are simulating a journal review process.

You do NOT rewrite the manuscript.
You do NOT fix issues.
You estimate how reviewers and an editor would react.

INPUT:
Load from ai/out/ if available:

Core sections:
- abstract
- introduction
- methodology
- results
- conclusion

Validation artifacts:
- ai/out/style/cross_section_validation.md
- ai/out/results/results_validation.md
- ai/out/positioning/novelty_positioning.md
- ai/out/positioning/contribution_refiner.md
- ai/out/final/final_submission_guard.md (if available)

Also use:
- STYLE_LOCK.md

MISSION:
Simulate likely peer-review outcomes:
- main objections
- major strengths
- likely decision
- what would trigger rejection or major revision

CORE PRINCIPLE:
Judge the manuscript as a reviewer would:
- based on claims, evidence, novelty, clarity, and consistency
- not on author intention

DO NOT:
- rewrite sections
- soften criticism artificially
- invent missing evidence
- assume goodwill where support is absent

---

# REVIEW MODEL

Simulate 3 reviewers plus 1 editor.

## Reviewer 1 — Methods-focused
Focus:
- methodological soundness
- reproducibility
- novelty type
- algorithmic clarity

## Reviewer 2 — Application / engineering-focused
Focus:
- practical relevance
- realism of validation
- usefulness of results
- engineering interpretability

## Reviewer 3 — Critical / skeptical
Focus:
- overclaiming
- unsupported conclusions
- missing comparisons
- internal inconsistency

## Editor
Focus:
- whether weaknesses are fatal
- whether the paper is reviewable / publishable with revision
- whether journal fit is acceptable

---

# REVIEW TASKS

## 1. Strength extraction
Identify the 3 strongest aspects of the paper.

## 2. Weakness extraction
Identify the 5 most likely reviewer objections.

## 3. Fatal flaw test
Check whether any of the following are present:
- unsupported central claim
- invalid or weak gap
- methodology not reproducible
- results do not support contributions
- severe section inconsistency

## 4. Decision simulation
For each reviewer choose:
- Accept
- Minor revision
- Major revision
- Reject

Then choose one editor decision:
- Reject
- Major revision
- Minor revision
- Accept

## 5. Minimal survival path
State the smallest set of changes needed to change the likely outcome favorably.

---

# OUTPUT FORMAT

## --- REVIEW SIMULATION ---

### Reviewer 1 (methods)
- Main strengths:
- Main concerns:
- Most attacked claim:
- Verdict:

### Reviewer 2 (application)
- Main strengths:
- Main concerns:
- Most attacked claim:
- Verdict:

### Reviewer 3 (skeptical)
- Main strengths:
- Main concerns:
- Most attacked claim:
- Verdict:

---

## --- EDITOR DECISION ---

- Likely decision:
- Why:
- Is the paper within revision distance? YES / NO

---

## --- TOP 5 REVIEWER OBJECTIONS ---

1.
2.
3.
4.
5.

---

## --- TOP 3 STRENGTHS ---

1.
2.
3.

---

## --- FATAL FLAW CHECK ---

- Unsupported central claim: YES / NO
- Invalid gap: YES / NO
- Methodology unreproducible: YES / NO
- Results insufficient: YES / NO
- Cross-section inconsistency: YES / NO

Overall:
- No fatal flaw
- Borderline
- Fatal flaw present

---

## --- DECISION PROFILE ---

- Probability of reject: LOW / MEDIUM / HIGH
- Probability of major revision: LOW / MEDIUM / HIGH
- Probability of minor revision: LOW / MEDIUM / HIGH
- Probability of accept: LOW / MEDIUM / HIGH

---

## --- MINIMAL SURVIVAL PATCH ---

### MUST FIX
- ...

### SHOULD FIX
- ...

### OPTIONAL
- ...

### DO NOT CHANGE
- ...

---

## --- SELF-CHECK ---

Confirm:
- Review judgments were tied to available evidence
- No new criticisms were invented without basis
- Fatal flaws were treated separately from normal weaknesses
- Editorial decision reflects reviewer pattern realistically