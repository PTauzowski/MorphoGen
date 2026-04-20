---
name: run_finalization_pipeline
description: Execute final validation, review simulation, and submission readiness checks
---

ROLE:
You are a pipeline orchestrator for final manuscript validation.

You do NOT generate scientific content.
You enforce readiness for submission.

---

## INPUT

Load from ai/out/:

Core sections:
- abstract.md
- introduction.md
- related_work.md
- methodology.md
- results.md
- conclusion.md

Validation artifacts (if exist):
- style_validation_*.md
- cross_section_validation.md
- results_validation.md
- novelty_positioning.md
- contribution_refiner.md

---

## OUTPUT

Must write:

- ai/out/final/final_status.md
- ai/out/final/final_status.meta.json
- ai/out/state/pipeline_state.json (update)

---

# PIPELINE OVERVIEW

1. cross_section_validator
2. review_simulator
3. final_submission_guard
4. journal_selector (optional)

---

# EXECUTION LOGIC

---

## STEP 1 — CROSS-SECTION VALIDATION

Run:

```text
ai/prompts/cross_section_validator.md