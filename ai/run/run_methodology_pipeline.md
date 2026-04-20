---
name: run_methodology_pipeline
description: Execute the full Methodology pipeline with code-grounded writing, technical audit, and reproducibility validation
---

ROLE:
Execute the full Methodology pipeline autonomously.

INPUT:
- Use the currently focused document as the manuscript
- Use the full project codebase (via VSCode access)
- Use references/ only if needed for attribution

ALWAYS LOAD:
- ai/STYLE_LOCK.md
- ai/WORKFLOW.md

MISSION:
Produce a Methodology section that is:
- technically correct
- fully consistent with implementation
- reproducible
- reviewer-proof

CORE PRINCIPLE:
The Methodology must be verifiable against code and sufficient for reimplementation.

DO NOT:
- invent algorithmic steps not present in code
- rely on abstract descriptions when concrete implementation exists
- skip audit stages
- rewrite entire sections unnecessarily
- assume “standard” steps without stating them

---

# PIPELINE

## STEP 0 — PRECONDITION CHECK

Confirm:
- manuscript is available in focused document
- codebase is accessible
- main implementation entry points exist

If code is not accessible:
- STOP
- report: “Methodology pipeline requires code access”

---

## STEP 1 — METHOD DISCOVERY (CODE-FIRST)

Inspect code:

Identify:
- main solver / runner (e.g., entry scripts)
- core algorithm functions
- data flow
- parameters and defaults
- special mechanisms (aggregation, filtering, constraints)

Output (internal):
- method structure map
- code trace map

DO NOT write Methodology yet.

---

## STEP 2 — WRITE METHODOLOGY

Run:
- ai/prompts/write_methodology.md

Goal:
- produce Methodology grounded in actual implementation
- include:
  - formulation
  - algorithm flow
  - key mechanisms
  - assumptions

Output:
- Methodology draft
- structure map
- code traceability map
- mismatch flags

---

## STEP 3 — STRUCTURE VALIDATION (NEW)

Run:

ai/prompts/structure_validator.md

Save:
- ai/out/style/structure_validation_methodology.md

---

### BLOCK CONDITIONS

If ANY of the following:

- CRITICAL structure issues exist (section boundary violation, micro-heading)
- fragmentation = HIGH (≥ 3 headed blocks under 5 lines)
- hierarchy depth > allowed (> 2 levels below `\section{}`)
- CRITICAL standard-component over-detail (standard component has its own `\subsection{}`, or > 10-line description of unmodified standard component)

→ status = FAILED
→ STOP pipeline
→ Write block reason to: ai/out/style/structure_validation_methodology.md

Rationale:
Structural and boundary problems must be corrected before semantic audit.
The write_methodology.md self-check should have prevented these; if they appear, the draft must be revised.

## STEP 3b — COMPANION SECTION CHECK

Goal:
Detect whether Methodology over-absorption has hollowed out the Experiments section.
This check runs even when the full-paper section_identity_validator is not available.

Procedure:

1. Locate the Experiments section (e.g., sec_experiments.tex or equivalent).
   If it does not exist → flag as WARNING: no Experiments section found.

2. Estimate content volume:
   - Thin: < 100 words of non-boilerplate content
   - Adequate: ≥ 100 words covering setup, baselines, and runs

3. Cross-reference with structure validation output (STEP 3):
   - Count how many subsections were flagged as boundary violations (ablation design, evaluation protocol, dataset splits, etc.)

4. Apply the hollowing rule:
   IF Experiments section is THIN
   AND structure validation found ≥ 1 boundary violation in Methodology
   → CRITICAL: Methodology absorption has hollowed the Experiments section
   → STOP pipeline
   → Write reason to: ai/out/style/structure_validation_methodology.md

   IF Experiments section is THIN but no boundary violations found:
   → WARNING only — note that Experiments may need expansion

Save:
- companion section status to: ai/out/style/structure_validation_methodology.md (append)

---

## STEP 4 — METHODOLOGY AUDIT

Run:
- ai/prompts/methodology_audit.md

Goal:
- verify:
  - logical correctness
  - completeness
  - code consistency
  - reproducibility

Output:
- audit report
- issue list
- reproducibility gaps
- reviewer risks

---

## STEP 5 — CONTRIBUTION ALIGNMENT CHECK

Run:
- ai/prompts/contribution_refiner.md (analysis mode only)

Goal:
- ensure:
  - contributions are actually implemented in Methodology
  - no “paper-only” contributions exist

Output:
- mismatches between claims and method

---

## STEP 6 — ISSUE FILTER

Run:
- ai/prompts/issue_filter.prompt.md

Input:
- issues from methodology_audit.md
- mismatches from contribution_refiner

Goal:
- classify:
  - ESSENTIAL
  - OPTIONAL
  - HARMFUL TO CHANGE

Output:
- safe change set

---

## STEP 7 — SAFE PATCH

Run:
- ai/prompts/safe_patch.prompt.md

Goal:
- apply only critical fixes:
  - missing definitions
  - missing steps
  - inconsistencies
  - reproducibility gaps

Rules:
- do not rewrite entire Methodology
- preserve structure
- maintain alignment with code

Output:
- revised Methodology
- applied changes
- skipped changes

---

## STEP 8 — REPRODUCIBILITY CHECK

Re-evaluate:

Ask:
- could an expert reimplement this method from text alone?

If NO:
- identify missing elements
- mark as NOT READY

If PARTIAL:
- list missing elements

---

## STEP 9 — REVIEWER SIMULATION

Run:
- ai/prompts/reviewer_simulation.prompt.md

Focus on:
- technical clarity
- missing steps
- suspicious assumptions
- unclear algorithm flow

Output:
- reviewer objections
- high-risk sentences

---

## STEP 10 — FINAL DECISION

Decide:

IF:
- no critical reproducibility gaps
- code and text are aligned
- logic is consistent

→ READY

IF:
- only minor gaps remain

→ BORDERLINE

IF:
- missing steps
- code mismatch
- unclear logic

→ NOT READY

---

# DECISION RULES

- Prefer adding missing steps over rewriting
- Prefer explicit definitions over implicit assumptions
- Prefer correctness over brevity
- If code contradicts text → code wins

---

# OUTPUT FORMAT

## --- FINAL METHODOLOGY ---
<final section>

---

## --- PIPELINE SUMMARY ---

- Code inspected: YES / NO
- Method discovered: YES / NO
- Methodology written: YES / NO
- Audit performed: PASS / PARTIAL / FAIL
- Safe patch applied: YES / NO
- Reproducibility check: PASS / PARTIAL / FAIL
- Reviewer simulation: YES / NO

---

## --- STATUS ---

Choose one:
- READY
- BORDERLINE
- NOT READY

---

## --- CODE ALIGNMENT STATUS ---

- Fully aligned:
- Partially aligned:
- Misaligned components:

---

## --- APPLIED CHANGES ---

- ...

---

## --- REMAINING GAPS ---

- ...

---

## --- TOP REVIEWER ATTACKS ---

1.
2.
3.

---

## --- REPRODUCIBILITY SUMMARY ---

- Can method be reimplemented: YES / PARTIAL / NO
- Missing elements (if any):

---

## --- NEXT ACTION ---

Choose one:
- proceed to Results pipeline
- fix reproducibility gaps
- align contributions with method
- investigate code–text inconsistencies

---

## --- FAILURE HANDLING ---

If:
- code is missing
- key algorithm cannot be identified

Then:
- STOP
- report missing elements explicitly
- do NOT generate Methodology

---

## --- SELF-CHECK ---

Confirm:
- Methodology reflects actual implementation
- All major algorithmic steps are described
- No invented steps were introduced
- Reproducibility was explicitly tested
- Critical issues were fixed with minimal edits