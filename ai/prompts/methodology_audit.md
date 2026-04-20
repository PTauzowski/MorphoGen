---
name: methodology_audit
description: Deep technical audit of Methodology for correctness, completeness, and reproducibility
---

ROLE:
Act as a senior reviewer in computational engineering auditing the Methodology section.

INPUT:
- Use the currently focused Methodology section (manuscript text is the primary source — audit operates on text regardless of how it was produced)
- Use the project codebase (full access via VSCode)
- Use if available (optional):
  - ai/out/methodology/methodology_draft.md — generator output for traceability
  - ai/out/positioning/contribution_refiner.md — for contribution alignment check

If generator outputs are absent, proceed using manuscript text and codebase directly.
Do NOT block on missing generator outputs.

ALWAYS LOAD:
- ai/STYLE_LOCK.md

MISSION:
Evaluate whether the Methodology:
- is logically correct
- is complete
- is consistent with implementation
- is reproducible
- would survive expert reviewer scrutiny

DO NOT:
- rewrite the Methodology
- improve style globally
- assume missing steps are “standard”
- ignore inconsistencies

CORE PRINCIPLE:
If a technically competent reader cannot reproduce the method from the description, the Methodology fails.

---

# AUDIT DIMENSIONS

## 1. METHOD LOGIC

Check:

- Is the method internally consistent?
- Does each step follow logically from the previous?
- Are there hidden transitions or jumps?

Detect:
- unexplained steps
- circular logic
- missing links

---

## 2. COMPLETENESS

Check whether the Methodology fully specifies:

- problem definition
- variables and notation
- objective function
- constraints
- solution procedure
- stopping criteria

Detect:
- missing definitions
- implicit assumptions
- underspecified steps

---

## 3. CODE CONSISTENCY

Compare Methodology with implementation:

For each major step:

- described in text?
- implemented in code?
- consistent?

Detect:
- described but not implemented
- implemented but not described
- mismatches in logic

---

## 4. MATHEMATICAL VALIDITY

Check:

- are equations correct?
- are variables consistently defined?
- does math correspond to actual algorithm?

Detect:
- symbolic inconsistencies
- unused equations
- misleading formalism

---

## 5. ALGORITHMIC CLARITY

Check:

- is computational flow clear?
- is iteration defined?
- is convergence defined?

Detect:
- vague algorithm description
- missing loop logic
- unclear dependencies

---

## 6. SPECIAL MECHANISMS

Focus on:

- aggregation strategies
- filtering
- stabilization
- coupling

Check:
- clearly defined?
- justified?
- correctly implemented?

These are often:
👉 the real contribution

---

## 7. PARAMETER TRANSPARENCY

Check:

- are key parameters defined?
- are their roles explained?
- are defaults or ranges given?

Detect:
- hidden hyperparameters
- unexplained constants

---

## 8. REPRODUCIBILITY

Ask:

- could an expert reimplement this from text?
- what information is missing?

List:
- missing elements required for reproduction

---

## 9. CONSISTENCY WITH CONTRIBUTIONS

Check:

- does Methodology actually implement the claimed contributions?
- or are contributions only conceptual?

Detect:
- contribution not realized in method
- mismatch between claim and implementation

---

## 10. REVIEWER ATTACK SIMULATION

Ask:

“What would a strict reviewer question here?”

Examples:
- “This step is unclear”
- “How is X computed?”
- “Where is this defined?”
- “Is this standard or new?”

---

## 11. STRUCTURAL QUALITY (STYLE_LOCK §8.7)

Check:

- Count all `\paragraph{}` heads. Are any used for blocks < 5 lines? → MAJOR per instance
- Is there a subsection → subsubsection → paragraph chain? → MAJOR
- Are consecutive headings separated by < 3 lines? → CRITICAL
- Total heading density: more than 1 heading per 5–8 lines of prose? → MAJOR

Classify overall:
- CLEAN: depth ≤ 2 levels, no blocks under 5 lines under a heading
- FRAGMENTED: ≥ 3 headed blocks under 5 lines
- OVER-SEGMENTED: more than 3 levels active in any chain

Flag each violation with exact location.

---

## 12. SECTION BOUNDARY DISCIPLINE (STYLE_LOCK §8.8)

For each subsection, apply the boundary test:
"Does this describe HOW THE METHOD WORKS, or HOW THE EXPERIMENTS WERE RUN?"

Flag as CRITICAL BOUNDARY VIOLATION if the subsection describes:
- ablation study design (factors, seeds, baselines)
- oracle or threshold analysis setup
- dataset splits
- transfer or zero-shot experiment setup
- benchmark selection rationale

Evaluation metrics — apply this distinction:
- NOT a boundary violation: metric definition integral to the method (e.g., why a specific metric is used, how it is computed in a non-obvious way)
- IS a boundary violation: standalone §Evaluation Metrics subsection that only names and defines standard metrics without connecting to the method design

Each flagged subsection must be marked for removal or relocation to Experiments.

---

## 13. STANDARD COMPONENT DISCIPLINE (STYLE_LOCK §8.9)

For each described component, classify:
- NOVEL / MODIFIED → full description is appropriate
- STANDARD (unmodified, published) → name + cite + 1 sentence maximum

Flag as OVER-DETAILED (MAJOR) if a standard component has:
- more than ~3–4 sentences of description
- its own `\paragraph{}` or `\subsubsection{}` heading
- internal equations or layer-level mechanics reproduced

Escalate to CRITICAL if:
- the standard component has its own `\subsection{}`
- OR its description exceeds 10 lines for a fully standard, unmodified component

List each violation with: component name, location, line count, severity.

---

# OUTPUT FORMAT

## --- METHODOLOGY AUDIT ---

### LOGIC
<assessment>

### COMPLETENESS
<assessment>

### CODE CONSISTENCY
<assessment>

### MATHEMATICAL VALIDITY
<assessment>

### ALGORITHMIC CLARITY
<assessment>

### SPECIAL MECHANISMS
<assessment>

### PARAMETER TRANSPARENCY
<assessment>

### REPRODUCIBILITY
<assessment>

### CONTRIBUTION CONSISTENCY
<assessment>

### REVIEWER RISK SUMMARY
<short synthesis>

### STRUCTURAL QUALITY
- Classification: CLEAN / FRAGMENTED / OVER-SEGMENTED
- Paragraph heads used: X (list any under 5 lines)
- Max depth: X levels
- Violations: list with location and severity

### SECTION BOUNDARY
- Boundary violations: list each misplaced subsection with severity CRITICAL

### STANDARD COMPONENT DISCIPLINE
- Over-detailed components: list each with severity (MAJOR / CRITICAL)

---

## --- CODE ALIGNMENT MAP ---

For each major component:

- Component:
- Described in text: YES / NO
- Implemented in code: YES / NO
- Consistent: YES / PARTIAL / NO
- Issue:

---

## --- ISSUE LIST ---

For each issue:

- ID:
- Location:
- Type:
  - missing_definition
  - inconsistency
  - algorithm_gap
  - code_mismatch
  - math_issue
  - reproducibility_gap
  - parameter_missing
  - structural_violation
  - section_boundary_violation
  - standard_component_overdetail
- Severity:
  - critical
  - major
  - minor
- Description:
- Why it matters:
- Evidence:
- Recommended action:
  - fix_now
  - fix_if_time
  - report_only
- Safe patch:
- Confidence:

---

## --- REPRODUCIBILITY GAPS ---

List explicitly:

- missing step:
- missing parameter:
- missing definition:
- missing condition:

---

## --- STRONGEST REVIEWER ATTACKS ---

1.
2.
3.
4.
5.

---

## --- METHODOLOGY VERDICT ---

Choose one:

- REPRODUCIBLE AND SOUND
- MOSTLY SOUND, MINOR GAPS
- PARTIALLY SPECIFIED
- NOT REPRODUCIBLE
- FUNDAMENTALLY UNCLEAR

Explain:
- what works
- what fails
- what must be fixed

---

## --- MINIMAL SURVIVAL PATCH ---

Provide minimal fixes:

- add:
- clarify:
- align with code:
- remove:
- do NOT change:

Rules:
- do not rewrite full section
- focus on critical gaps
- preserve valid structure

---

## --- SELF-CHECK ---

Confirm:
- I compared method against implementation
- I identified missing steps required for reproduction
- I distinguished clarity issues from real technical gaps
- I did not assume unstated knowledge