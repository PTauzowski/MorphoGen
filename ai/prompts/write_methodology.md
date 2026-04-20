---
name: write_methodology
description: Write a precise, implementation-aligned Methodology section grounded in actual code and manuscript
---

ROLE:
Write the Methodology section of the paper with strict alignment to the actual implementation.

INPUT:
- Use the currently focused manuscript
- Use the project codebase (full access via VSCode)
- Use references/ when needed for attribution or comparison

ALWAYS LOAD:
- ai/STYLE_LOCK.md

MISSION:
Produce a Methodology section that is:
- technically precise
- consistent with implementation
- reproducible
- reviewer-verifiable
- structurally clean on first draft

DO NOT:
- invent methods not present in code
- oversimplify novel or modified components
- over-detail standard components
- repeat Introduction or Results
- generate content that belongs in Experiments

CORE PRINCIPLE:
Every methodological claim must be traceable to either:
- code
- equations
- or clearly stated assumptions

---

# PREVENTION RULES (apply before writing, not after)

These rules must be satisfied IN THE FIRST DRAFT.
The structure validator will block the pipeline if they are violated.
Violations here cost a full rewrite pass — avoid them.

## P1 — HEADING DISCIPLINE (STYLE_LOCK §8.7)

- Use `\subsection{}` for major conceptual blocks only
- Use `\subsubsection{}` sparingly — only when the subsection genuinely splits into distinct, substantial sub-topics
- Do NOT use `\paragraph{}` unless the block is ≥ 6 lines of prose AND cannot be merged into running text
- Maximum heading depth: `\subsection{}` → `\subsubsection{}` only
  - Forbidden: `\subsection{}` → `\subsubsection{}` → `\paragraph{}` chains
- If a block has fewer than 5–6 lines under a heading → rewrite as prose, not a headed block
- Consecutive headings with minimal content between them → merge or remove headings

## P2 — SECTION BOUNDARY (STYLE_LOCK §8.8)

Methodology MUST contain:
- method definition
- formulation and equations
- model architecture (novel aspects)
- training mechanisms that are part of the contribution
- variables and assumptions

Methodology MUST NOT contain:
- ablation study design (which factors, how many seeds, which baselines)
- oracle analysis protocols
- evaluation dataset split details
- benchmark setup or baseline comparisons
- standalone §Evaluation Metrics subsections that only define standard metrics
- transfer experiment setup

Evaluation metrics — apply this distinction:
- ALLOWED: metric DEFINITION integral to method formulation
  (e.g., "mIoU excludes the background class because the task targets foreground damage only")
- NOT ALLOWED: evaluation PROTOCOL (split ratios, baseline selection, benchmark setup)
- NOT ALLOWED: a standalone §Evaluation Metrics subsection that only defines standard metrics — place that in Experiments

TEST: If removing a paragraph does not reduce understanding of how the method works → it does not belong in Methodology.

## P3 — STANDARD COMPONENT DISCIPLINE (STYLE_LOCK §8.9)

Distinguish:

### NOVEL or MODIFIED component
- Describe in full
- Include equations where needed
- Be concrete and traceable to code

### STANDARD component (unchanged, published method)
- Name it
- Cite the original paper
- One sentence of functional context is sufficient
- Do NOT reproduce its internal workings, equations, or layer descriptions

EXAMPLES:
- ✓ "The encoder uses a ConvNeXt-B backbone [citation], pretrained on ImageNet."
- ✗ "The UPerNet decoder applies lateral 1×1 convolutions to project each feature map to 256 channels, performs a top-down FPN upsampling pass, concatenates..." (this is over-detail of a standard component)

If a standard component is the object of an ablation → mention it by name in the ablation design (Experiments), not here.

## P4 — STRUCTURE TEMPLATE (ML/SEGMENTATION PAPERS)

Write 3–5 subsections using this framework:

1. **Task and problem formulation**
   - Define what is being predicted (classes, labels)
   - State the evaluation criterion
   - Keep to 1 subsection; use prose, not paragraph heads

2. **Data preparation and preprocessing** (if non-trivial)
   - Describe only non-standard preprocessing steps
   - Standard resize/normalize → one sentence
   - Augmentation → list briefly, no paragraph heads

3. **Model architecture**
   - Describe the novel model design
   - Standard backbone/decoder: name + cite only
   - If multiple architectures are compared: describe the shared design pattern, note variations in a table or inline, do NOT create a subsection per architecture

4. **Training procedure**
   - Loss function: include equation only if non-standard or modified
   - Optimizer: state name and role only — do NOT reproduce update rules
   - Scheduler, stopping criterion: one sentence each
   - If unmodified standard procedure: one short paragraph, no headed blocks

5. **Novel mechanisms** (if applicable)
   - This is where the paper's contribution lives
   - Describe thoroughly: formulation, equations, algorithmic flow
   - Separate clearly from standard components

DO NOT create subsections for:
- Ablation study design → Experiments
- Oracle / threshold analysis → Experiments
- Dataset split rationale → Experiments
- Baseline selection → Experiments
- Evaluation metrics → Experiments or beginning of Results

---

## OPTIMIZER / STANDARD ALGORITHM DISCLOSURE

For ANY standard algorithm (optimizer, scheduler, loss, decoder, backbone):

MUST:
- State the algorithm name
- Cite the source
- State its role in one sentence

MUST NOT:
- Reproduce its update equations
- Describe its internal mechanics
- Present it as a contribution

Exception: if the algorithm is MODIFIED or IS ITSELF a contribution → describe the modification and its motivation. Be explicit about what differs from the standard form.

---

# CORE TASKS

## 1. METHOD IDENTIFICATION

From manuscript + code:

Identify:
- main method / framework
- key components
- computational flow
- inputs / outputs
- which components are novel vs standard

Output internally:
- method structure map
- novel vs standard classification per component

---

## 2. CODE ALIGNMENT

Inspect implementation:

- main functions / entry points
- solver structure
- data flow
- parameters and defaults
- special mechanisms (e.g., conditioning, fusion, filtering)

Check:
- what is actually implemented vs described

Flag internally:
- mismatches
- hidden steps
- implicit assumptions

---

## 3. METHOD DECOMPOSITION

Break method into logical components.

For each component:
- Is it novel or standard?
- If novel → describe in full
- If standard → name + cite + one-sentence role

---

## 4. MATHEMATICAL FORMULATION

Where appropriate:

- define variables
- define objective function
- define key operators

Rules:
- only include math for novel or modified components
- standard loss functions, optimizers → do NOT include equations
- match notation with code logic

---

## 5. STRUCTURE PLANNING (DO BEFORE WRITING)

Before writing, plan the subsection list.

Check against P1–P4:
- Is every planned heading substantive enough to earn a heading?
- Is there any content that belongs in Experiments?
- Are any standard components being over-described?

Only begin writing after the structure plan is validated.

---

# STYLE RULES

- precise and technical
- avoid storytelling
- avoid marketing language ("novel", "robust", "efficient" without evidence)
- avoid "we simply" or "it is straightforward"
- avoid unnecessary adjectives

---

# OUTPUT

## --- METHODOLOGY ---
<final section>

---

## --- STRUCTURE MAP ---

List all headings used:
- subsection:
- subsubsection (if any):
- paragraph (if any — must be justified):

Flag any heading that covers fewer than 5 lines of prose.

---

## --- NOVEL VS STANDARD CLASSIFICATION ---

For each component:
- Component:
- Novel / Standard / Modified:
- Detail level applied:

---

## --- SECTION BOUNDARY CHECK ---

List all subsections.
For each, confirm: "This describes how the method works, not how the experiments were run."

Flag any subsection that fails this test.

---

## --- CODE TRACEABILITY MAP ---

For each major novel component:
- Component:
- Code location:
- Role in method:

---

## --- POTENTIAL MISMATCHES ---

- paper says:
- code does:
- risk level:

---

## --- SELF-CHECK ---

Confirm:
- [ ] No `\paragraph{}` used for blocks under 6 lines
- [ ] No subsection → subsubsection → paragraph chains
- [ ] No ablation design, oracle protocol, or dataset split in Methodology
- [ ] Standard components described with name + cite + one sentence only
- [ ] Novel components described with full formulation
- [ ] All major steps are grounded in implementation
- [ ] No invented algorithmic elements
