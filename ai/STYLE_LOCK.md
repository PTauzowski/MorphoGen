# SCIENTIFIC WRITING STYLE POLICY

This document defines non-negotiable writing constraints for all stages of the paper pipeline.

It applies to:
- all sections
- all prompts
- all rewrite/patch stages

---

# 1. GLOBAL PRIORITIES

1. Preserve scientific correctness
2. Preserve readability and narrative flow
3. Apply minimal necessary changes
4. Prefer reporting issues over rewriting when uncertain

---

# 2. GLOBAL HARD CONSTRAINTS

- Do NOT increase verbosity unless explicitly required
- Do NOT convert prose into checklist-style writing
- Do NOT rewrite entire sections when local edits suffice
- Do NOT introduce new claims without evidence
- Do NOT invent, approximate, or modify numerical results
- Inline-only outputs are INVALID — outputs must exist in ai/out/

---

# 3. GLOBAL CONSISTENCY RULES

## 3.1 Terminology Consistency
- Use one term per concept throughout the paper
- Do not alternate synonyms unless distinction is intentional
- Terms must match across Abstract, Introduction, Results, Conclusion

## 3.2 Abbreviation Rule
- Define every non-standard abbreviation on first use
- Do NOT redefine abbreviations later
- Avoid abbreviations used fewer than 3 times
- Prefer full term if readability suffers

## 3.3 Symbol & Notation Rule
- Define all symbols before use
- Use consistent notation across the manuscript
- Do not reuse symbols for different meanings unless explicitly reset

---

# 4. CLAIM–EVIDENCE DISCIPLINE

- Every claim must be supported by presented evidence
- Strong claims require strong validation
- If evidence is partial → claim must be softened
- Do NOT generalize beyond tested conditions

---

# 5. STANDARD METHOD SUPPRESSION

- Do NOT present standard methods as contributions
- Standard tools (e.g., FEM, CNNs, optimizers) must be:
  - briefly described
  - positioned as implementation components
- Novelty must be clearly separated from standard components

---

# 6. FIGURE & TABLE INTEGRATION

- Every figure/table must be interpreted in text
- Do NOT dump figures without explanation
- Do NOT repeat full table contents in prose
- Each figure must support a specific claim

---

# 7. SECTION ANTI-REDUNDANCY RULE

- Do NOT restate the same contribution across sections
- Each section has a distinct role:
  - Introduction → motivation & gap
  - Method → how it works
  - Results → what was observed
  - Conclusion → what is proven

---

# 8. SECTION-SPECIFIC RULES

---

## 8.1 ABSTRACT RULES

### Purpose:
Summarize the research question, method, and key findings.

### Structure:
- Background → Methods → Results → Conclusion
- 2–3 sentences per section with logical flow between them

### Constraints:
- Maximum 4 numeric values
- Report results at summary level only
- Do NOT include:
  - implementation details
  - training details
  - architecture internals
  - standard deviations or redundant metrics
  - citations or references (e.g., \cite{}, [1])
  - figure or table references (e.g., Fig., Table)
  - abbreviations not defined within the abstract itself


### Language:
- Focus on *what was tested and found*, not how
- Replace ML jargon with functional meaning

---

## 8.2 INTRODUCTION RULES

### Purpose:
Build the research argument.

### Requirements:
- Clearly state:
  - problem
  - gap
  - contribution
- Prefer quantitative narration where possible
- Avoid generic statements

### Forbidden:
- Literature dump without argument
- Overpromising results

---

## 8.3 RELATED WORK RULES

### Purpose:
Position the work in the literature.

### Requirements:
- Group by method/function, not author sequence
- Compare approaches, not papers
- Identify limitations that motivate your work

### Forbidden:
- Citation lists without synthesis
- Unsupported criticism

---

## 8.4 METHODOLOGY RULES

### Purpose:
Explain how the method works.

### Requirements:
- Describe actual implemented method
- Define all variables and assumptions
- Use equations only where needed

### Forbidden:
- Textbook explanations without direct relevance
- Over-detailed implementation of standard components → see §8.9
- Heading spam, `\paragraph{}` overuse, or deep heading chains → see §8.7
- Experimental setup, ablation design, dataset splits, or evaluation protocols → see §8.8

---

## 8.5 RESULTS RULES

### Purpose:
Present and interpret findings.

### Requirements:
- Report results relevant to claims
- Interpret results (not just present them)
- Compare against baselines or alternatives

### Forbidden:
- Dumping logs, raw outputs, or full experiment traces
- Mixing debug results with final results
- Claiming effects not directly supported

---

## 8.6 CONCLUSION RULES

### Purpose:
Close the argument.

### Requirements:
- Summarize validated contributions
- Reflect actual findings only
- Mention limitations proportionally

### Forbidden:
- Introducing new results
- Expanding claims beyond evidence

---

# 8.7 STRUCTURAL DISCIPLINE RULES

## Heading Usage

- Use headings only when they introduce a meaningful conceptual unit
- Do NOT use headings for blocks shorter than ~5–6 lines of prose
- Prefer continuous prose over multiple short headed blocks
- Avoid consecutive headings with minimal content between them

## Hierarchy Depth

- Maximum allowed depth: subsection → subsubsection
- Use of \paragraph{} is strongly discouraged
- Do NOT create subsection → subsubsection → paragraph chains

## Fragmentation

- Avoid splitting logically continuous explanations into multiple headed blocks
- If two consecutive blocks could be merged without loss of clarity → they MUST be merged

---

# 8.8 SECTION IDENTITY RULES

Each section must preserve its role in the global argument.

## Methodology vs Experiments

Methodology MUST:
- define the method
- describe formulation and algorithm

Methodology MUST NOT:
- describe evaluation protocols
- describe ablation study design
- define dataset splits
- describe benchmark setup

If removing a paragraph does not reduce understanding of the method,
that paragraph does NOT belong in Methodology.

---

## Results vs Discussion

Results MUST:
- report observed outcomes
- interpret them locally

Results MUST NOT:
- generalize beyond evidence
- provide literature-level interpretation

---

## Conclusion Discipline

Conclusion MUST:
- reflect only validated findings

Conclusion MUST NOT:
- introduce new results
- strengthen claims beyond Results

---

# 8.9 STANDARD COMPONENT DETAIL RULE (REFINED)

- Standard components must be:
  - named
  - referenced
  - briefly described only if needed for understanding

- Detailed explanation is allowed ONLY if:
  - the component is modified, OR
  - it directly interacts with the novel contribution

- If removing the description does not harm understanding of the contribution,
  the description MUST be shortened or removed.

---

# 8.10 TENSE RULES

Use the correct grammatical tense for each section. Tense signals the epistemic status of a claim.

| Section | Rule |
|---|---|
| Abstract — background/methods | Past tense |
| Abstract — implications/conclusions | Present tense |
| Introduction — established knowledge | Present tense |
| Introduction — own prior work, cited in gap | Past tense |
| Methodology | Past tense throughout |
| Results — findings | Past tense |
| Results — conclusions drawn from findings | Present tense |
| Discussion — interpretation of own results | Present tense |
| Discussion — own results being referenced | Past tense |
| Conclusion — validated contributions | Present tense |

Violations degrade scientific register and signal non-native or unchecked writing.

---

# 9. NUMERIC DISCIPLINE (GLOBAL)

- Numbers must be:
  - correct
  - non-redundant
  - essential

Do NOT:
- duplicate the same information in multiple forms
- report derived values if already implied
- overload sentences with metrics

---

# 10. JARGON CONTROL RULE

- Prefer functional explanation over internal ML terminology
- A non-AI domain expert must understand each sentence
- If >1 specialized term appears → simplify

---

# 11. SEVERITY DEFINITIONS

CRITICAL:
- changes scientific meaning
- introduces unsupported claims
- violates abstraction level
→ MUST be fixed

MAJOR:
- reduces clarity or accessibility
→ SHOULD be fixed

MINOR:
- stylistic improvement
→ OPTIONAL

---

# 12. SELF-CHECK BEFORE FINALIZING

- Did readability improve?
- Did scientific meaning remain unchanged?
- Did local edits remain local?
- Did verbosity increase unnecessarily?

If any answer is problematic → revise

---

# 13. PROHIBITED ACTIONS

DO NOT:
- invent numerical results
- average values unless explicitly defined
- approximate reported values
- silently modify claims

---

# 14. PIPELINE ALIGNMENT RULE

A section is valid ONLY if:
- it follows STYLE_LOCK.md
- AND it exists as a file in ai/out/

Inline-only compliance is INVALID.