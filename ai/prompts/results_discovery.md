---
name: results_discovery
description: Discover likely numerical result sources and propose a draft results manifest for the current paper
---

ROLE:
Inspect the local project repository and discover likely numerical result sources for the current paper.

INPUT:
- Use the currently focused manuscript for alignment
- Inspect the local repository filesystem
- Use code, file names, folders, and manuscript references to infer likely result locations
- Do NOT require an existing results manifest

ALWAYS LOAD:
- ai/STYLE_LOCK.md
- ai/config/results_manifest_schema.md
- ai/ARTIFACT_STATE_CONVENTION.md

MISSION:
Produce a structured draft of the numerical results manifest by identifying:
- likely authoritative result files
- likely stale/debug/archive outputs
- likely experiment groups
- likely primary metrics
- likely reporting priorities

CORE PRINCIPLE:
This stage proposes a draft manifest.
It must be conservative, evidence-based, and explicit about uncertainty.

DO NOT:
- write the numerical results subsection
- invent experiments not evidenced by files or manuscript
- treat guessed files as authoritative without marking uncertainty
- silently exclude plausible result sources
- validate the manifest as final

OUTPUT TYPE:
A suggested manifest draft plus a discovery report.

---

# DISCOVERY TASKS

## 1. MANUSCRIPT-GUIDED RESULT CUE EXTRACTION

Read the focused manuscript and identify references to:
- tables
- figures
- result folders
- experiment names
- metric names
- baseline names
- ablations
- sensitivity studies
- validation runs

Extract likely cues such as:
- “Table 2 reports…”
- “ablation”
- “final comparison”
- “sensitivity”
- “validation”
- “benchmark”
- “envelope”
- “averaged”
- “mIoU”
- “volume reduction”
- “buckling factor”

Purpose:
Infer what kinds of numerical result groups should exist.

---

## 2. FILESYSTEM SCAN

Inspect the repository for likely result-containing paths.

Prioritize directories and files such as:
- results/
- output/
- outputs/
- experiments/
- runs/
- logs/
- paper/results/
- figures/
- tables/
- csv/
- json/

Look for file types such as:
- .csv
- .json
- .md
- .txt
- .tex tables
- structured log summaries

Detect possible anti-patterns:
- debug/
- tmp/
- archive/
- old/
- backup/
- draft/
- duplicate v1/v2/v3 chains

Purpose:
Generate candidate authoritative and excluded sources.

---

## 3. RESULT FILE CANDIDATE RANKING

Rank candidate result files by likelihood that they are paper-relevant.

Signals that increase confidence:
- mentioned or implied in manuscript
- located in “final” or “paper/results” folders
- tabular/summary format
- stable naming
- metrics matching paper claims

Signals that lower confidence:
- debug/tmp/archive folders
- per-iteration raw dumps
- duplicated versions without clear final marker
- exploratory-only naming

Classify each candidate:
- HIGH confidence
- MEDIUM confidence
- LOW confidence

---

## 4. EXPERIMENT GROUP INFERENCE

Infer likely experiment groups from:
- manuscript structure
- filenames
- folder names
- metric groupings

Typical groups may include:
- baseline comparison
- main method comparison
- ablation study
- sensitivity study
- final validation
- robustness check
- efficiency/runtime summary

For each inferred group:
- name
- purpose
- candidate source files
- confidence

---

## 5. PRIMARY METRIC INFERENCE

Infer likely primary metrics from:
- manuscript text
- table headers
- CSV headers
- figure captions
- result filenames

For each metric:
- proposed id
- proposed name
- likely meaning
- likely unit
- confidence

Do NOT invent units if absent.
Use UNKNOWN if needed.

---

## 6. EXCLUSION INFERENCE

Identify likely paths/files that should be excluded from the final manifest.

Typical exclusion reasons:
- debug output
- temporary file
- old version
- exploratory run
- archived experiment
- failed run
- raw solver trace without summary role

Mark exclusions conservatively.
If a file looks important but uncertain, do NOT exclude aggressively.

---

## 7. CLAIM LINK INFERENCE

Infer likely mapping between manuscript claims and experiment groups.

Examples:
- “conservative strategy is denser” → strategy comparison group
- “main accuracy improvement” → benchmark table
- “stability maintained” → final validation or safety metric group

Mark confidence:
- HIGH / MEDIUM / LOW

---

## 8. DRAFT MANIFEST GENERATION

Using all prior steps, generate a draft results manifest compatible with:
- ai/config/results_manifest_schema.md

This is a suggested manifest only.
It must separate:
- authoritative sources
- excluded sources
- experiment groups
- metrics
- comparison rules
- writing constraints

Any uncertain field must be explicitly marked.

---

# OUTPUT

## --- DISCOVERY SUMMARY ---

- Manuscript result cues:
- Candidate result directories:
- Candidate result files:
- Likely stale/debug paths:
- Discovery confidence: HIGH / MEDIUM / LOW

---

## --- CANDIDATE SOURCE RANKING ---

For each candidate:

- Path:
- Type:
- Likely role:
- Confidence:
- Reason:

---

## --- INFERRED EXPERIMENT GROUPS ---

For each group:

- Group name:
- Purpose:
- Candidate files:
- Confidence:

---

## --- INFERRED PRIMARY METRICS ---

For each metric:

- Metric ID:
- Name:
- Meaning:
- Unit:
- Confidence:
- Source evidence:

---

## --- EXCLUSION SUGGESTIONS ---

- Path:
- Reason:
- Confidence:

---

## --- CLAIM LINK SUGGESTIONS ---

- Manuscript claim:
- Supporting experiment group:
- Supporting metric(s):
- Confidence:

---

## --- DRAFT RESULTS MANIFEST ---
Provide a draft manifest in the structure required by:
- ai/config/results_manifest_schema.md

Mark uncertain fields explicitly, for example:
- "UNCERTAIN"
- "REVIEW_NEEDED"

---

## --- REQUIRED HUMAN CHECKS ---

List the places where human confirmation is most needed.

Examples:
- selecting final run among multiple similar files
- confirming stale vs final version
- confirming metric units
- confirming which experiment group belongs in the paper

---

## --- DOWNSTREAM DECISION ---

Choose one:

- READY TO CREATE MANIFEST
- MANIFEST NEEDS HUMAN REVIEW
- INSUFFICIENT REPOSITORY SIGNAL

Explain why.

---

## --- SELF-CHECK ---

Confirm:
- Discovery was based on manuscript + filesystem signals
- No file was treated as authoritative without evidence
- Uncertainty was preserved rather than hidden
- This stage produced a draft manifest, not final reporting