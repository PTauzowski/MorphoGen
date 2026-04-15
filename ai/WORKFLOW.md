# AI Writing Workflow

## Abstract pipeline (default)

1. abstract_audit
2. abstract_filter
3. abstract_rewrite
4. jargon_detector
5. abstract_score

Decision:

* Score ≥ 90 → DONE
* 80–89 → repeat from step 2
* < 80 → major rewrite required

---

## Full paper audit

Use:

* agents/paper_wide_audit_agent.md

---

## When to use what

* Poor readability → jargon_detector
* Too long / detailed → abstract_filter
* Structural issues → abstract_audit
* Final decision → abstract_score
