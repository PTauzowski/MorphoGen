---
name: abstract_validate
description: Final validation of abstract
---

TASK:
Validate the abstract.

CHECK:

1. Structure:
- Background present?
- Method briefly stated?
- Results quantitative?
- Conclusion clear?
- 2–3 sentences per section (Background, Methods, Results, Conclusion)?

2. Density:
- too many numbers?
- redundant numbers?

3. Claims:
- any overstatement?
- any ambiguity?

4. Readability:
- smooth narrative?
- no checklist tone?

5. Consistency:
- matches results section?
- matches contribution?

FAIL CONDITIONS (override score — mark NOT READY if any apply):
- any structural section missing (Background, Methods, Results, Conclusion)
- any structural section has fewer than 2 sentences
- any structural section has more than 3 sentences
- contains a citation or reference (e.g., \cite{}, [1], (Author, Year))
- references a figure or table (e.g., Fig., Table, Figure)
- uses an abbreviation not defined within the abstract itself

OUTPUT:

--- VALIDATION REPORT ---
- STATUS:
  READY (≥90, no fail conditions)
  BORDERLINE (80–89, no fail conditions)
  NOT READY (<80 or any fail condition triggered)

--- ISSUES (if any) ---

--- FINAL VERDICT ---
Ready for submission / Needs revision