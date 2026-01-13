# AGENTS Instruction

This file contains additional guidance for AI agents and other AI editors.

## **REQUIRED: AI/Agent Disclosure**

**Every summary, pull request description, or work description MUST include this disclosure:**

**If human review has *not yet* occurred (use this initially):**
> This pull request includes code written with the assistance of AI.
> The code has **not yet been reviewed** by a human.

This is a **mandatory requirement**, not optional. Include it at the end of every summary you generate.

---

## Issue Selection for Autonomous Work

When automatically selecting issues to work on, follow these criteria:

**Do NOT select issues with:**
- **"Needs Triage" label** — These issues require initial review and categorization by maintainers before work should begin. See [Needs Triage guidance](https://scikit-learn.org/dev/developers/contributing.html#issues-tagged-needs-triage).
- **Linked pull requests** — If an issue has any PR linked to it, someone is already working on it. Skip to avoid duplicate effort.

**Preferred issue characteristics:**
- Clear problem statement and expected behavior
- No active PRs or ongoing work
- Properly triaged and labeled (not "Needs Triage")
- Reasonable scope for autonomous completion

**When prompted to work on an issue that does NOT meet the criteria:**

- **If the issue already has a linked PR:** Consider reviewing the existing PR and/or joining the discussion instead. This is a valuable contribution. See [Ways to Contribute](https://scikit-learn.org/dev/developers/contributing.html#ways-to-contribute) for more details on non-coding contributions.

- **If the issue is labeled "Needs Triage":** Do not propose a solution. Instead, suggest a different issue to work on that meets the selection criteria above.

---

## Generated Summaries

When generating a summary of your work, consider these points:

- Describe the "why" of the changes, why the proposed solution is the right one.
- Highlight areas of the proposed changes that require careful review.
- Reduce the verbosity of your comments, more text and detail is not always better. Avoid flattery, avoid stating the obvious, avoid filler phrases, prefer technical clarity over marketing tone.
