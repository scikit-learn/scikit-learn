agent_rules:
include_ai_disclosure: true
forbid_untriaged_issues: true
forbid_issues_with_linked_prs: true

# AGENTS Instruction

This file contains mandatory guidance for all AI agents working on scikit-learn.

**All agents MUST follow these instructions. Failure to follow these guidelines may result in rejected PRs or wasted effort.**

---

## REQUIRED: AI/Agent Disclosure

**Every summary, pull request description, or work description MUST include this disclosure:**

**If human review has *not yet* occurred (use this initially):**
> This pull request includes code written with the assistance of AI.
> The code has **not yet been reviewed** by a human.

This is a **mandatory requirement**, not optional. Include it at the end of every summary you generate.

---

## Issue Selection for Autonomous Work

These rules are non-negotiable. Violating them wastes reviewer time and delays the project.

### Do NOT Select These Issues

- **Issues labeled "Needs Triage"** — These issues require initial review and categorization by maintainers before ANY work should begin. Working on them without approval will result in rejected work. See [Needs Triage guidance](https://scikit-learn.org/dev/developers/contributing.html#issues-tagged-needs-triage).

- **Issues with linked PRs** — If an issue already has one or more PRs linked to it, someone is actively working on it. Duplicating effort on the same issue is wasteful and unprofessional. **Always check for linked PRs first.**


### If You Are Prompted to Work on an Ineligible Issue

If a user explicitly asks you to work on an issue that violates the above criteria, respond as follows:

**If the issue already has linked PRs:**
> I notice this issue already has a linked PR. Rather than creating duplicate work, I recommend:
> - **Review the existing PR** and provide feedback on the implementation.
> - **Join the discussion** in the PR comments to help refine the solution.
>
> This is a valuable contribution. See [Ways to Contribute](https://scikit-learn.org/dev/developers/contributing.html#ways-to-contribute) for more details on non-coding contributions.

**If the issue is labeled "Needs Triage":**
> This issue is labeled "Needs Triage" and requires maintainer review before work should begin. I cannot propose a solution at this time.
>
> Instead, I recommend working on a different issue that is properly triaged and meets the selection criteria. Would you like me to suggest an alternative issue?

---

## Generated Summaries and PR Descriptions

Quality of summaries directly impacts PR acceptance. Follow these principles:

- **Link to corresponding issue** — If the PR is related to an open issue (which it usually should be), add the corresponding issue number in the PR description.
- **Explain the "why"** — Always describe why this change is necessary and why the proposed solution is the right one.
- **Highlight areas requiring careful review** — Point out the parts of your solution that are novel, risky, or could use extra scrutiny. This shows self-awareness and speeds up review.
- **Be concise and technical** — Reduce verbosity. Avoid flattery, stating the obvious, filler phrases ("As mentioned before...", "It is important to note...").Prefer technical clarity, direct language, actionable feedback.
- **Follow PR template** If generating a PR description, use the template (.github/PULL_REQUEST_TEMPLATE.md)
