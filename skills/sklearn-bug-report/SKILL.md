---
name: sklearn-bug-report
description: >-
  Help compose a high-quality scikit-learn bug report. Use when the user
  encounters unexpected behavior, an error, or wrong results from scikit-learn
  and wants to file a GitHub issue. Guides triage, builds a minimal reproducer,
  verifies it, and generates text for each template field.
---

# scikit-learn Bug Report Helper

Help the user compose a complete, high-quality bug report for scikit-learn.
Do NOT submit the issue for the user. Produce the text they need to paste into
the GitHub issue form after reviewing it themselves.

## Step 1: Fetch the issue template

Fetch the current bug report template so the output always matches the latest
form fields. If working inside the scikit-learn repo, read the local file
first:

```bash
cat .github/ISSUE_TEMPLATE/bug_report.yml
```

If not in the repo or the file doesn't exist, fetch from GitHub:

```bash
curl -s https://raw.githubusercontent.com/scikit-learn/scikit-learn/main/.github/ISSUE_TEMPLATE/bug_report.yml
```

Parse the YAML to extract each `textarea` entry's `label`, `description`, and
`placeholder`. Use these as the authoritative section names and guidance for the
rest of the workflow. If both sources fail, fall back to these known fields:

1. Describe the bug and give evidence about its user-facing impact
2. Steps/Code to Reproduce
3. Expected Results
4. Actual Results
5. Versions
6. Interest in fixing the bug

## Step 2: Triage

Before drafting anything, determine what kind of problem the user has.
Ask clarifying questions if needed.

**Bug** -- unexpected behavior, wrong results, or an error from scikit-learn's
public API. Proceed to Step 3.

**Usage question** -- the user doesn't know how to use an estimator or
function. Redirect them:
- GitHub Discussions: https://github.com/scikit-learn/scikit-learn/discussions
- Stack Overflow: https://stackoverflow.com/questions/tagged/scikit-learn

**Feature request** -- the user wants behavior that doesn't exist yet.
Redirect to the feature request template:
https://github.com/scikit-learn/scikit-learn/issues/new?template=feature_request.yml

**Already reported** -- always search for existing issues as part of triage:

```bash
gh issue list --repo scikit-learn/scikit-learn --label "Bug" --state open \
  -S "<keywords from the user's description>" --limit 10
```

If a matching issue exists, show it to the user. If the existing issue covers
the same problem, recommend commenting on it rather than filing a duplicate.
If the user wants to proceed with a new issue anyway (e.g., different root
cause or additional context), continue and add "Related: #XXXX" near the end
of the description field.

## Step 3: Gather context

Collect information from the conversation and the user's environment:

1. **Error messages and tracebacks** the user has already shared.
2. **Code snippets** from the conversation or open files.
3. **Version info** -- run this in the user's environment:

```python
import sklearn; sklearn.show_versions()
```

Capture the full output. This is critical -- a bare version number like "1.8.0"
is not sufficient.

4. **Regression check** -- ask the user: "Did this work in a previous version
   of scikit-learn?" If yes, note which version. If the user doesn't know,
   omit this from the description rather than guessing.

## Step 4: Draft the report

Use the hybrid approach: draft everything you can from available context, then
ask the user targeted questions about any gaps.

### Field: Describe the bug and give evidence about its user-facing impact

Write 2-4 sentences covering:

- What the bug is (factual, concise).
- How the user discovered it (the context).
- Why it matters to them as a scikit-learn user (practical impact).
- Whether this is a regression (if known).

Do not use flattery, filler phrases, or marketing tone. Do not write "I believe"
or "it seems like". Be direct.

### Field: Steps/Code to Reproduce

Build a minimal, self-contained Python script. Requirements:

- All imports at the top.
- No external data files. Generate data with `numpy.random`,
  `sklearn.datasets.make_classification`, `sklearn.datasets.make_regression`,
  `sklearn.datasets.make_blobs`, or a few lines of inline Python.
- Only use the public scikit-learn API unless the bug requires private API
  access (justify this in the description if so).
- Set `random_state` where applicable so results are deterministic.
- Remove everything that is not needed to trigger the bug: no unrelated
  preprocessing, no scoring, no plotting, no train/test split unless relevant.
- Follow the project coding conventions (Black).
- The script should be copy-pasteable into a Python terminal and reproduce the
  problem with no modifications.

If the user's original code depends on external data, private datasets, or a
large notebook, help them replace it with synthetic data that triggers the same
issue. Walk through simplification interactively if needed.

### Field: Expected Results

Write a concrete statement of what should happen. Use exact values or behavior,
not vague phrases like "it should work." Examples:

- "No error is raised and `predict` returns an array of shape (100,)."
- "Both unsigned and signed integer targets produce the same
  `decision_function` values."

### Field: Actual Results

Include the complete output. If it's an error, include the full traceback
inside a `python-traceback` fenced code block. If it's wrong values, show the
actual values. Paste output as-is -- don't reformat line wrapping from numpy
or other libraries. Never truncate tracebacks.

If the user or the agent has insights into the root cause, include a
brief analysis.

### Field: Versions

Paste the full output of `sklearn.show_versions()` captured in Step 3.

### Field: Interest in fixing the bug

Ask the user if they would be interested in working on a fix. Do not prefill
this section for the user. This field should be filled out by the human.

## Step 5: Verify the reproducer

**This step is mandatory.** Before presenting the final report, run the
reproducer script from Step 4 in the user's environment. Write the reproducer
to a temporary file and execute it. Do not use `python -c` -- multi-line
scripts with quotes and parentheses break inside single-quoted shell strings.

Check the output:

- If the script **reproduces the bug** (raises the expected error, produces
  wrong values, etc.), proceed to Step 6.
- If the script **does not reproduce the bug** (runs clean, different error,
  different output), iterate: fix the script, adjust the data, or re-examine
  whether the bug exists in the installed version. Do not proceed until the
  reproducer demonstrably fails.
- If the bug cannot be reproduced, tell the user. It may have been fixed in
  their installed version, or the conditions to trigger it may be different
  from what was assumed. Help them investigate further.

## Step 6: Quality checklist

Before presenting the output, verify all of these:

- [ ] Reproducer was executed and confirmed to trigger the bug.
- [ ] Code is self-contained: all imports present, no external data.
- [ ] Data is synthetic (sklearn.datasets, numpy.random, or inline).
- [ ] Follow the project coding conventions (Black).
- [ ] Traceback is complete (not truncated).
- [ ] `sklearn.show_versions()` output is included in full.
- [ ] Impact statement is concrete and explains why the user cares.
- [ ] Existing issues were searched and no duplicate was found.
- [ ] Most related issues and PRs are cross-referenced with a short statements that explains how they are related but not duplicate.

## Step 7: Present the output

Present the completed text to the user with clear section markers matching the
GitHub form field names. Use this format:

---

**Review the text below, edit anything that needs correction, then paste each
section into the corresponding field on the
[scikit-learn bug report form](https://github.com/scikit-learn/scikit-learn/issues/new?template=bug_report.yml).
The form has separate text areas for each section -- do not paste everything
into one field.**

### Describe the bug and give evidence about its user-facing impact

(draft text here)

### Steps/Code to Reproduce

(draft code block here)

### Expected Results

(draft text here)

### Actual Results

(draft text/traceback here)

### Versions

(show_versions output here)

### Interest in fixing the bug

(do not write anything for this section. Remind the user to do that themselves)

---

Remind the user: "Please review and edit the text above before submitting.
AI-generated text can contain inaccuracies. Make sure the description
accurately reflects your experience and the reproducer matches your actual
environment."
