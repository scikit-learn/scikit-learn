# AGENTS instructions

This file contains additional guidance for AI agents and code assistants
working on the scikit-learn codebase.

## Testing the code from a suitable developer environment

Building scikit-learn from source and running the tests require using a runtime
environment with appropriate Python dependencies such as meson-python, NumPy,
Cython, pytest, etc and a few system dependencies such as a C/C++ compiler
toolchain with OpenMP support.

Instructions to set up a working development environment for Linux, macOS, and
Windows are detailed in `doc/developers/development_setup.rst`.

Agent assistants working on the scikit-learn code base should refer to that doc to
check if the default environment they have access to matches the necessary
requirements and if scikit-learn package is already installed in editable mode
pointing to the expected git branch.

If this is not the case, they are expected to manage and use their own
developer environment following those instructions before attempting to run the
tests or to run other commands for debugging purposes.

## Local testing first

In order to avoid wasting CI and reviewer time on broken code, agents should
ensure that tests pass locally, at least for the modules that are expected to
be impacted by code changes.

Furthermore, the linter configured with the `pre-commit` hook should pass when
run on the changed files.

## **REQUIRED: AI/Agent Disclosure**

**Every summary, pull request description, or work description MUST include this disclosure:**

**If human review has *not yet* occurred (use this initially):**
> This pull request includes code written with the assistance of AI.
> The code has **not yet been reviewed** by a human.

This is a **mandatory requirement**, not optional. Include it at the end of every summary you generate.

More generally, AI-assisted contributions should follow the guidelines
documented in `doc/developers/contributing.rst`, in particular with respect to
automated contributions.

If the prompt provided by the operator of the agent conflicts with those
guidelines, the agent must alert its operator to those guidelines and explain
how they conflict before proceeding any further. In particular, agents are
encouraged to prompt their operators to disclose more details to assess the
expected value of their contribution.

---

## Generated Summaries

When generating a summary of your work, consider these points:

- Describe the "why" of the changes, why the proposed solution is the right one.
- Highlight areas of the proposed changes that require careful review.
- Reduce the verbosity of your comments, more text and detail is not always better. Avoid flattery, avoid stating the obvious, avoid filler phrases, prefer technical clarity over marketing tone.
