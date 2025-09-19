.. _fix-ovr-decision:

Fix OVR Decision Function in MLP Optimizers
===========================================

This PR introduces vectorized implementations of the stochastic optimizers
(SGDOptimizer and AdamOptimizer) used in MLP to improve performance by
leveraging NumPy operations instead of Python loops.

Additionally, this PR:

- Resolves linting issues in `sklearn/preprocessing/_data.py` (unused variables prefixed with `_`).
- Ensures all pre-commit checks pass.
- Keeps all existing functionality intact and passes local tests.

Reference Issues/PRs
--------------------

- Fixes #30742 (related optimizer performance issue)
