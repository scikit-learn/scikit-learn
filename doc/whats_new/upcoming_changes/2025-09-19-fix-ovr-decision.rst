Fix OVR Decision and Vectorize SGD/Adam Optimizers
--------------------------------------------------

This PR introduces vectorized implementations of the stochastic optimizers
(SGDOptimizer and AdamOptimizer) used in MLP to improve performance by
leveraging NumPy operations instead of Python loops.

Additionally, it resolves linting issues in `sklearn/preprocessing/_data.py`
by prefixing unused variables with an underscore, ensuring all pre-commit
checks pass.

- Optimized MLP stochastic optimizers for better performance.
- Fixed linting issues in preprocessing module.
- All existing functionality remains intact.
