Adaptive repeated Lasso CV: benchmark proposal
============================================

This prototype is an opt-in benchmark using the public ``lasso_path`` and
``Lasso`` APIs. It does not add an estimator to scikit-learn's public API.
The purpose is to evaluate the statistical and computational tradeoffs before
proposing API changes. Dense, finite, single-output, unweighted regression
with an intercept is supported. The alpha grid is supplied in advance and
held fixed. Feature scaling is outside this prototype; when needed, it must
be fitted on each training fold.

Each repetition draws one independently shuffled K-fold partition. A warm
Lasso path is fitted within each fold after centering on training data only.
Fold MSEs are averaged equally, matching ``LassoCV``. Complete partition loss
curves are averaged before choosing the penalty. The procedure never averages
the selected penalties and never refits earlier partitions. It fits the final
selected penalty once on the full data.

The stopping statistic is each loss relative to the largest alpha. Its MCSE
must be below tolerance / (sqrt(2) * z) for every nonreference alpha, where
z uses error_rate / (m - 1) for a grid of m penalties. This Bonferroni allocation
is conservative and only supplies a normal approximation to simultaneous
reproducibility. It is not a finite-sample sequential confidence guarantee.
The minimum repetition count is burn-in, not proof of calibration. A budget
limit returns ``budget_exhausted``.

Risk precision and model selection are separate. The report includes the
single-partition winner histogram and the alphas within twice the requested
tolerance of the best observed mean loss. These are descriptive diagnostics;
the latter is not a confidence set. Independent algorithm runs are needed
to measure agreement of selected penalties and active sets.

Comparison methods
------------------

The benchmark compares single shuffled CV, fixed 10/50/100/500 repetitions,
adaptive repetition, deterministic unshuffled CV, and approximate leave-one-out
(ALO). ALO uses active-set leverage including the unpenalized intercept.
Its full-data path rescales alpha by (n - 1) / n so deletion corresponds to
holding scikit-learn's mean-loss penalty constant. The approximation is exact
when deleting an observation preserves the active set and coefficient signs,
and can fail when these change. Leverage near one produces a clear error.
ALO and K-fold CV target different validation criteria.

The unit tests compare full loss curves against ``LassoCV`` on identical folds
and against ``RepeatedKFold`` for a fixed repetition budget. ALO is checked
against exact leave-one-out fits in a stable-active-set example and the
intercept-only limit. Tests also count path fits and final refits, verify the
relative-risk MCSE and multiplicity adjustment, and cover early stopping,
budget exhaustion, invalid arguments, ties, and constant relative risks.

Run from the repository root with a source build installed::

    python benchmarks/bench_lasso_reproducibility.py --pairs 1000 --n-jobs 4
    pytest benchmarks/tests/test_lasso_reproducibility.py

Each pair comprises independent executions on the same fixed training data.
The report includes joint risk-curve exceedance with a Wilson 95% Monte Carlo
interval, alpha agreement, exact active-set agreement, held-out MSE, actual
path fit counts, budget exhaustion, and elapsed time. Deterministic methods
have perfect repeatability by construction; this does not establish predictive
superiority. The held-out data are never used for tuning or stopping.

This single-data experiment cannot establish that adaptive CV never worsens
predictions. Broader comparisons across sample sizes, correlations, sparsity,
and signal strengths are needed before an upstream estimator proposal. Timing
depends on hardware and concurrent work; path fit counts are more portable.

References
----------

* Ritzwoller and Romano (2026), https://doi.org/10.1257/aer.20231651.
* Wang et al. (2018), https://proceedings.mlr.press/v80/wang18m.html.

Recorded comparison
-------------------

With fixed data seed 1234, 1,000 independent pairs, tolerance 0.02, nominal
error 0.05, burn-in 20, and maximum 1,000 repetitions:

.. list-table:: Conditional reproducibility and predictive performance
   :header-rows: 1

   * - Method
     - Curve exceedance
     - Alpha agreement
     - Mean path fits/run
     - Test MSE
   * - Single shuffled
     - 99.0%
     - 34.5%
     - 5
     - 1.103509
   * - Fixed 10
     - 72.1%
     - 93.5%
     - 50
     - 1.101577
   * - Fixed 50
     - 21.5%
     - 100%
     - 250
     - 1.101473
   * - Fixed 100
     - 7.8%
     - 100%
     - 500
     - 1.101473
   * - Fixed 500
     - 0.0%
     - 100%
     - 2500
     - 1.101473
   * - Adaptive
     - 0.3%
     - 100%
     - 1468.123
     - 1.101473
   * - Unshuffled
     - 0.0%
     - 100%
     - 5
     - 1.111579
   * - ALO
     - 0.0%
     - 100%
     - 1
     - 1.111579

Adaptive CV used 293.6245 repetitions on average, without budget exhaustion.
Its exceedance interval is 0.102%--0.878%, well below the nominal 5%, consistent
with a conservative Bonferroni rule on highly correlated relative risks.
The observed zero exceedances for fixed 500 have an upper 95% Monte Carlo
limit of 0.383%; they do not establish zero probability.

Fifty fixed repetitions already produced identical selected penalties and
predictions across these runs, although their risk curves did not meet the
requested precision. Adaptive repetition therefore bought risk-curve precision
without improving prediction over fixed 50 in this example. ALO was much
cheaper and had slightly worse test MSE on this dataset. Neither observation
justifies a general predictive ranking or a new public estimator API.

``results/lasso-reproducibility.json`` contains the raw results, including
active-set agreement, Monte Carlo intervals, and timings. Four worker processes
were used with Python 3.12.12, NumPy 2.5.4, SciPy 1.18.1, and a source build of
scikit-learn 1.10.dev0 based on upstream revision
``84fc20f30f5150b415c11ac9a20123f6b2de876d``.

This pull request includes code written with the assistance of AI.
`Human review status <https://github.com/scikit-learn/scikit-learn/pull/35288>`_ is tracked in the pull request.
