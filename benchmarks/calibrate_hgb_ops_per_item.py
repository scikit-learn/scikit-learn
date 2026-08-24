"""Calibrate the "simple ops per item" constants used throughout the
HistGradientBoosting Cython code to decide whether a loop is worth
parallelizing (see `_use_threads_for_workload` in
sklearn/utils/_openmp_helpers.pyx, and its callers in
sklearn/ensemble/_hist_gradient_boosting/{_predictor,_gradient_boosting,
histogram,splitting}.pyx).

For each operation and a few representative workload shapes (varying
samples/features ratio, with/without missing values, with/without
categorical features, constant/variable hessians), this:

1. measures real, single-threaded (n_threads=1, so `_use_threads_for_workload`
   never fires) wall-clock time across several workload sizes;
2. fits a line to isolate the marginal per-item cost from fixed overhead;
3. converts that per-item cost to the same "instructions per second" unit
   used by `_use_threads_for_workload` (see
   `sklearn.utils._openmp_helpers._use_threads_for_workload`), so the
   numbers below are directly comparable to, and can replace, the
   hardcoded `ops_per_item` constants in the four files above.

It also runs a small sanity check: for one operation, it predicts (from the
measured ops_per_item and the machine's calibrated per-thread OpenMP
dispatch overhead, see `_omp_prange_dispatch_overhead`) the workload size at
which parallelizing starts to pay off for a few thread counts, then
empirically checks it by timing the real code path (which makes its own
threading decision) around that size.

Run directly: `python calibrate.py`
"""

import time

import numpy as np
from sklearn.ensemble._hist_gradient_boosting.common import (
    G_H_DTYPE,
    X_BINNED_DTYPE,
    Y_DTYPE,
    MonotonicConstraint,
)
from sklearn.ensemble._hist_gradient_boosting.histogram import HistogramBuilder
from sklearn.ensemble._hist_gradient_boosting.splitting import (
    Splitter,
    compute_node_value,
)
from sklearn.utils._openmp_helpers import _omp_prange_dispatch_overhead

from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower

rng = np.random.default_rng(0)

# Same convention as `_calibrate_min_instructions_per_thread`: the historical
# default of 2000 instructions was picked as the equivalent of ~1us of simple
# work.
INSTRUCTIONS_PER_SECOND = 2e9

# Current hardcoded constants, for side-by-side comparison in the report
# below. Keep this in sync by hand with the .pyx source if it changes.
CURRENT_VALUES = {
    "histogram populate (per sample)": {"constant hessian": 1, "variable hessian": 2},
    "histogram build (per feature x sample)": {
        "constant hessian": 1,
        "variable hessian": 2,
    },
    "histogram subtraction (per feature x bin)": {"": 3},
    # splitting.pyx uses a single flat constant regardless of feature kind
    # (deliberately conservative for the common, mostly-plain-numerical
    # case), so the same current value applies to all three shapes here.
    "find_node_split (per feature x bin)": {
        "plain": 5,
        "with missing values": 5,
        "categorical": 5,
    },
    "split_indices (per sample)": {
        "plain": 3,
        "with missing values": 7,
        "categorical": 9,
    },
    "predict (per sample)": {"binned": 25, "raw": 25},
    "update_raw_predictions (per sample)": {"": 1},
}


def best_of_seconds(fn, n_repeats=7):
    """Minimum wall-clock time over `n_repeats` calls to `fn()`.

    Uses the minimum, not the mean/median: external interference (CPU
    frequency scaling, other processes, scheduler preemption, ...) can only
    ever slow a measurement down, never speed it up below the true cost, so
    the minimum across samples is the one least contaminated by such noise
    (the same reasoning `timeit`, and this repo's own
    `_calibrate_min_instructions_per_thread`, use).
    """
    best = None
    for _ in range(n_repeats):
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        if best is None or dt < best:
            best = dt
    return best


def fit_marginal_ops_per_item(predictors, seconds):
    """Least-squares fit of `seconds ~ intercept + predictors @ coef`.

    `predictors` has shape (n_points, n_terms) and holds, for each
    measurement, the work-item counts (e.g. n_samples, n_features*n_samples)
    that `_use_threads_for_workload` would be given. Returns `coef`,
    converted to instructions/item, dropping the fitted intercept (which
    absorbs fixed, size-independent overhead such as allocations).
    """
    design = np.column_stack([np.ones(len(seconds)), predictors])
    coef, *_ = np.linalg.lstsq(design, seconds, rcond=None)
    return coef[1:] * INSTRUCTIONS_PER_SECOND


def report(title, values_by_shape):
    print(f"\n{title}")
    current = CURRENT_VALUES.get(title, {})
    for shape, ops_per_item in values_by_shape.items():
        # A near-zero true cost can fit slightly negative under measurement
        # noise; clip it, since a negative "instructions per item" isn't
        # physically meaningful.
        ops_per_item = max(0.0, ops_per_item)
        current_value = current.get(shape)
        suffix = f"  (current: {current_value})" if current_value is not None else ""
        label = f"  {shape}:" if shape else ":"
        print(f"{label:<28} {ops_per_item:8.2f} instructions/item{suffix}")


# -----------------------------------------------------------------------
# Histogram building (populate + per-feature-per-sample build), from
# HistogramBuilder.compute_histograms_brute.
# -----------------------------------------------------------------------
def bench_histogram_brute():
    # `compute_histograms_brute`'s time is, to good approximation,
    # `intercept + populate_ops * n_samples + build_ops * n_features * n_samples`
    # (plus a small, n_samples-independent allocation cost). Jointly fitting
    # all 3 terms from a 2D (n_samples, n_features) sweep is ill-conditioned
    # here: times span 3 orders of magnitude, so the fit is dominated by the
    # largest ones and can't reliably separate the two terms that both scale
    # with n_samples. Instead, fix n_samples large (so per-call overhead is
    # negligible) and sweep only n_features: at fixed n_samples, time becomes
    # a simple 1D line in n_features, whose intercept and slope map directly
    # to `populate_ops` and `build_ops`.
    n_bins = 256
    n_samples_total = 400_000
    n_samples = n_samples_total // 2  # non-root node; see comment below
    n_features_grid = np.array([2, 10, 30, 60, 100, 200])

    results = {}
    for hessians_are_constant, label in [
        (True, "constant hessian"),
        (False, "variable hessian"),
    ]:
        seconds = []
        for n_features in n_features_grid:
            X_binned = np.asfortranarray(
                rng.integers(0, n_bins - 1, size=(n_samples_total, n_features)),
                dtype=X_BINNED_DTYPE,
            )
            gradients = rng.random(n_samples_total).astype(G_H_DTYPE)
            hessians = (
                np.ones(1, dtype=G_H_DTYPE)
                if hessians_are_constant
                else rng.random(n_samples_total).astype(G_H_DTYPE)
            )
            hb = HistogramBuilder(
                X_binned, n_bins, gradients, hessians, hessians_are_constant, 1
            )
            # Use half the samples as a non-root node, to exercise the
            # "populate" reordering step (skipped for the root node, where
            # sample_indices already covers the whole array).
            sample_indices = np.arange(n_samples, dtype=np.uint32)
            dt = best_of_seconds(lambda: hb.compute_histograms_brute(sample_indices))
            seconds.append(dt)

        design = np.column_stack([np.ones_like(n_features_grid), n_features_grid])
        (intercept, slope), *_ = np.linalg.lstsq(design, np.array(seconds), rcond=None)
        populate_ops = intercept * INSTRUCTIONS_PER_SECOND / n_samples
        build_ops = slope * INSTRUCTIONS_PER_SECOND / n_samples
        results.setdefault("histogram populate (per sample)", {})[label] = populate_ops
        results.setdefault("histogram build (per feature x sample)", {})[label] = (
            build_ops
        )

    for title, values in results.items():
        report(title, values)


# -----------------------------------------------------------------------
# Histogram subtraction, from HistogramBuilder.compute_histograms_subtraction.
# -----------------------------------------------------------------------
def bench_histogram_subtraction():
    predictors = []
    seconds = []
    for n_bins in (64, 256):
        for n_features in (4, 20, 100, 400):
            hist_dtype = HistogramBuilder(
                np.asfortranarray(np.zeros((2, n_features), dtype=X_BINNED_DTYPE)),
                n_bins,
                np.zeros(2, dtype=G_H_DTYPE),
                np.ones(1, dtype=G_H_DTYPE),
                True,
                1,
            ).compute_histograms_brute(np.arange(2, dtype=np.uint32))
            parent = np.asarray(hist_dtype).copy()
            sibling = np.asarray(hist_dtype).copy()

            def run():
                HistogramBuilder(
                    np.asfortranarray(np.zeros((2, n_features), dtype=X_BINNED_DTYPE)),
                    n_bins,
                    np.zeros(2, dtype=G_H_DTYPE),
                    np.ones(1, dtype=G_H_DTYPE),
                    True,
                    1,
                ).compute_histograms_subtraction(parent.copy(), sibling)

            dt = best_of_seconds(run)
            predictors.append([n_features * n_bins])
            seconds.append(dt)

    (sub_ops,) = fit_marginal_ops_per_item(np.array(predictors), np.array(seconds))
    report("histogram subtraction (per feature x bin)", {"": sub_ops})


# -----------------------------------------------------------------------
# Splitter set up shared by find_node_split and split_indices benchmarks.
# -----------------------------------------------------------------------
N_CATEGORIES = 20  # number of distinct categories used for shape="categorical"


def _bins_scanned_per_feature(shape, n_bins):
    """Bins actually scanned per feature by find_node_split/split_indices.

    Categorical splits only scan the `n_bins_non_missing` bins actually in
    use (N_CATEGORIES here), not the nominal `n_bins` — using the latter as
    the work-item count would understate the true per-item cost.
    """
    return N_CATEGORIES if shape == "categorical" else n_bins


def _make_splitter_and_histograms(n_samples, n_features, n_bins=256, shape="plain"):
    missing_values_bin_idx = n_bins - 1

    if shape == "categorical":
        col = rng.integers(0, N_CATEGORIES, size=n_samples)
    else:
        col = rng.integers(0, n_bins - 1, size=n_samples)

    if shape == "with missing values":
        missing_mask = rng.random(n_samples) < 0.1
        col = col.copy()
        col[missing_mask] = missing_values_bin_idx

    X_binned = np.asfortranarray(
        np.tile(col.reshape(-1, 1), (1, n_features)), dtype=X_BINNED_DTYPE
    )
    n_bins_non_missing = np.full(
        n_features,
        N_CATEGORIES if shape == "categorical" else n_bins - 1,
        dtype=np.uint32,
    )
    has_missing_values = np.full(
        n_features, shape == "with missing values", dtype=np.uint8
    )
    is_categorical = np.full(n_features, shape == "categorical", dtype=np.uint8)
    monotonic_cst = np.full(n_features, MonotonicConstraint.NO_CST, dtype=np.int8)

    gradients = rng.standard_normal(n_samples).astype(G_H_DTYPE)
    hessians = rng.random(n_samples).astype(G_H_DTYPE)

    hb = HistogramBuilder(X_binned, n_bins, gradients, hessians, False, 1)
    histograms = hb.compute_histograms_brute(np.arange(n_samples, dtype=np.uint32))

    splitter = Splitter(
        X_binned,
        n_bins_non_missing,
        missing_values_bin_idx,
        has_missing_values,
        is_categorical,
        monotonic_cst,
        l2_regularization=0.0,
        hessians_are_constant=False,
        n_threads=1,
    )
    sum_gradients = float(gradients.sum())
    sum_hessians = float(hessians.sum())
    value = compute_node_value(sum_gradients, sum_hessians, -np.inf, np.inf, 0.0)
    return splitter, histograms, sum_gradients, sum_hessians, value


# -----------------------------------------------------------------------
# find_node_split, from Splitter.find_node_split.
# -----------------------------------------------------------------------
def bench_find_node_split():
    n_samples = 5_000
    values = {}
    for shape in ("plain", "with missing values", "categorical"):
        # For "categorical", only n_features varies the true workload
        # (n_bins doesn't: only N_CATEGORIES bins are ever scanned per
        # feature), so sweep n_features more finely there instead of also
        # sweeping n_bins.
        n_bins_grid = (256,) if shape == "categorical" else (64, 256)
        n_features_grid = (
            (4, 8, 20, 50, 100) if shape == "categorical" else (4, 20, 100)
        )

        predictors = []
        seconds = []
        for n_bins in n_bins_grid:
            bins_per_feature = _bins_scanned_per_feature(shape, n_bins)
            for n_features in n_features_grid:
                splitter, histograms, sum_g, sum_h, value = (
                    _make_splitter_and_histograms(
                        n_samples, n_features, n_bins=n_bins, shape=shape
                    )
                )
                dt = best_of_seconds(
                    lambda: splitter.find_node_split(
                        n_samples, histograms, sum_g, sum_h, value
                    )
                )
                predictors.append([n_features * bins_per_feature])
                seconds.append(dt)

        (ops,) = fit_marginal_ops_per_item(np.array(predictors), np.array(seconds))
        values[shape] = ops

    report("find_node_split (per feature x bin)", values)


# -----------------------------------------------------------------------
# split_indices, from Splitter.split_indices.
# -----------------------------------------------------------------------
def bench_split_indices():
    n_features = 8
    values = {}
    for shape in ("plain", "with missing values", "categorical"):
        predictors = []
        seconds = []
        for n_samples in (2_000, 20_000, 200_000, 1_000_000):
            splitter, histograms, sum_g, sum_h, value = _make_splitter_and_histograms(
                n_samples, n_features, shape=shape
            )
            split_info = splitter.find_node_split(
                n_samples, histograms, sum_g, sum_h, value
            )

            def run():
                splitter.split_indices(
                    split_info, np.arange(n_samples, dtype=np.uint32)
                )

            dt = best_of_seconds(run)
            predictors.append([n_samples])
            seconds.append(dt)

        (ops,) = fit_marginal_ops_per_item(np.array(predictors), np.array(seconds))
        values[shape] = ops

    report("split_indices (per sample)", values)


# -----------------------------------------------------------------------
# Prediction (binned and raw), from TreePredictor.predict_binned/.predict.
# -----------------------------------------------------------------------
def bench_predict():
    n_bins = 256
    n_train, n_features = 5_000, 10
    X_binned_train = np.asfortranarray(
        rng.integers(0, n_bins - 1, size=(n_train, n_features)), dtype=X_BINNED_DTYPE
    )
    gradients = rng.standard_normal(n_train).astype(G_H_DTYPE)
    hessians = np.ones(1, dtype=G_H_DTYPE)
    grower = TreeGrower(
        X_binned_train,
        gradients,
        hessians,
        n_bins=n_bins,
        max_leaf_nodes=31,
        min_samples_leaf=5,
        n_threads=1,
    )
    grower.grow()
    predictor = grower.make_predictor(
        binning_thresholds=[np.arange(n_bins - 1) for _ in range(n_features)]
    )

    known_cat_bitsets = np.zeros((0, 8), dtype=np.uint32)
    f_idx_map = np.zeros(0, dtype=np.uint32)

    values = {}
    for label in ("binned", "raw"):
        predictors = []
        seconds = []
        for n_samples in (2_000, 20_000, 200_000, 1_000_000):
            # Input data is generated once, outside the timed closure, so
            # its own (non-negligible, and very different between the two
            # dtypes/generators) cost doesn't contaminate the measurement.
            if label == "binned":
                X = np.asfortranarray(
                    rng.integers(0, n_bins - 1, size=(n_samples, n_features)),
                    dtype=X_BINNED_DTYPE,
                )
                run_once = lambda: predictor.predict_binned(X, n_bins - 1, 1)
            else:
                X = rng.random((n_samples, n_features))
                run_once = lambda: predictor.predict(X, known_cat_bitsets, f_idx_map, 1)

            dt = best_of_seconds(run_once)
            predictors.append([n_samples])
            seconds.append(dt)

        (ops,) = fit_marginal_ops_per_item(np.array(predictors), np.array(seconds))
        values[label] = ops

    report("predict (per sample)", values)


# -----------------------------------------------------------------------
# update_raw_predictions, from _gradient_boosting._update_raw_predictions.
# Its inner per-sample loop is a `cdef inline` helper with no Python-visible
# entry point of its own, so this goes through the public
# `_update_raw_predictions` wrapper instead, with a minimal duck-typed
# stand-in for the `grower` object it expects (only `.splitter.partition`
# and `.finalized_leaves[*].{partition_start,partition_stop,value}` are
# actually read).
# -----------------------------------------------------------------------
def bench_update_raw_predictions():
    from sklearn.ensemble._hist_gradient_boosting._gradient_boosting import (
        _update_raw_predictions,
    )

    class _FakeSplitter:
        def __init__(self, partition):
            self.partition = partition

    class _FakeLeaf:
        def __init__(self, start, stop, value):
            self.partition_start = start
            self.partition_stop = stop
            self.value = value

    class _FakeGrower:
        def __init__(self, partition, leaves):
            self.splitter = _FakeSplitter(partition)
            self.finalized_leaves = leaves

    predictors = []
    seconds = []
    for n_samples in (2_000, 20_000, 200_000, 1_000_000, 5_000_000):
        partition = np.arange(n_samples, dtype=np.uint32)
        # A handful of leaves, as a real fitted tree would have; the cost is
        # linear in n_samples regardless of how many leaves it's split into.
        n_leaves = 8
        edges = np.linspace(0, n_samples, n_leaves + 1, dtype=np.uint32)
        leaves = [_FakeLeaf(edges[i], edges[i + 1], float(i)) for i in range(n_leaves)]
        grower = _FakeGrower(partition, leaves)
        raw_predictions = np.zeros(n_samples, dtype=Y_DTYPE)

        dt = best_of_seconds(
            lambda: _update_raw_predictions(raw_predictions, grower, 1)
        )
        predictors.append([n_samples])
        seconds.append(dt)

    (ops,) = fit_marginal_ops_per_item(np.array(predictors), np.array(seconds))
    report("update_raw_predictions (per sample)", {"": ops})


# -----------------------------------------------------------------------
# Sanity check: does the measured ops_per_item, combined with this
# machine's calibrated MIN_INSTRUCTIONS_PER_THREAD, predict a crossover
# workload size that matches where threading the real code path actually
# starts paying off?
# -----------------------------------------------------------------------
def validate_crossover():
    # Measured separately above (bench_find_node_split, shape="plain");
    # hardcode a representative value here so this section can run
    # standalone if the earlier ones are skipped.
    ops_per_item = 5

    print(
        f"\nfind_node_split crossover check (ops_per_item={ops_per_item}, n_bins=256):"
    )
    for n_threads in (2, 4, 8):
        # `_use_threads_for_workload` parallelizes iff
        #   t_serial / n_threads + overhead(n_threads) < t_serial
        # i.e. (
        #   solving for t_serial,
        #   then workload_size = t_serial * INSTRUCTIONS_PER_SECOND
        # ):
        #   workload_size > INSTRUCTIONS_PER_SECOND * overhead(n_threads)
        #                   * n_threads / (n_threads - 1)
        overhead = _omp_prange_dispatch_overhead(n_threads)
        predicted_crossover_workload = (
            INSTRUCTIONS_PER_SECOND * overhead * n_threads / (n_threads - 1)
        )
        predicted_crossover = predicted_crossover_workload / (ops_per_item * 256)
        print(
            f"  n_threads={n_threads}: dispatch overhead={overhead * 1e6:.2f}us, "
            f"predicted crossover at ~{predicted_crossover:.0f} features"
        )

        for n_features_test in (
            max(1, int(predicted_crossover * 0.5)),
            max(1, int(predicted_crossover * 2)),
        ):
            splitter, histograms, sum_g, sum_h, value = _make_splitter_and_histograms(
                5_000, n_features_test, n_bins=256, shape="plain"
            )
            splitter.n_threads = 1
            t_seq = best_of_seconds(
                lambda: splitter.find_node_split(
                    5_000, histograms, sum_g, sum_h, value
                ),
                n_repeats=20,
            )
            splitter.n_threads = n_threads
            t_par = best_of_seconds(
                lambda: splitter.find_node_split(
                    5_000, histograms, sum_g, sum_h, value
                ),
                n_repeats=20,
            )
            speedup = t_seq / t_par
            print(
                f"    n_features={n_features_test:>4}: "
                f"sequential={t_seq * 1e3:.3f}ms, "
                f"n_threads={n_threads}={t_par * 1e3:.3f}ms, "
                f"speedup={speedup:.2f}x"
            )


def main():
    print("Calibrating HGB 'ops per item' constants on this machine...")
    bench_histogram_brute()
    bench_histogram_subtraction()
    bench_find_node_split()
    bench_split_indices()
    bench_predict()
    bench_update_raw_predictions()
    validate_crossover()


if __name__ == "__main__":
    main()
