from __future__ import annotations

import pickle

import numpy as np

from sklearn.preprocessing import TargetEncoder

from .common import Benchmark, clear_tmp, get_estimator_path


class TargetEncoderTransformSmallBatch(Benchmark):
    """
    Measure TargetEncoder.fit and .transform latency and peak memory across
    batch sizes around the small-batch threshold.

    Parameters:
      - dtype: int64 vs object (object injects None/np.nan at transform time)
      - target_type: binary vs multiclass
      - n_rows: batch size (only affects transform benchmarks; fit always uses
        the full training set)

    Transform timings/memory:
      - time_transform_default / peakmem_transform_default:
            whatever branch TargetEncoder.transform chooses
      - time_transform_forced_vectorized / peakmem_transform_forced_vectorized:
            force the vectorized branch regardless of n_rows

    Fit timings/memory:
      - time_fit / peakmem_fit:
            full fit including _init_int_lookup (current branch behaviour)
      - time_fit_no_lut / peakmem_fit_no_lut:
            fit with _init_int_lookup replaced by a no-op, isolating its cost
    """

    param_names = ["dtype", "target_type", "n_rows"]
    params = (
        ["int64", "object"],
        ["binary", "multiclass"],
        [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096],
    )

    _n_cats: int = 100_000

    def setup_cache(self):
        clear_tmp()

        n_fit = 500_000 if Benchmark.data_size == "large" else 200_000
        n_cats = 100_000
        type(self)._n_cats = n_cats

        rng = np.random.default_rng(0)

        X_fit_int = rng.integers(0, n_cats, size=(n_fit, 1), dtype=np.int64)
        X_fit_obj = X_fit_int.astype(object)

        # Exercise missing handling as well:
        X_fit_obj[:100, 0] = None
        X_fit_obj[100:200, 0] = np.nan

        y_binary = rng.integers(0, 2, size=n_fit, dtype=np.int64)
        y_multi = rng.integers(0, 4, size=n_fit, dtype=np.int64)

        for dtype in self.params[0]:
            X_fit = X_fit_int if dtype == "int64" else X_fit_obj

            for target_type in self.params[1]:
                y = y_binary if target_type == "binary" else y_multi

                enc = TargetEncoder(target_type=target_type, random_state=0).fit(
                    X_fit, y
                )

                est_path = get_estimator_path(
                    self,
                    Benchmark.save_dir,
                    (dtype, target_type),
                    Benchmark.save_estimators,
                )
                with est_path.open("wb") as f:
                    pickle.dump((enc, X_fit.copy(), y.copy()), f)

    def setup(self, dtype: str, target_type: str, n_rows: int):
        est_path = get_estimator_path(
            self,
            Benchmark.save_dir,
            (dtype, target_type),
            Benchmark.save_estimators,
        )

        with est_path.open("rb") as f:
            self.enc, self.X_fit, self.y = pickle.load(f)

        with est_path.open("rb") as f:
            self.enc_vec, _, _ = pickle.load(f)

        if hasattr(self.enc_vec, "_small_batch_threshold"):
            self.enc_vec._small_batch_threshold = -1  # always take vectorized path

        self.target_type = target_type
        rng = np.random.default_rng(0)
        n_cats = type(self)._n_cats
        max_rows = max(self.params[2])

        X = rng.integers(0, 2 * n_cats, size=(max_rows, 1), dtype=np.int64)[
            :n_rows
        ].copy()

        if dtype == "object":
            X = X.astype(object)
            if n_rows >= 2:
                X[0, 0] = None
                X[-1, 0] = np.nan

        self.X = X

        self.enc.transform(self.X)
        self.enc_vec.transform(self.X)

    def time_transform_default(self, dtype: str, target_type: str, n_rows: int):
        self.enc.transform(self.X)

    def time_transform_forced_vectorized(
        self, dtype: str, target_type: str, n_rows: int
    ):
        self.enc_vec.transform(self.X)

    def peakmem_transform_default(self, dtype: str, target_type: str, n_rows: int):
        self.enc.transform(self.X)

    def peakmem_transform_forced_vectorized(
        self, dtype: str, target_type: str, n_rows: int
    ):
        self.enc_vec.transform(self.X)

    def time_fit(self, dtype: str, target_type: str, n_rows: int):
        TargetEncoder(target_type=self.target_type, random_state=0).fit(
            self.X_fit, self.y
        )

    def time_fit_no_lut(self, dtype: str, target_type: str, n_rows: int):
        enc = TargetEncoder(target_type=self.target_type, random_state=0)
        enc._init_int_lookup = lambda: None
        enc.fit(self.X_fit, self.y)

    def peakmem_fit(self, dtype: str, target_type: str, n_rows: int):
        TargetEncoder(target_type=self.target_type, random_state=0).fit(
            self.X_fit, self.y
        )

    def peakmem_fit_no_lut(self, dtype: str, target_type: str, n_rows: int):
        enc = TargetEncoder(target_type=self.target_type, random_state=0)
        enc._init_int_lookup = lambda: None
        enc.fit(self.X_fit, self.y)
