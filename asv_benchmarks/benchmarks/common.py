import itertools
import json
import os
import pickle
import timeit
from abc import ABC, abstractmethod
from multiprocessing import cpu_count
from pathlib import Path

import numpy as np


def get_from_config():
    """Get benchmarks configuration from the config.json file"""
    current_path = Path(__file__).resolve().parent

    config_path = current_path / "config.json"
    with open(config_path, "r") as config_file:
        config_file = "".join(line for line in config_file if line and "//" not in line)
        config = json.loads(config_file)

    profile = os.getenv("SKLBENCH_PROFILE", config["profile"])

    n_jobs_vals_env = os.getenv("SKLBENCH_NJOBS")
    if n_jobs_vals_env:
        n_jobs_vals = json.loads(n_jobs_vals_env)
    else:
        n_jobs_vals = config["n_jobs_vals"]
    if not n_jobs_vals:
        n_jobs_vals = list(range(1, 1 + cpu_count()))

    cache_path = current_path / "cache"
    cache_path.mkdir(exist_ok=True)
    (cache_path / "estimators").mkdir(exist_ok=True)
    (cache_path / "tmp").mkdir(exist_ok=True)

    save_estimators = os.getenv("SKLBENCH_SAVE_ESTIMATORS", config["save_estimators"])
    save_dir = os.getenv("ASV_COMMIT", "new")[:8]

    if save_estimators:
        (cache_path / "estimators" / save_dir).mkdir(exist_ok=True)

    base_commit = os.getenv("SKLBENCH_BASE_COMMIT", config["base_commit"])

    bench_predict = os.getenv("SKLBENCH_PREDICT", config["bench_predict"])
    bench_transform = os.getenv("SKLBENCH_TRANSFORM", config["bench_transform"])

    return (
        profile,
        n_jobs_vals,
        save_estimators,
        save_dir,
        base_commit,
        bench_predict,
        bench_transform,
    )


def get_estimator_path(benchmark, directory, params, save=False):
    """Get path of pickled fitted estimator"""
    path = Path(__file__).resolve().parent / "cache"
    path = (path / "estimators" / directory) if save else (path / "tmp")

    filename = (
        benchmark.__class__.__name__
        + "_estimator_"
        + "_".join(list(map(str, params)))
        + ".pkl"
    )

    return path / filename


def clear_tmp():
    """Clean the tmp directory"""
    path = Path(__file__).resolve().parent / "cache" / "tmp"
    for child in path.iterdir():
        child.unlink()


class Benchmark(ABC):
    """Abstract base class for all the benchmarks"""

    timer = timeit.default_timer  # wall time
    processes = 1
    timeout = 500

    (
        profile,
        n_jobs_vals,
        save_estimators,
        save_dir,
        base_commit,
        bench_predict,
        bench_transform,
    ) = get_from_config()

    if profile == "fast":
        warmup_time = 0
        repeat = 1
        number = 1
        min_run_count = 1
        data_size = "small"
    elif profile == "regular":
        warmup_time = 1
        repeat = (3, 100, 30)
        data_size = "small"
    elif profile == "large_scale":
        warmup_time = 1
        repeat = 3
        number = 1
        data_size = "large"

    @property
    @abstractmethod
    def params(self):
        pass


class Estimator(ABC):
    """Abstract base class for all benchmarks of estimators"""

    @abstractmethod
    def make_data(self, params):
        """Return the dataset for a combination of parameters"""
        # The datasets are cached using joblib.Memory so it's fast and can be
        # called for each repeat
        pass

    @abstractmethod
    def make_estimator(self, params):
        """Return an instance of the estimator for a combination of parameters"""
        pass

    def skip(self, params):
        """Return True if the benchmark should be skipped for these params"""
        return False

    def setup_cache(self):
        """Pickle a fitted estimator for all combinations of parameters"""
        # This is run once per benchmark class.

        clear_tmp()

        param_grid = list(itertools.product(*self.params))

        for params in param_grid:
            if self.skip(params):
                continue

            estimator = self.make_estimator(params)
            X, _, y, _ = self.make_data(params)

            estimator.fit(X, y)

            est_path = get_estimator_path(
                self, Benchmark.save_dir, params, Benchmark.save_estimators
            )
            with est_path.open(mode="wb") as f:
                pickle.dump(estimator, f)

    def setup(self, *params):
        """Generate dataset and load the fitted estimator"""
        # This is run once per combination of parameters and per repeat so we
        # need to avoid doing expensive operations there.

        if self.skip(params):
            raise NotImplementedError

        self.X, self.X_val, self.y, self.y_val = self.make_data(params)

        est_path = get_estimator_path(
            self, Benchmark.save_dir, params, Benchmark.save_estimators
        )
        with est_path.open(mode="rb") as f:
            self.estimator = pickle.load(f)

        self.make_scorers()

    def time_fit(self, *args):
        self.estimator.fit(self.X, self.y)

    def peakmem_fit(self, *args):
        self.estimator.fit(self.X, self.y)

    def track_train_score(self, *args):
        if hasattr(self.estimator, "predict"):
            y_pred = self.estimator.predict(self.X)
        else:
            y_pred = None
        return float(self.train_scorer(self.y, y_pred))

    def track_test_score(self, *args):
        if hasattr(self.estimator, "predict"):
            y_val_pred = self.estimator.predict(self.X_val)
        else:
            y_val_pred = None
        return float(self.test_scorer(self.y_val, y_val_pred))


class Predictor(ABC):
    """Abstract base class for benchmarks of estimators implementing predict"""

    if Benchmark.bench_predict:

        def time_predict(self, *args):
            self.estimator.predict(self.X)

        def peakmem_predict(self, *args):
            self.estimator.predict(self.X)

        if Benchmark.base_commit is not None:

            def track_same_prediction(self, *args):
                est_path = get_estimator_path(self, Benchmark.base_commit, args, True)
                with est_path.open(mode="rb") as f:
                    estimator_base = pickle.load(f)

                y_val_pred_base = estimator_base.predict(self.X_val)
                y_val_pred = self.estimator.predict(self.X_val)

                return np.allclose(y_val_pred_base, y_val_pred)

    @property
    @abstractmethod
    def params(self):
        pass


class Transformer(ABC):
    """Abstract base class for benchmarks of estimators implementing transform"""

    if Benchmark.bench_transform:

        def time_transform(self, *args):
            self.estimator.transform(self.X)

        def peakmem_transform(self, *args):
            self.estimator.transform(self.X)

        if Benchmark.base_commit is not None:

            def track_same_transform(self, *args):
                est_path = get_estimator_path(self, Benchmark.base_commit, args, True)
                with est_path.open(mode="rb") as f:
                    estimator_base = pickle.load(f)

                X_val_t_base = estimator_base.transform(self.X_val)
                X_val_t = self.estimator.transform(self.X_val)

                return np.allclose(X_val_t_base, X_val_t)

    @property
    @abstractmethod
    def params(self):
        pass
