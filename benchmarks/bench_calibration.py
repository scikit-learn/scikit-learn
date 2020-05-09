import numpy as np
import itertools
import time
import seaborn as sns
import gc

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.calibration import CalibratedClassifierCV
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier


class CalibratedClassifierCVBenchmark(object):
    """
    Utilities to benchmark CalibratedClassifierCV.

    :param n_samples: number of samples of the dataset to use for the benchmark
    :param n_trials: number of trials of each test for the benchmark
    """

    def __init__(self, n_samples=100000, n_features=10, n_trials=20):
        self.n_samples = n_samples
        self.n_features = n_features
        self.n_trials = n_trials
        self.X, self.y = make_classification(
            n_samples=self.n_samples,
            n_features=self.n_features,
            n_informative=3,
            n_redundant=0,
            random_state=42,
        )

    def to_benchmark(self, *args, **kwargs):
        raise NotImplementedError

    def benchmark(self, func):
        """
        Return a proxy function to run benchmark
        on the original one.

        :param func: function
        :return:
        """

        def proxy(*args, **key_args):
            times = []
            gc.collect()
            for _ in range(self.n_trials):
                t1 = time.time()
                _ = func(*args, **key_args)
                t2 = time.time()
                times.append(t2 - t1)

            mean = np.mean(times)
            std = np.std(times)
            print("{} trials: {:.4f} ± {:.4f} s"
                  .format(self.n_trials, mean, std))
            print(times)
            return args, key_args, times

        return proxy

    def _run_benchmark(self):
        res = []
        for values in itertools.product(*self.params):
            args = dict(zip(self.param_names, values))
            print(self.__class__.__name__, args, end=" ")
            res.append(self.benchmark(self.to_benchmark)(**args))
        print()
        return res

    def _plot_res(self, benchmark_class, res):
        plt.figure(figsize=(21, 13))

        df = pd.DataFrame()
        for args, key_args, times in res:
            d = pd.DataFrame(dict(times=times))
            key_args_string = "\n".join(list(
                map(lambda x: f"{x[0]}: {x[1]}", key_args.items())))

            d['args'] = key_args_string
            df = pd.concat([df, d])

        ax = sns.boxplot(x="args", y="times", data=df, whis=0.4)

        ax.set_ylabel("Execution time (in sec.)")
        ax.xaxis.set_tick_params(rotation=45)
        plt.title(f"{benchmark_class.__doc__ }"
                  f" — dataset-shape: ({self.n_samples}, "
                  f"{self.n_features}) — {self.n_trials} trials")
        ax.yaxis.grid(True)

        plt.savefig(f"{benchmark_class.__name__.lower()}.png")

        plt.show()
        plt.close()

    def run_all(self):
        """
        Run all the benchmark defined in classes
        extending CalibratedClassifierCVBenchmark.

        :return:
        """
        for benchmark_class in self.__class__.__subclasses__():
            res = benchmark_class()._run_benchmark()
            self._plot_res(benchmark_class, res)
            print()


class BenchmarkNJobsSingleThreadAlgo(CalibratedClassifierCVBenchmark):
    """ Time vs single-threaded algo. and threads numbers."""

    ESTIMATORS = {
        "LogisticReg.": LogisticRegression(),
        "CART": DecisionTreeClassifier(),
    }

    params = [list(ESTIMATORS), [1, 2, 4, 8]]

    param_names = ["estimator_name", "n_jobs"]

    def to_benchmark(self, estimator_name, n_jobs):
        clf = self.ESTIMATORS[estimator_name]
        clf_c = CalibratedClassifierCV(clf, n_jobs=n_jobs)
        clf_c.fit(self.X, self.y)


class BenchmarkNJobsMultiThreadAlgo(CalibratedClassifierCVBenchmark):
    """ Time vs multi-threaded algo. and threads numbers."""

    ESTIMATORS = {
        "ExtraTrees": ExtraTreesClassifier(),
        "RandomForest": RandomForestClassifier(),
    }

    params = [list(ESTIMATORS), [1, 2, 4, 8]]

    param_names = ["estimator_name", "n_jobs"]

    def to_benchmark(self, estimator_name, n_jobs):
        clf = self.ESTIMATORS[estimator_name]
        clf_c = CalibratedClassifierCV(clf, n_jobs=n_jobs)
        clf_c.fit(self.X, self.y)


class BenchmarkCV(CalibratedClassifierCVBenchmark):
    """ Time vs number of folds."""

    params = [[2, 3, 4, 5]]

    param_names = ["cv"]

    def to_benchmark(self, cv):
        clf = LogisticRegression()
        clf_c = CalibratedClassifierCV(clf, cv=cv)
        clf_c.fit(self.X, self.y)


class BenchmarkNJobsCV(CalibratedClassifierCVBenchmark):
    """ Time vs number of threads and number of folds."""

    params = [[2, 3, 4, 5], [1, 2, 4, 8]]

    param_names = ["cv", "n_jobs"]

    def to_benchmark(self, cv, n_jobs):
        clf = LogisticRegression()
        clf_c = CalibratedClassifierCV(clf, n_jobs=n_jobs, cv=cv)
        clf_c.fit(self.X, self.y)


if __name__ == "__main__":
    CalibratedClassifierCVBenchmark().run_all()
