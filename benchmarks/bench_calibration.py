import numpy as np
import itertools
import time

from sklearn.calibration import CalibratedClassifierCV
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier


class CalibratedClassifierCVBenchmark(object):
    """
    Utilities to benchmark CalibratedClassifierCV.
    """

    def __init__(self):
        self.X, self.y = make_classification(n_samples=10000, n_features=10,
                                             n_informative=3, n_redundant=0,
                                             random_state=42)

    def to_benchmark(self, *args, **kwargs):
        raise NotImplementedError

    def benchmark(self, func, n_trials = 10):
        """
        Return a proxy function to run benchmark
        on the original one.

        :param func: function
        :param n_trials: number of execution for the benchmark
        :return:
        """

        def proxy(*args, **key_args):

            times = []
            for _ in range(n_trials):
                t1 = time.time()
                r = func(*args, **key_args)
                t2 = time.time()
                times.append(t2 - t1)

            mean = np.mean(times)
            std = np.std(times)
            print("{} trials: {:.4f} ± {:.4f} s".format(n_trials, mean, std))
            return r

        return proxy

    def _run_benchmark(self):
        for values in itertools.product(*self.params):
            args = dict(zip(self.param_names, values))
            print(self.__class__.__name__, args, end=" ")
            self.benchmark(self.to_benchmark)(**args)
        print()

    def run_all(self):
        """
        Run all the benchmark defined in classes
        extending CalibratedClassifierCVBenchmark.

        :return:
        """
        for benchmark_class in self.__class__.__subclasses__():
            benchmark_class()._run_benchmark()
            print()


class BenchmarkNJobsAlgo(CalibratedClassifierCVBenchmark):

    ESTIMATORS = {
        'LogisticReg.': LogisticRegression(),
        'CART': DecisionTreeClassifier(),
        'ExtraTrees': ExtraTreesClassifier(),
        'RandomForest': RandomForestClassifier(),
    }

    params = [list(ESTIMATORS), [1, 2, 4, 8]]

    param_names = ['estimator_name', 'n_jobs']

    def to_benchmark(self, estimator_name, n_jobs):
        clf = self.ESTIMATORS[estimator_name]
        clf_c = CalibratedClassifierCV(clf, n_jobs=n_jobs)
        clf_c.fit(self.X, self.y)


class BenchmarkCV(CalibratedClassifierCVBenchmark):

    params = [[2, 3, 4, 5]]

    param_names = ['cv']

    def to_benchmark(self, cv):
        clf = LogisticRegression()
        clf_c = CalibratedClassifierCV(clf, cv=cv)
        clf_c.fit(self.X, self.y)


class BenchmarkNJobsCV(CalibratedClassifierCVBenchmark):

    params = [[2, 3, 4, 5], [1, 2, 4, 8]]

    param_names = ['cv', 'n_jobs']

    def to_benchmark(self, cv, n_jobs):
        clf = LogisticRegression()
        clf_c = CalibratedClassifierCV(clf, n_jobs=n_jobs, cv=cv)
        clf_c.fit(self.X, self.y)


if __name__ == "__main__":
    CalibratedClassifierCVBenchmark().run_all()