from sklearn.neighbors import KNeighborsClassifier, RadiusNeighborsRegressor

from .common import Benchmark, Estimator, Predictor
from .datasets import _20newsgroups_lowdim_dataset, _synth_regression_dataset
from .utils import make_gen_classif_scorers, make_gen_reg_scorers


class KNeighborsClassifierBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for KNeighborsClassifier.
    """

    param_names = ["algorithm", "dimension", "n_jobs"]
    params = (["brute", "kd_tree", "ball_tree"], ["low", "high"], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        algorithm, dimension, n_jobs = params

        if Benchmark.data_size == "large":
            n_components = 40 if dimension == "low" else 200
        else:
            n_components = 10 if dimension == "low" else 50

        data = _20newsgroups_lowdim_dataset(n_components=n_components)

        return data

    def make_estimator(self, params):
        algorithm, dimension, n_jobs = params

        estimator = KNeighborsClassifier(algorithm=algorithm, n_jobs=n_jobs)

        return estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)


class RadiusNeighborsRegressorBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for RadiusNeighborsRegressor
    """

    param_names = ["algorithm", "dimension", "n_jobs"]
    params = (
        ["brute", "kd_tree", "ball_tree"],
        ["very-low", "low", "high"],
        Benchmark.n_jobs_vals,
    )
    dim_to_number = {"very-low": 3, "low": 20, "high": 200}

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        algorithm, dimension, n_jobs = params

        n_features = self.dim_to_number[dimension]

        data = _synth_regression_dataset(n_samples=10000, n_features=n_features)

        return data

    def make_estimator(self, params):
        algorithm, dimension, n_jobs = params

        estimator = RadiusNeighborsRegressor(
            algorithm=algorithm, n_jobs=n_jobs, radius=0.05
        )

        return estimator

    def make_scorers(self):
        make_gen_reg_scorers(self)
