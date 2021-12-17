from sklearn.decomposition import PCA, DictionaryLearning, MiniBatchDictionaryLearning

from .common import Benchmark, Estimator, Transformer
from .datasets import _olivetti_faces_dataset, _mnist_dataset
from .utils import make_pca_scorers, make_dict_learning_scorers


class PCABenchmark(Transformer, Estimator, Benchmark):
    """
    Benchmarks for PCA.
    """

    param_names = ["svd_solver"]
    params = (["full", "arpack", "randomized"],)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        return _mnist_dataset()

    def make_estimator(self, params):
        (svd_solver,) = params

        estimator = PCA(n_components=32, svd_solver=svd_solver, random_state=0)

        return estimator

    def make_scorers(self):
        make_pca_scorers(self)


class DictionaryLearningBenchmark(Transformer, Estimator, Benchmark):
    """
    Benchmarks for DictionaryLearning.
    """

    param_names = ["fit_algorithm", "n_jobs"]
    params = (["lars", "cd"], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        return _olivetti_faces_dataset()

    def make_estimator(self, params):
        fit_algorithm, n_jobs = params

        estimator = DictionaryLearning(
            n_components=15,
            fit_algorithm=fit_algorithm,
            alpha=0.1,
            max_iter=20,
            tol=1e-16,
            random_state=0,
            n_jobs=n_jobs,
        )

        return estimator

    def make_scorers(self):
        make_dict_learning_scorers(self)


class MiniBatchDictionaryLearningBenchmark(Transformer, Estimator, Benchmark):
    """
    Benchmarks for MiniBatchDictionaryLearning
    """

    param_names = ["fit_algorithm", "n_jobs"]
    params = (["lars", "cd"], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        return _olivetti_faces_dataset()

    def make_estimator(self, params):
        fit_algorithm, n_jobs = params

        estimator = MiniBatchDictionaryLearning(
            n_components=15,
            fit_algorithm=fit_algorithm,
            alpha=0.1,
            batch_size=3,
            random_state=0,
            n_jobs=n_jobs,
        )

        return estimator

    def make_scorers(self):
        make_dict_learning_scorers(self)
