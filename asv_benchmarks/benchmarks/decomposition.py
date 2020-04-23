from sklearn.decomposition import (PCA, DictionaryLearning,
                                   MiniBatchDictionaryLearning)

from .common import Benchmark, Estimator, Transformer
from .datasets import _olivetti_faces_dataset, _mnist_dataset
from .utils import make_pca_scorers, make_dict_learning_scorers


class PCA_bench(Transformer, Estimator, Benchmark):
    """
    Benchmarks for PCA.
    """

    param_names = ['svd_solver']
    params = (['full', 'arpack', 'randomized'],)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        svd_solver, = params

        data = _mnist_dataset()

        estimator = PCA(n_components=32,
                        svd_solver=svd_solver,
                        random_state=0)

        return data, estimator

    def make_scorers(self):
        make_pca_scorers(self)


class DictionaryLearning_bench(Transformer, Estimator, Benchmark):
    """
    Benchmarks for DictionaryLearning.
    """

    param_names = ['fit_algorithm', 'n_jobs']
    params = (['lars', 'cd'], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        fit_algorithm, n_jobs = params

        data = _olivetti_faces_dataset()

        estimator = DictionaryLearning(n_components=15,
                                       fit_algorithm=fit_algorithm,
                                       alpha=0.1,
                                       max_iter=20,
                                       tol=1e-16,
                                       random_state=0,
                                       n_jobs=n_jobs)

        return data, estimator

    def make_scorers(self):
        make_dict_learning_scorers(self)


class MiniBatchDictionaryLearning_bench(Transformer, Estimator, Benchmark):
    """
    Benchmarks for MiniBatchDictionaryLearning
    """

    param_names = ['fit_algorithm', 'n_jobs']
    params = (['lars', 'cd'], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        fit_algorithm, n_jobs = params

        data = _olivetti_faces_dataset()

        estimator = MiniBatchDictionaryLearning(n_components=15,
                                                fit_algorithm=fit_algorithm,
                                                alpha=0.1,
                                                batch_size=3,
                                                random_state=0,
                                                n_jobs=n_jobs)

        return data, estimator

    def make_scorers(self):
        make_dict_learning_scorers(self)
