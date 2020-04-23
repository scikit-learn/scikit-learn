from sklearn.neighbors import KNeighborsClassifier

from .common import Benchmark, Estimator, Predictor
from .datasets import _20newsgroups_lowdim_dataset
from .utils import make_gen_classif_scorers


class KNeighborsClassifier_bench(Predictor, Estimator, Benchmark):
    """
    Benchmarks for KNeighborsClassifier.
    """

    param_names = ['algorithm', 'dimension', 'n_jobs']
    params = (['brute', 'kd_tree', 'ball_tree'],
              ['low', 'high'],
              Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        algorithm, dimension, n_jobs = params

        if Benchmark.data_size == 'large':
            n_components = 40 if dimension == 'low' else 200
        else:
            n_components = 10 if dimension == 'low' else 50

        data = _20newsgroups_lowdim_dataset(n_components=n_components)

        estimator = KNeighborsClassifier(algorithm=algorithm,
                                         n_jobs=n_jobs)

        return data, estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)
