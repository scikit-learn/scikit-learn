from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

from .common import Benchmark, Estimator, Predictor
from .datasets import (_20newsgroups_highdim_dataset,
                       _20newsgroups_lowdim_dataset)
from .utils import make_gen_classif_scorers


class RandomForestClassifierBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for RandomForestClassifier.
    """

    param_names = ['representation', 'n_jobs']
    params = (['dense', 'sparse'], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        representation, n_jobs = params

        if representation == 'sparse':
            data = _20newsgroups_highdim_dataset()
        else:
            data = _20newsgroups_lowdim_dataset()

        return data

    def make_estimator(self, params):
        representation, n_jobs = params

        n_estimators = 500 if Benchmark.data_size == 'large' else 100

        estimator = RandomForestClassifier(n_estimators=n_estimators,
                                           min_samples_split=10,
                                           max_features='log2',
                                           n_jobs=n_jobs,
                                           random_state=0)

        return estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)


class GradientBoostingClassifierBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for GradientBoostingClassifier.
    """

    param_names = ['representation']
    params = (['dense', 'sparse'],)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        representation, = params

        if representation == 'sparse':
            data = _20newsgroups_highdim_dataset()
        else:
            data = _20newsgroups_lowdim_dataset()

        return data

    def make_estimator(self, params):
        representation, = params

        n_estimators = 100 if Benchmark.data_size == 'large' else 10

        estimator = GradientBoostingClassifier(n_estimators=n_estimators,
                                               max_features='log2',
                                               subsample=0.5,
                                               random_state=0)

        return estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)
