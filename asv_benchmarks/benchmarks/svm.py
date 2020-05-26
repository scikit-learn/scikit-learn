from sklearn.svm import SVC

from .common import Benchmark, Estimator, Predictor
from .datasets import _synth_classification_dataset
from .utils import optimal_cache_size, make_gen_classif_scorers


class SVCBenchmark(Predictor, Estimator, Benchmark):
    """Benchmarks for SVC."""

    param_names = ['kernel']
    params = (['linear', 'poly', 'rbf', 'sigmoid'],)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        kernel, = params

        data = _synth_classification_dataset()

        estimator = SVC(cache_size=optimal_cache_size(data[0].shape[1]),
                        max_iter=100,
                        tol=1e-16,
                        kernel=kernel,
                        random_state=0,
                        gamma='scale')

        return data, estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)
