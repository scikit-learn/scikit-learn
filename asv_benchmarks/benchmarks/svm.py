from sklearn.svm import SVC

from .common import Benchmark, Estimator, Predictor
from .datasets import _synth_classification_dataset
from .utils import make_gen_classif_scorers


class SVCBenchmark(Predictor, Estimator, Benchmark):
    """Benchmarks for SVC."""

    param_names = ['kernel']
    params = (['linear', 'poly', 'rbf', 'sigmoid'],)

    def setup_cache(self):
        super().setup_cache()

    def make_data(self, params):
        return _synth_classification_dataset()

    def make_estimator(self, params):
        kernel, = params

        estimator = SVC(max_iter=100,
                        tol=1e-16,
                        kernel=kernel,
                        random_state=0,
                        gamma='scale')

        return estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)
