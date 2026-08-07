import inspect

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C

from .common import Benchmark, Estimator
from .datasets import _synth_regression_dataset
from .utils import make_gen_reg_scorers


class GaussianProcessRegressorBenchmark(Estimator, Benchmark):
    """
    Benchmarks for GaussianProcessRegressor.
    """

    param_names = ["n_restarts", "n_jobs"]
    params = ([0, 1, 2, 4, 8], Benchmark.n_jobs_vals)

    def make_data(self, params):
        return _synth_regression_dataset(n_samples=500, n_features=5, n_informative=5)

    def make_estimator(self, params):
        n_restarts, n_jobs = params

        kernel = C(0.1, (1e-2, 1e2)) * RBF(
            length_scale=1.0, length_scale_bounds=(1e-3, 1e3)
        ) + C(1e-5, (1e-5, 1e2))

        kwargs = {
            "kernel": kernel,
            "n_restarts_optimizer": n_restarts,
            "normalize_y": True,
            "random_state": 0,
        }

        if "n_jobs" in inspect.signature(GaussianProcessRegressor.__init__).parameters:
            kwargs["n_jobs"] = n_jobs

        return GaussianProcessRegressor(**kwargs)

    def make_scorers(self):
        make_gen_reg_scorers(self)
