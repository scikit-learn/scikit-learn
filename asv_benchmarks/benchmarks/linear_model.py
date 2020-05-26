from sklearn.linear_model import (LogisticRegression, Ridge, ElasticNet, Lasso,
                                  LinearRegression, SGDRegressor)

from .common import Benchmark, Estimator, Predictor
from .datasets import (_20newsgroups_highdim_dataset,
                       _20newsgroups_lowdim_dataset,
                       _synth_regression_dataset,
                       _synth_regression_sparse_dataset)
from .utils import make_gen_classif_scorers, make_gen_reg_scorers


class LogisticRegressionBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for LogisticRegression.
    """

    param_names = ['representation', 'solver', 'n_jobs']
    params = (['dense', 'sparse'], ['lbfgs', 'saga'], Benchmark.n_jobs_vals)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, solver, n_jobs = params

        penalty = 'l2' if solver == 'lbfgs' else 'l1'

        if Benchmark.data_size == 'large':
            if representation == 'sparse':
                data = _20newsgroups_highdim_dataset(n_samples=10000)
            else:
                data = _20newsgroups_lowdim_dataset(n_components=1e3)
        else:
            if representation == 'sparse':
                data = _20newsgroups_highdim_dataset(n_samples=2500)
            else:
                data = _20newsgroups_lowdim_dataset()

        estimator = LogisticRegression(solver=solver,
                                       penalty=penalty,
                                       multi_class='multinomial',
                                       tol=0.01,
                                       n_jobs=n_jobs,
                                       random_state=0)

        return data, estimator

    def make_scorers(self):
        make_gen_classif_scorers(self)


class RidgeBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for Ridge.
    """

    param_names = ['representation', 'solver']
    params = (['dense', 'sparse'],
              ['auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag', 'saga'])

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, solver = params

        if representation == 'sparse' and solver == 'svd':
            return

        if representation == 'dense':
            data = _synth_regression_dataset(n_samples=500000, n_features=100)
        else:
            data = _synth_regression_sparse_dataset(n_samples=100000,
                                                    n_features=10000,
                                                    density=0.005)

        estimator = Ridge(solver=solver,
                          fit_intercept=False,
                          random_state=0)

        return data, estimator

    def setup_(self, params):
        representation, solver = params

        if representation == 'sparse' and solver == 'svd':
            raise NotImplementedError

    def make_scorers(self):
        make_gen_reg_scorers(self)


class LinearRegressionBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for Linear Reagression.
    """

    param_names = ['representation']
    params = (['dense', 'sparse'],)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, = params

        if representation == 'dense':
            data = _synth_regression_dataset(n_samples=1000000, n_features=100)
        else:
            data = _synth_regression_sparse_dataset(n_samples=10000,
                                                    n_features=100000,
                                                    density=0.01)

        estimator = LinearRegression()

        return data, estimator

    def make_scorers(self):
        make_gen_reg_scorers(self)


class SGDRegressorBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmark for SGD
    """

    param_names = ['representation']
    params = (['dense', 'sparse'],)

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, = params

        if representation == 'dense':
            data = _synth_regression_dataset(n_samples=100000, n_features=200)
        else:
            data = _synth_regression_sparse_dataset(n_samples=100000,
                                                    n_features=1000,
                                                    density=0.01)

        estimator = SGDRegressor(max_iter=1000,
                                 tol=1e-16,
                                 random_state=0)

        return data, estimator

    def make_scorers(self):
        make_gen_reg_scorers(self)


class ElasticNetBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for ElasticNet.
    """

    param_names = ['representation', 'precompute']
    params = (['dense', 'sparse'], [True, False])

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, precompute = params

        if representation == 'sparse' and precompute is False:
            return

        if representation == 'dense':
            data = _synth_regression_dataset(n_samples=1000000, n_features=100)
        else:
            data = _synth_regression_sparse_dataset(n_samples=50000,
                                                    n_features=5000,
                                                    density=0.01)

        estimator = ElasticNet(precompute=precompute,
                               alpha=0.001,
                               random_state=0)

        return data, estimator

    def setup_(self, params):
        representation, precompute = params

        if representation == 'sparse' and precompute is False:
            raise NotImplementedError

    def make_scorers(self):
        make_gen_reg_scorers(self)


class LassoBenchmark(Predictor, Estimator, Benchmark):
    """
    Benchmarks for Lasso.
    """

    param_names = ['representation', 'precompute']
    params = (['dense', 'sparse'], [True, False])

    def setup_cache(self):
        super().setup_cache()

    def setup_cache_(self, params):
        representation, precompute = params

        if representation == 'sparse' and precompute is False:
            return

        if representation == 'dense':
            data = _synth_regression_dataset(n_samples=1000000, n_features=100)
        else:
            data = _synth_regression_sparse_dataset(n_samples=50000,
                                                    n_features=5000,
                                                    density=0.01)

        estimator = Lasso(precompute=precompute,
                          alpha=0.001,
                          random_state=0)

        return data, estimator

    def setup_(self, params):
        representation, precompute = params

        if representation == 'sparse' and precompute is False:
            raise NotImplementedError

    def make_scorers(self):
        make_gen_reg_scorers(self)
