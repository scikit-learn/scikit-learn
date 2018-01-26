from functools import partial
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.decomposition import _update_dict
from sklearn.decomposition.dict_learning_utils import _general_update_dict
from sklearn.linear_model.coordescendant import (L2INF_CONSTRAINT,
                                                 coordescendant)
from sklearn.linear_model.python_wrappers import _py_proj_l2
from .coordescendant_slow import coordescendant_slow
from .test_coordescendant import _make_test_data


def test_sklearn_update_dict():
    for n_samples in [4, 10]:
        for n_features in [n_samples // 2, n_samples * 2]:
            for n_targets in [n_features // 2, n_features * 2]:
                X, Y, gram, cov, W = _make_test_data(
                    n_samples, n_features, n_targets, real=True)
                for zero in [True, False]:
                    gram_ = gram.copy()
                    cov_ = cov.copy()
                    if zero:
                        cov_ *= 0.
                        gram_ *= 0.

                W1 = W.copy()
                _update_dict(W1.T, cov_.T, gram_.T, random_state=0)

                for penalty_model in [_py_proj_l2, L2INF_CONSTRAINT]:
                    for solver in [coordescendant, coordescendant_slow]:
                        W2 = W.copy()
                        updater = partial(
                            _general_update_dict, penalty_model=penalty_model,
                            emulate_sklearn=True, solver=solver, random_state=0)
                        _update_dict(W2.T, cov_.T, gram_.T, random_state=0,
                                     updater=updater)
                    assert_array_almost_equal(W1, W2)
