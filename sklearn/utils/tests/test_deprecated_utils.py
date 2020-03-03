import pytest
import types
import numpy as np
import warnings

from sklearn.dummy import DummyClassifier
from sklearn.utils import all_estimators
from sklearn.utils.estimator_checks import choose_check_classifiers_labels
from sklearn.utils.estimator_checks import NotAnArray
from sklearn.utils.estimator_checks import enforce_estimator_tags_y
from sklearn.utils.estimator_checks import is_public_parameter
from sklearn.utils.estimator_checks import pairwise_estimator_convert_X
from sklearn.utils.estimator_checks import set_checking_parameters
from sklearn.utils.optimize import newton_cg
from sklearn.utils.random import random_choice_csc
from sklearn.utils import safe_indexing


# This file tests the utils that are deprecated


# TODO: remove in 0.24
def test_choose_check_classifiers_labels_deprecated():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        choose_check_classifiers_labels(None, None, None)


# TODO: remove in 0.24
def test_enforce_estimator_tags_y():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        enforce_estimator_tags_y(DummyClassifier(), np.array([0, 1]))


# TODO: remove in 0.24
def test_notanarray():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        NotAnArray([1, 2])


# TODO: remove in 0.24
def test_is_public_parameter():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        is_public_parameter('hello')


# TODO: remove in 0.24
def test_pairwise_estimator_convert_X():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        pairwise_estimator_convert_X([[1, 2]], DummyClassifier())


# TODO: remove in 0.24
def test_set_checking_parameters():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        set_checking_parameters(DummyClassifier())


# TODO: remove in 0.24
def test_newton_cg():
    rng = np.random.RandomState(0)
    A = rng.normal(size=(10, 10))
    x0 = np.ones(10)

    def func(x):
        Ax = A.dot(x)
        return .5 * (Ax).dot(Ax)

    def grad(x):
        return A.T.dot(A.dot(x))

    def grad_hess(x):
        return grad(x), lambda x: A.T.dot(A.dot(x))

    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        newton_cg(grad_hess, func, grad, x0)


# TODO: remove in 0.24
def test_random_choice_csc():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        random_choice_csc(10, [[2]])


# TODO: remove in 0.24
def test_safe_indexing():
    with pytest.warns(FutureWarning,
                      match="removed in version 0.24"):
        safe_indexing([1, 2], 0)


# TODO: remove in 0.24
def test_partial_dependence_no_shadowing():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/15842
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        from sklearn.inspection.partial_dependence import partial_dependence as _  # noqa

        # Calling all_estimators() also triggers a recursive import of all
        # submodules, including deprecated ones.
        all_estimators()

    from sklearn.inspection import partial_dependence
    assert isinstance(partial_dependence, types.FunctionType)


# TODO: remove in 0.24
def test_dict_learning_no_shadowing():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/15842
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        from sklearn.decomposition.dict_learning import dict_learning as _  # noqa

        # Calling all_estimators() also triggers a recursive import of all
        # submodules, including deprecated ones.
        all_estimators()

    from sklearn.decomposition import dict_learning
    assert isinstance(dict_learning, types.FunctionType)
