from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn._pprint import _Formatter
from sklearn.utils.testing import assert_raise_message
from sklearn.base import BaseEstimator


def test_indent_est_param():

    pipeline = make_pipeline(StandardScaler(), LogisticRegression())

    f = _Formatter(indent_est='name')
    expected_repr = """
Pipeline(
         steps=[('standardscaler',
                 StandardScaler(copy=True, with_mean=True, with_std=True))
                ('logisticregression',
                 LogisticRegression(C=1.0, class_weight=None, dual=False,
                                    fit_intercept=True, intercept_scaling=1,
                                    max_iter=100, multi_class='ovr', n_jobs=1,
                                    penalty='l2', random_state=None,
                                    solver='liblinear', tol=0.0001, verbose=0,
                                    warm_start=False))],
         memory=None)"""
    expected_repr = expected_repr[1:]  # Remove first \n
    assert expected_repr == f(pipeline)

    f = _Formatter(indent_est='step')
    expected_repr = """
Pipeline(
    steps=[('standardscaler',
            StandardScaler(copy=True, with_mean=True, with_std=True))
           ('logisticregression',
            LogisticRegression(C=1.0, class_weight=None, dual=False,
                fit_intercept=True, intercept_scaling=1, max_iter=100,
                multi_class='ovr', n_jobs=1, penalty='l2', random_state=None,
                solver='liblinear', tol=0.0001, verbose=0, warm_start=False))],
    memory=None)"""
    expected_repr = expected_repr[1:]  # Remove first \n
    assert expected_repr == f(pipeline)

    f = _Formatter(indent_est='wrong')
    assert_raise_message(ValueError, "Invalid indent_est parameter.", f,
                         pipeline)


def test_changed_only_param():

    pipeline = make_pipeline(StandardScaler(), LogisticRegression(C=999))
    f = _Formatter(changed_only=True)
    expected_repr = """
Pipeline(
    steps=[('standardscaler', StandardScaler())
           ('logisticregression', LogisticRegression(C=999))],
    )"""
    expected_repr = expected_repr[1:]  # Remove first \n
    assert expected_repr == f(pipeline)


def test_visual():
    # Not really tests. Just here to make reviewing easier.

    # 84 chararacters long estimator name
    class Veeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeerylong(BaseEstimator): # noqa
        def __init__(self, something=True):
            pass

        def fit(X):
            pass

    VeryLongName =  Veeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeerylong  # noqa

    f = _Formatter()

    expected_repr = """
Veeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeerylong(
    something=None)"""
    expected_repr = expected_repr[1:]  # Remove first \n
    assert expected_repr == f(VeryLongName())

    pipeline = make_pipeline(StandardScaler(), VeryLongName())
    expected_repr = """
Pipeline(
    steps=[('standardscaler',
            StandardScaler(copy=True, with_mean=True, with_std=True))
           (
            'veeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeerylong',
            Veeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeerylong(
                something=None))],
    memory=None)"""
    expected_repr = expected_repr[1:]  # Remove first \n
    assert expected_repr == f(pipeline)


    estimator = RFE(RFE(RFE(RFE(RFE(RFE(RFE(LogisticRegression)))))))
    expected_repr = """
RFE(
    estimator=RFE(
        estimator=RFE(
            estimator=RFE(
                estimator=RFE(
                    estimator=RFE(
                        estimator=RFE(
                            estimator=<class 'sklearn.linear_model.logistic.LogisticRegression'>,
                            n_features_to_select=None, step=1, verbose=0),
                        n_features_to_select=None, step=1, verbose=0),
                    n_features_to_select=None, step=1, verbose=0),
                n_features_to_select=None, step=1, verbose=0),
            n_features_to_select=None, step=1, verbose=0),
        n_features_to_select=None, step=1, verbose=0),
    n_features_to_select=None, step=1, verbose=0)"""
    expected_repr = expected_repr[1:]  # Remove first \n
    assert expected_repr == f(estimator)
