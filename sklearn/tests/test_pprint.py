from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn._pprint import _Formatter
from sklearn.utils.testing import assert_raise_message
from sklearn.base import BaseEstimator
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF


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
    repr_ = f(pipeline)
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    assert max(len(line) for line in repr_.split('\n')) <= f.width

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
    repr_ = f(pipeline)
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    assert max(len(line) for line in repr_.split('\n')) <= f.width

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
    repr_ = f(VeryLongName())
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    # assert max(len(line) for line in repr_.split('\n')) <= f.width

    # Pipeline with very long estimator name
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
    repr_ = f(pipeline)
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    #  assert max(len(line) for line in repr_.split('\n')) <= f.width

    # Deeply nested estimator
    estimator = RFE(RFE(RFE(RFE(RFE(RFE(RFE(LogisticRegression())))))))
    expected_repr = ("""
RFE(
    estimator=RFE(
        estimator=RFE(
            estimator=RFE(
                estimator=RFE(
                    estimator=RFE(
                        estimator=RFE(
                            estimator=LogisticRegression(C=1.0, class_weight=None,"""  # noqa
                     """
                                dual=False, fit_intercept=True,
                                intercept_scaling=1, max_iter=100,
                                multi_class='ovr', n_jobs=1, penalty='l2',
                                random_state=None, solver='liblinear',
                                tol=0.0001, verbose=0, warm_start=False),
                            n_features_to_select=None, step=1, verbose=0),
                        n_features_to_select=None, step=1, verbose=0),
                    n_features_to_select=None, step=1, verbose=0),
                n_features_to_select=None, step=1, verbose=0),
            n_features_to_select=None, step=1, verbose=0),
        n_features_to_select=None, step=1, verbose=0),
    n_features_to_select=None, step=1, verbose=0)""")

    expected_repr = expected_repr[1:]  # Remove first \n
    repr_ = f(estimator)
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    #  assert max(len(line) for line in repr_.split('\n')) <= f.width

    # Grid Search
    param_grid = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                   'C': [1, 10, 100, 1000]},
                  {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

    gs = GridSearchCV(SVC(), param_grid, cv=5)
    expected_repr = """
GridSearchCV(cv=5, error_score='raise-deprecating',
    estimator=SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
        decision_function_shape='ovr', degree=3, gamma='auto_deprecated',
        kernel='rbf', max_iter=-1, probability=False, random_state=None,
        shrinking=True, tol=0.001, verbose=False),
    fit_params=None, iid='warn', n_jobs=1,
    param_grid=[
        {'kernel': ['rbf'], 'gamma': [0.001, 0.0001],
            'C': [1, 10, 100, 1000]},
        {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}],
    pre_dispatch='2*n_jobs', refit=True, return_train_score='warn',
    scoring=None, verbose=0)"""

    expected_repr = expected_repr[1:]  # Remove first \n
    repr_ = f(gs)
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    #  assert max(len(line) for line in repr_.split('\n')) <= f.width

    # Grid Search with a pipeline inside
    pipeline = Pipeline([
        ('reduce_dim', PCA()),
        ('classify', LinearSVC())
    ])
    N_FEATURES_OPTIONS = [2, 4, 8]
    C_OPTIONS = [1, 10, 100, 1000]
    param_grid = [
        {
            'reduce_dim': [PCA(iterated_power=7), NMF()],
            'reduce_dim__n_components': N_FEATURES_OPTIONS,
            'classify__C': C_OPTIONS
        },
        {
            'reduce_dim': [SelectKBest(chi2)],
            'reduce_dim__k': N_FEATURES_OPTIONS,
            'classify__C': C_OPTIONS
        },
    ]
    gs = GridSearchCV(pipeline, cv=3, n_jobs=1, param_grid=param_grid)
    expected_repr = """
GridSearchCV(cv=3, error_score='raise-deprecating',
    estimator=Pipeline(
        steps=[('reduce_dim',
                PCA(copy=True, iterated_power='auto', n_components=None,
                    random_state=None, svd_solver='auto', tol=0.0,
                    whiten=False))
               ('classify',
                LinearSVC(C=1.0, class_weight=None, dual=True,
                    fit_intercept=True, intercept_scaling=1,
                    loss='squared_hinge', max_iter=1000, multi_class='ovr',
                    penalty='l2', random_state=None, tol=0.0001, verbose=0))],
        memory=None),
    fit_params=None, iid='warn', n_jobs=1,
    param_grid=[
        {
            'reduce_dim': [
                PCA(copy=True, iterated_power=7, n_components=None,
                    random_state=None, svd_solver='auto', tol=0.0,
                    whiten=False),
                NMF(alpha=0.0, beta_loss='frobenius', init=None, l1_ratio=0.0,
                    max_iter=200, n_components=None, random_state=None,
                    shuffle=False, solver='cd', tol=0.0001, verbose=0)],
            'reduce_dim__n_components': [2, 4, 8],
            'classify__C': [1, 10, 100, 1000]},
        {'reduce_dim': [SelectKBest(k=10, score_func=chi2)],
            'reduce_dim__k': [2, 4, 8], 'classify__C': [1, 10, 100, 1000]}],
    pre_dispatch='2*n_jobs', refit=True, return_train_score='warn',
    scoring=None, verbose=0)"""

    expected_repr = expected_repr[1:]  # Remove first \n
    repr_ = f(gs)
    try:
        assert expected_repr == repr_
    except AssertionError:  # Older versions, different param orders
        assert set(expected_repr) == set(repr_)
    #  assert max(len(line) for line in repr_.split('\n')) <= f.width
