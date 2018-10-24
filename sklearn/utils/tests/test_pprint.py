import re

from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn import set_config


# Ignore flake8 (lots of line too long issues)
# flake8: noqa

def test_basic():
    # Basic pprint test
    lr = LogisticRegression()
    expected = """
LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
                   intercept_scaling=1, max_iter=100, multi_class='warn',
                   n_jobs=None, penalty='l2', random_state=None, solver='warn',
                   tol=0.0001, verbose=0, warm_start=False)"""

    expected = expected[1:]  # remove first \n
    assert lr.__repr__() == expected


def test_changed_only():
    # Make sure the changed_only param is correctly used
    set_config(print_changed_only=True)
    lr = LogisticRegression(C=99)
    expected = """LogisticRegression(C=99)"""
    assert lr.__repr__() == expected
    set_config(print_changed_only=False)


def test_pipeline():
    # Render a pipeline object
    pipeline = make_pipeline(StandardScaler(), LogisticRegression(C=999))
    expected = """
Pipeline(memory=None,
         steps=[('standardscaler',
                 StandardScaler(copy=True, with_mean=True, with_std=True)),
                ('logisticregression',
                 LogisticRegression(C=999, class_weight=None, dual=False,
                                    fit_intercept=True, intercept_scaling=1,
                                    max_iter=100, multi_class='warn',
                                    n_jobs=None, penalty='l2',
                                    random_state=None, solver='warn',
                                    tol=0.0001, verbose=0, warm_start=False))])"""
    expected = expected[1:]  # remove first \n
    assert pipeline.__repr__() == expected


def test_deeply_nested():
    # Render a deeply nested estimator
    rfe = RFE(RFE(RFE(RFE(RFE(RFE(RFE(LogisticRegression())))))))
    expected = """
RFE(estimator=RFE(estimator=RFE(estimator=RFE(estimator=RFE(estimator=RFE(estimator=RFE(estimator=LogisticRegression(C=1.0,
                                                                                                                     class_weight=None,
                                                                                                                     dual=False,
                                                                                                                     fit_intercept=True,
                                                                                                                     intercept_scaling=1,
                                                                                                                     max_iter=100,
                                                                                                                     multi_class='warn',
                                                                                                                     n_jobs=None,
                                                                                                                     penalty='l2',
                                                                                                                     random_state=None,
                                                                                                                     solver='warn',
                                                                                                                     tol=0.0001,
                                                                                                                     verbose=0,
                                                                                                                     warm_start=False),
                                                                                        n_features_to_select=None,
                                                                                        step=1,
                                                                                        verbose=0),
                                                                          n_features_to_select=None,
                                                                          step=1,
                                                                          verbose=0),
                                                            n_features_to_select=None,
                                                            step=1, verbose=0),
                                              n_features_to_select=None, step=1,
                                              verbose=0),
                                n_features_to_select=None, step=1, verbose=0),
                  n_features_to_select=None, step=1, verbose=0),
    n_features_to_select=None, step=1, verbose=0)"""

    expected = expected[1:]  # remove first \n
    assert rfe.__repr__() == expected


def test_gridsearch():
    # render a gridsearch
    param_grid = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                   'C': [1, 10, 100, 1000]},
                  {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
    gs = GridSearchCV(SVC(), param_grid, cv=5)

    expected = """
GridSearchCV(cv=5, error_score='raise-deprecating',
             estimator=SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
                           decision_function_shape='ovr', degree=3,
                           gamma='auto_deprecated', kernel='rbf', max_iter=-1,
                           probability=False, random_state=None, shrinking=True,
                           tol=0.001, verbose=False),
             iid='warn', n_jobs=None,
             param_grid=[{'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001],
                          'kernel': ['rbf']},
                         {'C': [1, 10, 100, 1000], 'kernel': ['linear']}],
             pre_dispatch='2*n_jobs', refit=True, return_train_score=False,
             scoring=None, verbose=0)"""

    expected = expected[1:]  # remove first \n
    assert gs.__repr__() == expected


def test_gridsearch_pipeline():
    # render a pipeline inside a gridsearch
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
        }
    ]
    gspipline = GridSearchCV(pipeline, cv=3, n_jobs=1, param_grid=param_grid)
    expected = """
GridSearchCV(cv=3, error_score='raise-deprecating',
             estimator=Pipeline(memory=None,
                                steps=[('reduce_dim',
                                        PCA(copy=True, iterated_power='auto',
                                            n_components=None,
                                            random_state=None,
                                            svd_solver='auto', tol=0.0,
                                            whiten=False)),
                                       ('classify',
                                        LinearSVC(C=1.0, class_weight=None,
                                                  dual=True, fit_intercept=True,
                                                  intercept_scaling=1,
                                                  loss='squared_hinge',
                                                  max_iter=1000,
                                                  multi_class='ovr',
                                                  penalty='l2',
                                                  random_state=None, tol=0.0001,
                                                  verbose=0))]),
             iid='warn', n_jobs=1,
             param_grid=[{'classify__C': [1, 10, 100, 1000],
                          'reduce_dim': [PCA(copy=True, iterated_power=7,
                                             n_components=None,
                                             random_state=None,
                                             svd_solver='auto', tol=0.0,
                                             whiten=False),
                                         NMF(alpha=0.0, beta_loss='frobenius',
                                             init=None, l1_ratio=0.0,
                                             max_iter=200, n_components=None,
                                             random_state=None, shuffle=False,
                                             solver='cd', tol=0.0001,
                                             verbose=0)],
                          'reduce_dim__n_components': [2, 4, 8]},
                         {'classify__C': [1, 10, 100, 1000],
                          'reduce_dim': [SelectKBest(k=10,
                                                     score_func=<function chi2 at some_address>)],
                          'reduce_dim__k': [2, 4, 8]}],
             pre_dispatch='2*n_jobs', refit=True, return_train_score=False,
             scoring=None, verbose=0)"""

    expected = expected[1:]  # remove first \n
    repr_ = gspipline.__repr__()
    # Remove address of '<function chi2 at 0x.....>' for reproducibility
    repr_ = re.sub('function chi2 at 0x.*>',
                   'function chi2 at some_address>', repr_)
    assert repr_ == expected
