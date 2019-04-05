"""
Testing Recursive feature elimination
"""

import pytest
import numpy as np
from scipy import sparse

from sklearn.feature_selection.rfe import RFE, RFECV
from sklearn.datasets import load_iris, make_friedman1
from sklearn.svm import SVC, SVR
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import get_scorer, make_scorer, zero_one_loss
from sklearn.model_selection import cross_val_score, GroupKFold
from sklearn.utils import check_random_state
from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_array_almost_equal, assert_greater,
                                   assert_raises, ignore_warnings)


class MockClassifier:
    """
    Dummy classifier to test recursive feature elimination
    """

    def __init__(self, foo_param=0):
        self.foo_param = foo_param

    def fit(self, X, Y):
        assert len(X) == len(Y)
        self.coef_ = np.ones(X.shape[1], dtype=np.float64)
        return self

    def predict(self, T):
        return T.shape[0]

    predict_proba = predict
    decision_function = predict
    transform = predict

    def score(self, X=None, Y=None):
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score

    def get_params(self, deep=True):
        return {'foo_param': self.foo_param}

    def set_params(self, **params):
        return self


def test_rfe_features_importance():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    clf = RandomForestClassifier(n_estimators=20,
                                 random_state=generator, max_depth=2)
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    rfe.fit(X, y)
    assert_equal(len(rfe.ranking_), X.shape[1])

    clf_svc = SVC(kernel="linear")
    rfe_svc = RFE(estimator=clf_svc, n_features_to_select=4, step=0.1)
    rfe_svc.fit(X, y)

    # Check if the supports are equal
    assert_array_equal(rfe.get_support(), rfe_svc.get_support())


def test_rfe():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    X_sparse = sparse.csr_matrix(X)
    y = iris.target

    # dense model
    clf = SVC(kernel="linear")
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    rfe.fit(X, y)
    X_r = rfe.transform(X)
    clf.fit(X_r, y)
    assert_equal(len(rfe.ranking_), X.shape[1])

    # sparse model
    clf_sparse = SVC(kernel="linear")
    rfe_sparse = RFE(estimator=clf_sparse, n_features_to_select=4, step=0.1)
    rfe_sparse.fit(X_sparse, y)
    X_r_sparse = rfe_sparse.transform(X_sparse)

    assert_equal(X_r.shape, iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])

    assert_array_almost_equal(rfe.predict(X), clf.predict(iris.data))
    assert_equal(rfe.score(X, y), clf.score(iris.data, iris.target))
    assert_array_almost_equal(X_r, X_r_sparse.toarray())


def test_rfe_mockclassifier():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    # dense model
    clf = MockClassifier()
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    rfe.fit(X, y)
    X_r = rfe.transform(X)
    clf.fit(X_r, y)
    assert_equal(len(rfe.ranking_), X.shape[1])
    assert_equal(X_r.shape, iris.data.shape)


def test_rfecv():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = list(iris.target)   # regression test: list should be supported

    # Test using the score function
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=5)
    rfecv.fit(X, y)
    # non-regression test for missing worst feature:
    assert_equal(len(rfecv.grid_scores_), X.shape[1])
    assert_equal(len(rfecv.ranking_), X.shape[1])
    X_r = rfecv.transform(X)

    # All the noisy variable were filtered out
    assert_array_equal(X_r, iris.data)

    # same in sparse
    rfecv_sparse = RFECV(estimator=SVC(kernel="linear"), step=1, cv=5)
    X_sparse = sparse.csr_matrix(X)
    rfecv_sparse.fit(X_sparse, y)
    X_r_sparse = rfecv_sparse.transform(X_sparse)
    assert_array_equal(X_r_sparse.toarray(), iris.data)

    # Test using a customized loss function
    scoring = make_scorer(zero_one_loss, greater_is_better=False)
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=5,
                  scoring=scoring)
    ignore_warnings(rfecv.fit)(X, y)
    X_r = rfecv.transform(X)
    assert_array_equal(X_r, iris.data)

    # Test using a scorer
    scorer = get_scorer('accuracy')
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=5,
                  scoring=scorer)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)
    assert_array_equal(X_r, iris.data)

    # Test fix on grid_scores
    def test_scorer(estimator, X, y):
        return 1.0
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=5,
                  scoring=test_scorer)
    rfecv.fit(X, y)
    assert_array_equal(rfecv.grid_scores_, np.ones(len(rfecv.grid_scores_)))
    # In the event of cross validation score ties, the expected behavior of
    # RFECV is to return the FEWEST features that maximize the CV score.
    # Because test_scorer always returns 1.0 in this example, RFECV should
    # reduce the dimensionality to a single feature (i.e. n_features_ = 1)
    assert_equal(rfecv.n_features_, 1)

    # Same as the first two tests, but with step=2
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=2, cv=5)
    rfecv.fit(X, y)
    assert_equal(len(rfecv.grid_scores_), 6)
    assert_equal(len(rfecv.ranking_), X.shape[1])
    X_r = rfecv.transform(X)
    assert_array_equal(X_r, iris.data)

    rfecv_sparse = RFECV(estimator=SVC(kernel="linear"), step=2, cv=5)
    X_sparse = sparse.csr_matrix(X)
    rfecv_sparse.fit(X_sparse, y)
    X_r_sparse = rfecv_sparse.transform(X_sparse)
    assert_array_equal(X_r_sparse.toarray(), iris.data)

    # Verifying that steps < 1 don't blow up.
    rfecv_sparse = RFECV(estimator=SVC(kernel="linear"), step=.2, cv=5)
    X_sparse = sparse.csr_matrix(X)
    rfecv_sparse.fit(X_sparse, y)
    X_r_sparse = rfecv_sparse.transform(X_sparse)
    assert_array_equal(X_r_sparse.toarray(), iris.data)


def test_rfecv_mockclassifier():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = list(iris.target)   # regression test: list should be supported

    # Test using the score function
    rfecv = RFECV(estimator=MockClassifier(), step=1, cv=5)
    rfecv.fit(X, y)
    # non-regression test for missing worst feature:
    assert_equal(len(rfecv.grid_scores_), X.shape[1])
    assert_equal(len(rfecv.ranking_), X.shape[1])


def test_rfecv_verbose_output():
    # Check verbose=1 is producing an output.
    from io import StringIO
    import sys
    sys.stdout = StringIO()

    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = list(iris.target)

    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=5, verbose=1)
    rfecv.fit(X, y)

    verbose_output = sys.stdout
    verbose_output.seek(0)
    assert_greater(len(verbose_output.readline()), 0)


def test_rfecv_grid_scores_size():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = list(iris.target)   # regression test: list should be supported

    # Non-regression test for varying combinations of step and
    # min_features_to_select.
    for step, min_features_to_select in [[2, 1], [2, 2], [3, 3]]:
        rfecv = RFECV(estimator=MockClassifier(), step=step,
                      min_features_to_select=min_features_to_select, cv=5)
        rfecv.fit(X, y)

        score_len = np.ceil(
            (X.shape[1] - min_features_to_select) / step) + 1
        assert len(rfecv.grid_scores_) == score_len
        assert len(rfecv.ranking_) == X.shape[1]
        assert rfecv.n_features_ >= min_features_to_select


@pytest.mark.filterwarnings('ignore: The default value of cv')  # 0.22
def test_rfe_estimator_tags():
    rfe = RFE(SVC(kernel='linear'))
    assert_equal(rfe._estimator_type, "classifier")
    # make sure that cross-validation is stratified
    iris = load_iris()
    score = cross_val_score(rfe, iris.data, iris.target)
    assert_greater(score.min(), .7)


def test_rfe_min_step():
    n_features = 10
    X, y = make_friedman1(n_samples=50, n_features=n_features, random_state=0)
    n_samples, n_features = X.shape
    estimator = SVR(kernel="linear")

    # Test when floor(step * n_features) <= 0
    selector = RFE(estimator, step=0.01)
    sel = selector.fit(X, y)
    assert_equal(sel.support_.sum(), n_features // 2)

    # Test when step is between (0,1) and floor(step * n_features) > 0
    selector = RFE(estimator, step=0.20)
    sel = selector.fit(X, y)
    assert_equal(sel.support_.sum(), n_features // 2)

    # Test when step is an integer
    selector = RFE(estimator, step=5)
    sel = selector.fit(X, y)
    assert_equal(sel.support_.sum(), n_features // 2)


def test_number_of_subsets_of_features():
    # In RFE, 'number_of_subsets_of_features'
    # = the number of iterations in '_fit'
    # = max(ranking_)
    # = 1 + (n_features + step - n_features_to_select - 1) // step
    # After optimization #4534, this number
    # = 1 + np.ceil((n_features - n_features_to_select) / float(step))
    # This test case is to test their equivalence, refer to #4534 and #3824

    def formula1(n_features, n_features_to_select, step):
        return 1 + ((n_features + step - n_features_to_select - 1) // step)

    def formula2(n_features, n_features_to_select, step):
        return 1 + np.ceil((n_features - n_features_to_select) / float(step))

    # RFE
    # Case 1, n_features - n_features_to_select is divisible by step
    # Case 2, n_features - n_features_to_select is not divisible by step
    n_features_list = [11, 11]
    n_features_to_select_list = [3, 3]
    step_list = [2, 3]
    for n_features, n_features_to_select, step in zip(
            n_features_list, n_features_to_select_list, step_list):
        generator = check_random_state(43)
        X = generator.normal(size=(100, n_features))
        y = generator.rand(100).round()
        rfe = RFE(estimator=SVC(kernel="linear"),
                  n_features_to_select=n_features_to_select, step=step)
        rfe.fit(X, y)
        # this number also equals to the maximum of ranking_
        assert_equal(np.max(rfe.ranking_),
                     formula1(n_features, n_features_to_select, step))
        assert_equal(np.max(rfe.ranking_),
                     formula2(n_features, n_features_to_select, step))

    # In RFECV, 'fit' calls 'RFE._fit'
    # 'number_of_subsets_of_features' of RFE
    # = the size of 'grid_scores' of RFECV
    # = the number of iterations of the for loop before optimization #4534

    # RFECV, n_features_to_select = 1
    # Case 1, n_features - 1 is divisible by step
    # Case 2, n_features - 1 is not divisible by step
    n_features_to_select = 1
    n_features_list = [11, 10]
    step_list = [2, 2]
    for n_features, step in zip(n_features_list, step_list):
        generator = check_random_state(43)
        X = generator.normal(size=(100, n_features))
        y = generator.rand(100).round()
        rfecv = RFECV(estimator=SVC(kernel="linear"), step=step, cv=5)
        rfecv.fit(X, y)
        assert_equal(rfecv.grid_scores_.shape[0],
                     formula1(n_features, n_features_to_select, step))
        assert_equal(rfecv.grid_scores_.shape[0],
                     formula2(n_features, n_features_to_select, step))

    # Advanced step options
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target
    rfe = RFE(estimator=SVC(kernel="linear"), n_features_to_select=10,
              tune_step_at=5)
    assert_raises(ValueError, rfe.fit, X, y)

    n_features_list = [300] * 14
    n_features_to_select_list = [50] * 14
    step_list = [100, 100, 100, 100, 100, 100,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
    tune_step_at_list = [105, 105, 105, 0.325, 0.325, 0.325,
                         105, 105, 105, 105, 0.325, 0.325, 0.325, 0.325]
    tuning_step_list = [10, 0.125, 0.125, 10, 0.125, 0.125,
                        10, 10, 0.125, 0.125, 10, 10, 0.125, 0.125]
    reducing_step_list = [False, False, True, False, False, True,
                          False, True, False, True, False, True, False, True]
    n_remaining_feature_steps_list = [
        [300, 200, 105, 95, 85, 75, 65, 55, 50],
        [300, 200, 105, 92, 79, 66, 53, 50],
        [300, 200, 105, 92, 81, 71, 63, 56, 50],
        [300, 200, 100, 97, 87, 77, 67, 57, 50],
        [300, 200, 100, 97, 85, 73, 61, 50],
        [300, 200, 100, 97, 85, 75, 66, 58, 51, 50],
        [300, 240, 180, 120, 105, 95, 85, 75, 65, 55, 50],
        [300, 240, 192, 154, 124, 105, 95, 85, 75, 65, 55, 50],
        [300, 240, 180, 120, 105, 92, 79, 66, 53, 50],
        [300, 240, 192, 154, 124, 105, 92, 81, 71, 63, 56, 50],
        [300, 240, 180, 120, 97, 87, 77, 67, 57, 50],
        [300, 240, 192, 154, 124, 100, 97, 87, 77, 67, 57, 50],
        [300, 240, 180, 120, 97, 85, 73, 61, 50],
        [300, 240, 192, 154, 124, 100, 97, 85, 75, 66, 58, 51, 50],
    ]
    for (n_features, n_features_to_select, step,
         tune_step_at, tuning_step, reducing_step,
         n_remaining_feature_steps) in zip(
             n_features_list, n_features_to_select_list, step_list,
             tune_step_at_list, tuning_step_list, reducing_step_list,
             n_remaining_feature_steps_list):
        generator = check_random_state(43)
        X = generator.normal(size=(100, n_features))
        y = generator.rand(100).round()
        rfe = RFE(estimator=SVC(kernel="linear"),
                  n_features_to_select=n_features_to_select, step=step,
                  tune_step_at=tune_step_at, tuning_step=tuning_step,
                  reducing_step=reducing_step)
        rfe.fit(X, y)
        assert_array_equal(rfe.n_remaining_feature_steps_,
                           n_remaining_feature_steps)


@pytest.mark.filterwarnings('ignore: The default value of cv')  # 0.22
def test_rfe_cv_n_jobs():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    rfecv = RFECV(estimator=SVC(kernel='linear'))
    rfecv.fit(X, y)
    rfecv_ranking = rfecv.ranking_
    rfecv_grid_scores = rfecv.grid_scores_

    rfecv.set_params(n_jobs=2)
    rfecv.fit(X, y)
    assert_array_almost_equal(rfecv.ranking_, rfecv_ranking)
    assert_array_almost_equal(rfecv.grid_scores_, rfecv_grid_scores)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
def test_rfe_cv_groups():
    generator = check_random_state(0)
    iris = load_iris()
    number_groups = 4
    groups = np.floor(np.linspace(0, number_groups, len(iris.target)))
    X = iris.data
    y = (iris.target > 0).astype(int)

    est_groups = RFECV(
        estimator=RandomForestClassifier(random_state=generator),
        step=1,
        scoring='accuracy',
        cv=GroupKFold(n_splits=2)
    )
    est_groups.fit(X, y, groups=groups)
    assert est_groups.n_features_ > 0
