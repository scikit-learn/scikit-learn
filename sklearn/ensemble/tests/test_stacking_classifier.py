"""Testing for the StackingClassifier"""
from sklearn.utils.testing import assert_almost_equal, assert_array_equal
from sklearn.utils.testing import assert_equal, assert_array_almost_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn import datasets
from sklearn.datasets import make_multilabel_classification
from sklearn.svm import SVC
from sklearn.multiclass import OneVsRestClassifier
import numpy as np
from sklearn.neighbors import KNeighborsClassifier


# Load the iris dataset and randomly permute it
iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target


def test_estimator_init():
    clf = LogisticRegression(random_state=1)
    eclf = StackingClassifier(estimators=[], meta_estimator=clf)
    msg = ('Invalid `estimators` attribute, `estimators` should be'
           ' a list of (string, estimator) tuples')
    assert_raise_message(AttributeError, msg, eclf.fit, X, y)

    eclf = StackingClassifier(estimators=[('lr', clf)], meta_estimator=None)
    msg = ('Invalid `meta_estimator` attribute, '
           '`meta_estimator` should be a classifier')
    assert_raise_message(AttributeError, msg, eclf.fit, X, y)

    clf_no_proba = SVC()
    eclf = StackingClassifier(estimators=[('SVC', clf_no_proba)],
                              meta_estimator=clf, method='predict_proba')
    msg = 'Underlying estimator SVC does not support predict_proba'
    assert_raise_message(ValueError, msg, eclf.fit, X, y)


def test_notfitted():
    eclf = StackingClassifier(estimators=[('lr1', LogisticRegression()),
                                          ('lr2', LogisticRegression())],
                              meta_estimator=LogisticRegression())
    msg = ("This StackingClassifier instance is not fitted yet. Call \'fit\'"
           " with appropriate arguments before using this method.")
    assert_raise_message(NotFittedError, msg, eclf.predict, X)


def test_multilabel():
    """Check if error is raised for multilabel classification."""
    X, y = make_multilabel_classification(n_classes=2, n_labels=1,
                                          allow_unlabeled=False,
                                          random_state=123)
    clf = OneVsRestClassifier(SVC(kernel='linear'))
    clf_meta = OneVsRestClassifier(SVC(kernel='linear'))

    eclf = StackingClassifier(estimators=[('ovr1', clf)],
                              meta_estimator=clf_meta,
                              method='decision_function')

    try:
        eclf.fit(X, y)
    except NotImplementedError:
        return


def test_cv():
    """Test cross-validation option cv=1"""
    clf1 = LogisticRegression(random_state=1)
    clfm = GaussianNB()
    eclf = StackingClassifier(estimators=[('lr', clf1)], cv=1,
                              meta_estimator=clfm).fit(X, y)
    s = clf1.fit(X, y).predict_proba(X)[:, 1:]
    assert_array_equal(clfm.fit(s, y).predict(s), eclf.predict(X))


def test_gridsearch():
    """Check GridSearch support."""
    clf1 = LogisticRegression(random_state=1)
    clf2 = RandomForestClassifier(random_state=1)
    clf3 = GaussianNB()
    clf_meta = LogisticRegression(random_state=1)
    eclf = StackingClassifier(estimators=[
                ('lr', clf1), ('rf', clf2), ('nb', clf3)],
                meta_estimator=clf_meta)

    params = {'lr__C': [1.0, 100.0],
              'method': ['predict', 'predict_proba', 'predict_log_proba'],
              'meta_estimator__C': [10.0, 20.0]}

    grid = GridSearchCV(estimator=eclf, param_grid=params, cv=5)

    grid.fit(iris.data, iris.target)


def test_parallel_predict():
    """Check parallel backend of StackingClassifier on iris."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    clfm = GaussianNB()

    eclf1 = StackingClassifier(
        estimators=[('lr', clf1), ('rf', clf2), ('nb', clf3)],
        meta_estimator=clfm,
        n_jobs=1).fit(X, y)
    eclf2 = StackingClassifier(
        estimators=[('lr', clf1), ('rf', clf2), ('nb', clf3)],
        meta_estimator=clfm,
        n_jobs=2).fit(X, y)

    assert_array_equal(eclf1.predict(X), eclf2.predict(X))
    assert_array_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))


def test_sample_weight():
    """Tests sample_weight parameter of StackingClassifier"""
    decimal = 10

    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    clfm = GaussianNB()
    eclf1 = StackingClassifier(
        estimators=[('lr', clf1), ('rf', clf2), ('nb', clf3)],
        meta_estimator=clfm
    ).fit(X, y, sample_weight=np.ones((len(y),)))
    eclf2 = StackingClassifier(
        estimators=[('lr', clf1), ('rf', clf2), ('nb', clf3)],
        meta_estimator=clfm
    ).fit(X, y)
    assert_array_almost_equal(eclf1.predict(X), eclf2.predict(X),
                              decimal=decimal)
    assert_array_almost_equal(eclf1.predict_proba(X), eclf2.predict_proba(X),
                              decimal=decimal)

    sample_weight = np.random.RandomState(123).uniform(size=(len(y),))
    clf4 = KNeighborsClassifier()
    eclf3 = StackingClassifier(
        estimators=[('lr', clf1), ('knn', clf4)],
        meta_estimator=clfm)
    msg = ('Underlying estimator `knn` does not support `sample_weight`.')
    assert_raise_message(ValueError, msg, eclf3.fit, X, y,
                         sample_weight=sample_weight)
    eclf4 = StackingClassifier(
        estimators=[('lr', clf1), ('nb', clf3)],
        meta_estimator=clf4)
    msg = ("Underlying meta estimator "
           "<class 'sklearn.neighbors.classification.KNeighborsClassifier'> "
           "does not support `sample_weight`.")
    assert_raise_message(ValueError, msg, eclf4.fit, X, y,
                         sample_weight=sample_weight)


def test_classify_iris():
    """Check classification by majority label on dataset iris."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    clfm = GaussianNB()
    eclf = StackingClassifier(estimators=[('lg', clf1),
                                          ('rf', clf2),
                                          ('nb', clf3)],
                              meta_estimator=clfm)
    scores = cross_val_score(eclf, X, y, cv=5, scoring='accuracy')
    assert_almost_equal(scores.mean(), 0.95, decimal=2)


def test_predict_on_toy_problem():
    """Manually check predicted class labels for toy dataset."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()

    X = np.array([[-1.1, -1.5],
                  [-1.2, -1.4],
                  [-3.4, -2.2],
                  [1.1, 1.2],
                  [2.1, 1.4],
                  [3.1, 2.3]])

    y = np.array([1, 1, 1, 2, 2, 2])

    assert_equal(all(clf1.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))
    assert_equal(all(clf2.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))
    assert_equal(all(clf3.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))

    eclf = StackingClassifier(estimators=[('rf', clf2), ('nb', clf3)],
                              meta_estimator=clf1)
    assert_equal(all(eclf.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))

    eclf = StackingClassifier(estimators=[('lr', clf1), ('nb', clf3)],
                              meta_estimator=clf2)
    assert_equal(all(eclf.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))

    eclf = StackingClassifier(estimators=[('lr', clf1), ('rf', clf2)],
                              meta_estimator=clf3)
    assert_equal(all(eclf.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))
