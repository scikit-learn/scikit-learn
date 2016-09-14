"""Testing for the StackingClassifier"""

from sklearn.utils.testing import assert_raise_message
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn import datasets
from sklearn.datasets import make_multilabel_classification
from sklearn.svm import SVC
from sklearn.multiclass import OneVsRestClassifier


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


def test_gridsearch():
    """Check GridSearch support."""
    clf1 = LogisticRegression(random_state=1)
    clf2 = RandomForestClassifier(random_state=1)
    clf3 = GaussianNB()
    clf_meta = LogisticRegression(random_state=1)
    eclf = StackingClassifier(estimators=[
                ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                meta_estimator=clf_meta)

    params = {'lr__C': [1.0, 100.0],
              'method': ['predict', 'predict_proba', 'predict_log_proba'],
              'meta_estimator__C': [10.0, 20.0]}

    grid = GridSearchCV(estimator=eclf, param_grid=params, cv=5)
    grid.fit(iris.data, iris.target)


def test_parallel_predict():
    """Check parallel backend of StackingClassifier on toy dataset."""
    pass


def test_sample_weight():
    """Tests sample_weight parameter of StackingClassifier"""
    pass


def test_predict_on_toy_problem():
    """Manually check predicted class labels for toy dataset."""
    pass


def test_predict_proba_on_toy_problem():
    """Calculate predicted probabilities on toy dataset."""
    pass
