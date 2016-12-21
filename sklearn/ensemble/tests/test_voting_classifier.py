"""Testing for the VotingClassifier"""

import numpy as np
from sklearn.utils.testing import assert_almost_equal, assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn import datasets
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_multilabel_classification
from sklearn.svm import SVC
from sklearn.multiclass import OneVsRestClassifier
from sklearn.neighbors import KNeighborsClassifier


# Load the iris dataset and randomly permute it
iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target


def test_estimator_init():
    eclf = VotingClassifier(estimators=[])
    msg = ('Invalid `estimators` attribute, `estimators` should be'
           ' a list of (string, estimator) tuples')
    assert_raise_message(AttributeError, msg, eclf.fit, X, y)

    clf = LogisticRegression(random_state=1)

    eclf = VotingClassifier(estimators=[('lr', clf)], voting='error')
    msg = ('Voting must be \'soft\' or \'hard\'; got (voting=\'error\')')
    assert_raise_message(ValueError, msg, eclf.fit, X, y)

    eclf = VotingClassifier(estimators=[('lr', clf)], weights=[1, 2])
    msg = ('Number of classifiers and weights must be equal'
           '; got 2 weights, 1 estimators')
    assert_raise_message(ValueError, msg, eclf.fit, X, y)


def test_predictproba_hardvoting():
    eclf = VotingClassifier(estimators=[('lr1', LogisticRegression()),
                                        ('lr2', LogisticRegression())],
                            voting='hard')
    msg = "predict_proba is not available when voting='hard'"
    assert_raise_message(AttributeError, msg, eclf.predict_proba, X)


def test_notfitted():
    eclf = VotingClassifier(estimators=[('lr1', LogisticRegression()),
                                        ('lr2', LogisticRegression())],
                            voting='soft')
    msg = ("This VotingClassifier instance is not fitted yet. Call \'fit\'"
           " with appropriate arguments before using this method.")
    assert_raise_message(NotFittedError, msg, eclf.predict_proba, X)


def test_majority_label_iris():
    """Check classification by majority label on dataset iris."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    eclf = VotingClassifier(estimators=[
                ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                voting='hard')
    scores = cross_val_score(eclf, X, y, cv=5, scoring='accuracy')
    assert_almost_equal(scores.mean(), 0.95, decimal=2)


def test_tie_situation():
    """Check voting classifier selects smaller class label in tie situation."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    eclf = VotingClassifier(estimators=[('lr', clf1), ('rf', clf2)],
                            voting='hard')
    assert_equal(clf1.fit(X, y).predict(X)[73], 2)
    assert_equal(clf2.fit(X, y).predict(X)[73], 1)
    assert_equal(eclf.fit(X, y).predict(X)[73], 1)


def test_weights_iris():
    """Check classification by average probabilities on dataset iris."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    eclf = VotingClassifier(estimators=[
                            ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                            voting='soft',
                            weights=[1, 2, 10])
    scores = cross_val_score(eclf, X, y, cv=5, scoring='accuracy')
    assert_almost_equal(scores.mean(), 0.93, decimal=2)


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

    eclf = VotingClassifier(estimators=[
                            ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                            voting='hard',
                            weights=[1, 1, 1])
    assert_equal(all(eclf.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))

    eclf = VotingClassifier(estimators=[
                            ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                            voting='soft',
                            weights=[1, 1, 1])
    assert_equal(all(eclf.fit(X, y).predict(X)), all([1, 1, 1, 2, 2, 2]))


def test_predict_proba_on_toy_problem():
    """Calculate predicted probabilities on toy dataset."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    X = np.array([[-1.1, -1.5], [-1.2, -1.4], [-3.4, -2.2], [1.1, 1.2]])
    y = np.array([1, 1, 2, 2])

    clf1_res = np.array([[0.59790391, 0.40209609],
                         [0.57622162, 0.42377838],
                         [0.50728456, 0.49271544],
                         [0.40241774, 0.59758226]])

    clf2_res = np.array([[0.8, 0.2],
                         [0.8, 0.2],
                         [0.2, 0.8],
                         [0.3, 0.7]])

    clf3_res = np.array([[0.9985082, 0.0014918],
                         [0.99845843, 0.00154157],
                         [0., 1.],
                         [0., 1.]])

    t00 = (2*clf1_res[0][0] + clf2_res[0][0] + clf3_res[0][0]) / 4
    t11 = (2*clf1_res[1][1] + clf2_res[1][1] + clf3_res[1][1]) / 4
    t21 = (2*clf1_res[2][1] + clf2_res[2][1] + clf3_res[2][1]) / 4
    t31 = (2*clf1_res[3][1] + clf2_res[3][1] + clf3_res[3][1]) / 4

    eclf = VotingClassifier(estimators=[
                            ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                            voting='soft',
                            weights=[2, 1, 1])
    eclf_res = eclf.fit(X, y).predict_proba(X)

    assert_almost_equal(t00, eclf_res[0][0], decimal=1)
    assert_almost_equal(t11, eclf_res[1][1], decimal=1)
    assert_almost_equal(t21, eclf_res[2][1], decimal=1)
    assert_almost_equal(t31, eclf_res[3][1], decimal=1)

    try:
        eclf = VotingClassifier(estimators=[
                                ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                                voting='hard')
        eclf.fit(X, y).predict_proba(X)

    except AttributeError:
        pass
    else:
        raise AssertionError('AttributeError for voting == "hard"'
                             ' and with predict_proba not raised')


def test_multilabel():
    """Check if error is raised for multilabel classification."""
    X, y = make_multilabel_classification(n_classes=2, n_labels=1,
                                          allow_unlabeled=False,
                                          random_state=123)
    clf = OneVsRestClassifier(SVC(kernel='linear'))

    eclf = VotingClassifier(estimators=[('ovr', clf)], voting='hard')

    try:
        eclf.fit(X, y)
    except NotImplementedError:
        return


def test_gridsearch():
    """Check GridSearch support."""
    clf1 = LogisticRegression(random_state=1)
    clf2 = RandomForestClassifier(random_state=1)
    clf3 = GaussianNB()
    eclf = VotingClassifier(estimators=[
                ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
                voting='soft')

    params = {'lr__C': [1.0, 100.0],
              'voting': ['soft', 'hard'],
              'weights': [[0.5, 0.5, 0.5], [1.0, 0.5, 0.5]]}

    grid = GridSearchCV(estimator=eclf, param_grid=params, cv=5)
    grid.fit(iris.data, iris.target)


def test_parallel_predict():
    """Check parallel backend of VotingClassifier on toy dataset."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    X = np.array([[-1.1, -1.5], [-1.2, -1.4], [-3.4, -2.2], [1.1, 1.2]])
    y = np.array([1, 1, 2, 2])

    eclf1 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
        voting='soft',
        n_jobs=1).fit(X, y)
    eclf2 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
        voting='soft',
        n_jobs=2).fit(X, y)

    assert_array_equal(eclf1.predict(X), eclf2.predict(X))
    assert_array_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))


def test_sample_weight():
    """Tests sample_weight parameter of VotingClassifier"""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = SVC(probability=True, random_state=123)
    eclf1 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('svc', clf3)],
        voting='soft').fit(X, y, sample_weight=np.ones((len(y),)))
    eclf2 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('svc', clf3)],
        voting='soft').fit(X, y)
    assert_array_equal(eclf1.predict(X), eclf2.predict(X))
    assert_array_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))

    sample_weight = np.random.RandomState(123).uniform(size=(len(y),))
    eclf3 = VotingClassifier(estimators=[('lr', clf1)], voting='soft')
    eclf3.fit(X, y, sample_weight)
    clf1.fit(X, y, sample_weight)
    assert_array_equal(eclf3.predict(X), clf1.predict(X))
    assert_array_equal(eclf3.predict_proba(X), clf1.predict_proba(X))

    clf4 = KNeighborsClassifier()
    eclf3 = VotingClassifier(estimators=[
        ('lr', clf1), ('svc', clf3), ('knn', clf4)],
        voting='soft')
    msg = ('Underlying estimator \'knn\' does not support sample weights.')
    assert_raise_message(ValueError, msg, eclf3.fit, X, y, sample_weight)


def test_estimator_weights_format():
    # Test estimator weights inputs as list and array
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    eclf1 = VotingClassifier(estimators=[
                ('lr', clf1), ('rf', clf2)],
                weights=[1, 2],
                voting='soft')
    eclf2 = VotingClassifier(estimators=[
                ('lr', clf1), ('rf', clf2)],
                weights=np.array((1, 2)),
                voting='soft')
    eclf1.fit(X, y)
    eclf2.fit(X, y)
    assert_array_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))
