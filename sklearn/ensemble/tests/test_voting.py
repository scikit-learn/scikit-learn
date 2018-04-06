"""Testing for the VotingClassifier and AverageRegressor"""

import numpy as np
from sklearn.utils.testing import assert_almost_equal, assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal, assert_true, assert_false
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns_message
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
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin
from sklearn.ensemble import AverageRegressor
from sklearn.utils.validation import check_random_state
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest


# ===========================================================================
# tests for sklearn.ensemble.VotingClassifier
# ===========================================================================

# Load the iris dataset and randomly permute it
iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target


def test_estimator_init_voting_classifier():
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

    eclf = VotingClassifier(estimators=[('lr', clf), ('lr', clf)],
                            weights=[1, 2])
    msg = "Names provided are not unique: ['lr', 'lr']"
    assert_raise_message(ValueError, msg, eclf.fit, X, y)

    eclf = VotingClassifier(estimators=[('lr__', clf)])
    msg = "Estimator names must not contain __: got ['lr__']"
    assert_raise_message(ValueError, msg, eclf.fit, X, y)

    eclf = VotingClassifier(estimators=[('estimators', clf)])
    msg = "Estimator names conflict with constructor arguments: ['estimators']"
    assert_raise_message(ValueError, msg, eclf.fit, X, y)


def test_predictproba_hardvoting():
    eclf = VotingClassifier(estimators=[('lr1', LogisticRegression()),
                                        ('lr2', LogisticRegression())],
                            voting='hard')
    msg = "predict_proba is not available when voting='hard'"
    assert_raise_message(AttributeError, msg, eclf.predict_proba, X)


def test_notfitted_voting_classifier():
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


def test_gridsearch_voting_classifier():
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


def test_parallel_fit():
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
    assert_array_almost_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))


def test_sample_weight_voting_classifier():
    """Tests sample_weight parameter of VotingClassifier"""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = SVC(gamma='scale', probability=True, random_state=123)
    eclf1 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('svc', clf3)],
        voting='soft').fit(X, y, sample_weight=np.ones((len(y),)))
    eclf2 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('svc', clf3)],
        voting='soft').fit(X, y)
    assert_array_equal(eclf1.predict(X), eclf2.predict(X))
    assert_array_almost_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))

    sample_weight = np.random.RandomState(123).uniform(size=(len(y),))
    eclf3 = VotingClassifier(estimators=[('lr', clf1)], voting='soft')
    eclf3.fit(X, y, sample_weight)
    clf1.fit(X, y, sample_weight)
    assert_array_equal(eclf3.predict(X), clf1.predict(X))
    assert_array_almost_equal(eclf3.predict_proba(X), clf1.predict_proba(X))

    clf4 = KNeighborsClassifier()
    eclf3 = VotingClassifier(estimators=[
        ('lr', clf1), ('svc', clf3), ('knn', clf4)],
        voting='soft')
    msg = ('Underlying estimator \'knn\' does not support sample weights.')
    assert_raise_message(ValueError, msg, eclf3.fit, X, y, sample_weight)


def test_sample_weight_kwargs_voting_classifier():
    """Check that VotingClassifier passes sample_weight as kwargs"""
    class MockClassifier(BaseEstimator, ClassifierMixin):
        """Mock Classifier to check that sample_weight is received as kwargs"""
        def fit(self, X, y, *args, **sample_weight):
            assert_true('sample_weight' in sample_weight)

    clf = MockClassifier()
    eclf = VotingClassifier(estimators=[('mock', clf)], voting='soft')

    # Should not raise an error.
    eclf.fit(X, y, sample_weight=np.ones((len(y),)))


def test_set_params_voting_classifier():
    """set_params should be able to set estimators"""
    clf1 = LogisticRegression(random_state=123, C=1.0)
    clf2 = RandomForestClassifier(random_state=123, max_depth=None)
    clf3 = GaussianNB()
    eclf1 = VotingClassifier([('lr', clf1), ('rf', clf2)], voting='soft',
                             weights=[1, 2])
    assert_true('lr' in eclf1.named_estimators)
    assert_true(eclf1.named_estimators.lr is eclf1.estimators[0][1])
    assert_true(eclf1.named_estimators.lr is eclf1.named_estimators['lr'])
    eclf1.fit(X, y)
    assert_true('lr' in eclf1.named_estimators_)
    assert_true(eclf1.named_estimators_.lr is eclf1.estimators_[0])
    assert_true(eclf1.named_estimators_.lr is eclf1.named_estimators_['lr'])

    eclf2 = VotingClassifier([('lr', clf1), ('nb', clf3)], voting='soft',
                             weights=[1, 2])
    eclf2.set_params(nb=clf2).fit(X, y)
    assert_false(hasattr(eclf2, 'nb'))

    assert_array_equal(eclf1.predict(X), eclf2.predict(X))
    assert_array_almost_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))
    assert_equal(eclf2.estimators[0][1].get_params(), clf1.get_params())
    assert_equal(eclf2.estimators[1][1].get_params(), clf2.get_params())

    eclf1.set_params(lr__C=10.0)
    eclf2.set_params(nb__max_depth=5)

    assert_true(eclf1.estimators[0][1].get_params()['C'] == 10.0)
    assert_true(eclf2.estimators[1][1].get_params()['max_depth'] == 5)
    assert_equal(eclf1.get_params()["lr__C"],
                 eclf1.get_params()["lr"].get_params()['C'])


def test_set_estimator_none():
    """VotingClassifier set_params should be able to set estimators as None"""
    # Test predict
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    eclf1 = VotingClassifier(estimators=[('lr', clf1), ('rf', clf2),
                                         ('nb', clf3)],
                             voting='hard', weights=[1, 0, 0.5]).fit(X, y)

    eclf2 = VotingClassifier(estimators=[('lr', clf1), ('rf', clf2),
                                         ('nb', clf3)],
                             voting='hard', weights=[1, 1, 0.5])
    eclf2.set_params(rf=None).fit(X, y)
    assert_array_equal(eclf1.predict(X), eclf2.predict(X))

    assert_true(dict(eclf2.estimators)["rf"] is None)
    assert_true(len(eclf2.estimators_) == 2)
    assert_true(all([not isinstance(est, RandomForestClassifier) for est in
                     eclf2.estimators_]))
    assert_true(eclf2.get_params()["rf"] is None)

    eclf1.set_params(voting='soft').fit(X, y)
    eclf2.set_params(voting='soft').fit(X, y)
    assert_array_equal(eclf1.predict(X), eclf2.predict(X))
    assert_array_almost_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))
    msg = ('All estimators are None. At least one is required'
           ' to be a classifier!')
    assert_raise_message(
        ValueError, msg, eclf2.set_params(lr=None, rf=None, nb=None).fit, X, y)

    # Test soft voting transform
    X1 = np.array([[1], [2]])
    y1 = np.array([1, 2])
    eclf1 = VotingClassifier(estimators=[('rf', clf2), ('nb', clf3)],
                             voting='soft', weights=[0, 0.5]).fit(X1, y1)

    eclf2 = VotingClassifier(estimators=[('rf', clf2), ('nb', clf3)],
                             voting='soft', weights=[1, 0.5])
    eclf2.set_params(rf=None).fit(X1, y1)
    assert_array_almost_equal(eclf1.transform(X1),
                              np.array([[[0.7, 0.3], [0.3, 0.7]],
                                        [[1., 0.], [0., 1.]]]))
    assert_array_almost_equal(eclf2.transform(X1),
                              np.array([[[1., 0.],
                                         [0., 1.]]]))
    eclf1.set_params(voting='hard')
    eclf2.set_params(voting='hard')
    assert_array_equal(eclf1.transform(X1), np.array([[0, 0], [1, 1]]))
    assert_array_equal(eclf2.transform(X1), np.array([[0], [1]]))


def test_estimator_weights_format_voting_classifier():
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
    assert_array_almost_equal(eclf1.predict_proba(X), eclf2.predict_proba(X))


def test_transform():
    """Check transform method of VotingClassifier on toy dataset."""
    clf1 = LogisticRegression(random_state=123)
    clf2 = RandomForestClassifier(random_state=123)
    clf3 = GaussianNB()
    X = np.array([[-1.1, -1.5], [-1.2, -1.4], [-3.4, -2.2], [1.1, 1.2]])
    y = np.array([1, 1, 2, 2])

    eclf1 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
        voting='soft').fit(X, y)
    eclf2 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
        voting='soft',
        flatten_transform=True).fit(X, y)
    eclf3 = VotingClassifier(estimators=[
        ('lr', clf1), ('rf', clf2), ('gnb', clf3)],
        voting='soft',
        flatten_transform=False).fit(X, y)

    warn_msg = ("'flatten_transform' default value will be "
                "changed to True in 0.21. "
                "To silence this warning you may"
                " explicitly set flatten_transform=False.")
    res = assert_warns_message(DeprecationWarning, warn_msg,
                               eclf1.transform, X)
    assert_array_equal(res.shape, (3, 4, 2))
    assert_array_equal(eclf2.transform(X).shape, (4, 6))
    assert_array_equal(eclf3.transform(X).shape, (3, 4, 2))
    assert_array_almost_equal(res.swapaxes(0, 1).reshape((4, 6)),
                              eclf2.transform(X))
    assert_array_almost_equal(
            eclf3.transform(X).swapaxes(0, 1).reshape((4, 6)),
            eclf2.transform(X)
    )

# ===========================================================================
# tests for sklearn.ensemble.AverageRegressor
# ===========================================================================

rng = check_random_state(0)
# load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

X_train, X_test, y_train, y_test = train_test_split(boston.data[:150],
                                                    boston.target[:150],
                                                    random_state=rng)


def test_estimator_init_average_regressor():
    """Test various init configurations"""

    lr = LinearRegression()
    # test init with empty list of estimators
    ensemble = AverageRegressor(estimators=[])
    msg = ('Invalid `estimators` attribute, `estimators` should be'
           ' a list of (string, estimator) tuples')
    assert_raise_message(AttributeError, msg, ensemble.fit, X_train, y_train)

    # test init with all estimators set to None
    ensemble = AverageRegressor(estimators=[('r1', None),
                                            ('r2', None)])
    msg = ('All estimators are None. At least one is '
           'required to be a regressor!')
    assert_raise_message(ValueError, msg, ensemble.fit, X_train, y_train)

    # test init with shape of estimators different than that of weights
    ensemble = AverageRegressor(estimators=[('lr', lr)], weights=[1, 2])
    msg = ('Number of regressors and weights must be equal'
           '; got 2 weights, 1 estimators')
    assert_raise_message(ValueError, msg, ensemble.fit, X_train, y_train)

    # valide that estimators names are unique
    ensemble = AverageRegressor(estimators=[('lr', lr), ('lr', lr)],
                                weights=[1, 2])
    msg = "Names provided are not unique: ['lr', 'lr']"
    assert_raise_message(ValueError, msg, ensemble.fit, X_train, y_train)

    ensemble = AverageRegressor(estimators=[('lr__', lr)])
    msg = "Estimator names must not contain __: got ['lr__']"
    assert_raise_message(ValueError, msg, ensemble.fit, X_train, y_train)

    ensemble = AverageRegressor(estimators=[('estimators', lr)])
    msg = "Estimator names conflict with constructor arguments: ['estimators']"
    assert_raise_message(ValueError, msg, ensemble.fit, X_train, y_train)


def test_notfitted_average_regressor():
    """Test raise error when predicting with non-fitted estimator"""
    ensemble = AverageRegressor(estimators=[('lr1', LinearRegression()),
                                            ('lr2', LinearRegression())])
    msg = ("This AverageRegressor instance is not fitted yet. Call \'fit\'"
           " with appropriate arguments before using this method.")
    assert_raise_message(NotFittedError, msg, ensemble.predict, X_train)


def test_single_estimator():
    """ Check singleton ensembles."""

    clf1 = AverageRegressor(estimators=[
        ('knn', KNeighborsRegressor())]).fit(X_train, y_train)

    clf2 = KNeighborsRegressor().fit(X_train, y_train)

    assert_array_almost_equal(clf1.predict(X_test), clf2.predict(X_test))


def test_average_prediction():
    """ Check AverageRegressor with multiple models"""
    reg1 = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                            random_state=rng,
                            n_estimators=100)
    reg2 = RandomForestRegressor(random_state=rng,
                                 n_estimators=100)
    reg3 = GradientBoostingRegressor(random_state=rng,
                                     n_estimators=100)

    ensemble = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)])

    ensemble.fit(X_train, y_train)
    scores = cross_val_score(ensemble,
                             X_train,
                             y_train,
                             cv=5,
                             scoring='r2')

    assert_almost_equal(scores.mean(), 0.85, decimal=1)


def test_weights_average_regression_boston():
    """Check weighted average regression prediction on boston dataset."""

    reg1 = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                            random_state=rng,
                            n_estimators=100)
    reg2 = RandomForestRegressor(random_state=rng,
                                 n_estimators=100)
    reg3 = GradientBoostingRegressor(random_state=rng,
                                     n_estimators=100)

    ensemble = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)],
        weights=[1, 2, 10])

    scores = cross_val_score(ensemble,
                             X_train,
                             y_train,
                             cv=5,
                             scoring='r2')
    assert_almost_equal(scores.mean(), 0.87, decimal=1)


def test_parallel_regression():
    """ Check parallel regression. """

    estimators = [('dtr', DecisionTreeRegressor())]
    ensemble = AverageRegressor(estimators=estimators,
                                n_jobs=3).fit(X_train, y_train)

    ensemble.set_params(n_jobs=1)
    y1 = ensemble.predict(X_test)
    ensemble.set_params(n_jobs=2)
    y2 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = AverageRegressor(estimators=estimators,
                                n_jobs=1).fit(X_train, y_train)

    y3 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y3)


def test_gridsearch_average_regressor():
    """Check GridSearch support."""

    reg1 = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                            random_state=rng,
                            n_estimators=10)
    reg2 = RandomForestRegressor(random_state=rng,
                                 n_estimators=10)
    reg3 = GradientBoostingRegressor(random_state=rng,
                                     n_estimators=10)

    ensemble = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)]
                    )

    params = {'weights': [[0.5, 0.5, 0.5], [1.0, 0.5, 0.5]]}

    grid = GridSearchCV(estimator=ensemble, param_grid=params, cv=5)
    grid.fit(X_train, y_train)


def test_estimator_weights_format_average_regressor():
    """ Test estimator weights inputs as list and array """

    reg1 = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                            random_state=rng,
                            n_estimators=10)
    reg2 = RandomForestRegressor(random_state=rng,
                                 n_estimators=10)
    reg3 = GradientBoostingRegressor(random_state=rng,
                                     n_estimators=10)

    ensemble1 = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)],
        weights=[1, 2, 10])

    ensemble2 = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)],
        weights=np.array([1, 2, 10]))

    ensemble1.fit(X_train, y_train)
    ensemble2.fit(X_train, y_train)

    assert_array_almost_equal(ensemble1.predict(X_test),
                              ensemble2.predict(X_test))


def test_sample_weight_kwargs_average_regressor():
    """Check that AverageRegressor passes sample_weight as kwargs"""
    class MockRegressor(BaseEstimator, RegressorMixin):
        """Mock Regressor to check that sample_weight is received as kwargs"""
        def fit(self, X, y, *args, **sample_weight):
            assert_true('sample_weight' in sample_weight)

    clf = MockRegressor()
    eclf = AverageRegressor(estimators=[('mock', clf)])

    # Should not raise an error.
    eclf.fit(X_train, y_train, sample_weight=np.ones((len(y_train),)))


def test_sample_weight_average_regressor():
    """Tests sample_weight parameter of AverageRegressor"""
    reg1 = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                            random_state=rng,
                            n_estimators=10)
    reg2 = RandomForestRegressor(random_state=rng,
                                 n_estimators=10)
    reg3 = GradientBoostingRegressor(random_state=rng,
                                     n_estimators=10)

    ensemble1 = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)]
                    ).fit(X_train,
                          y_train,
                          sample_weight=np.ones((len(y_train),))
                          )
    ensemble2 = AverageRegressor(
        estimators=[('br', reg1),
                    ('rf', reg2),
                    ('gbr', reg3)]
                    ).fit(X_train,
                          y_train
                          )

    assert_array_equal(ensemble1.predict(X_test), ensemble2.predict(X_test))

    sample_weight = np.random.RandomState(123).uniform(size=(len(y_train),))

    ensemble3 = AverageRegressor(
        estimators=[('gbr', reg3)])

    ensemble3.fit(X_train, y_train, sample_weight)
    reg3.fit(X_train, y_train, sample_weight)
    assert_array_equal(ensemble3.predict(X_test), reg3.predict(X_test))

    knnR = KNeighborsRegressor()
    ensemble4 = AverageRegressor(estimators=[
        ('rf', reg2), ('knn', knnR)])

    msg = ('Underlying estimator \'knn\' does not support sample weights.')
    assert_raise_message(
        ValueError, msg, ensemble4.fit, X_train, y_train, sample_weight)


def test_average_regression_with_pipeline():

    pipline = make_pipeline(SelectKBest(k=1), DecisionTreeRegressor())
    estimator = AverageRegressor([('pp', pipline)])

    estimator.fit(X_train, y_train)
    assert_true(hasattr(estimator, 'estimators_'))


def test_set_params_average_regressor():
    """set_params should be able to set estimators"""
    reg1 = LinearRegression()
    reg2 = RandomForestRegressor(random_state=rng,
                                 n_estimators=10)
    reg3 = GradientBoostingRegressor(random_state=rng,
                                     n_estimators=10)

    ensmble = AverageRegressor(estimators=[
        ('lr', reg1), ('rf', reg2)], weights=[1, 2])

    assert_true('lr' in ensmble.named_estimators)
    assert_true(ensmble.named_estimators.lr is ensmble.estimators[0][1])
    assert_true(ensmble.named_estimators.lr is ensmble.named_estimators['lr'])

    ensmble.fit(X_train, y_train)
    assert_true('lr' in ensmble.named_estimators_)
    assert_true(ensmble.named_estimators_.lr is ensmble.estimators_[0])
    assert_true(ensmble.named_estimators_.lr is
                ensmble.named_estimators_['lr'])

    ensmble2 = AverageRegressor([
        ('lr', reg1), ('gbr', reg3)], weights=[1, 2])
    ensmble2.set_params(gbr=reg2).fit(X_train, y_train)
    assert_false(hasattr(ensmble2, 'gbr'))

    assert_array_equal(ensmble.predict(X_train), ensmble2.predict(X_train))
    assert_equal(ensmble2.estimators[0][1].get_params(), reg1.get_params())
    assert_equal(ensmble2.estimators[1][1].get_params(), reg2.get_params())

    ensmble2.set_params(gbr__max_depth=5)
    assert_true(ensmble2.estimators[1][1].get_params()['max_depth'] == 5)

    assert_equal(ensmble2.get_params()["gbr__max_depth"],
                 ensmble2.get_params()["gbr"].get_params()['max_depth'])
