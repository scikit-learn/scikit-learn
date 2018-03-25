"""
Testing for the AverageRegressor module (sklearn.ensemble.AverageRegressor).
"""

# Author: Mohamed Ali Jamaoui m.ali.jamaoui@gmail.com
#
# License: BSD 3 clause

import numpy as np
from sklearn import datasets
from sklearn.ensemble import AverageRegressor
from sklearn.utils.validation import check_random_state
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.exceptions import NotFittedError
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest


rng = check_random_state(0)
# load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                    boston.target,
                                                    random_state=rng)


def test_estimator_init():
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


def test_notfitted():
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

    ensemble.fit(X_train, y_train)
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


def test_gridsearch():
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


def test_estimator_weights_format():
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


def test_sample_weight_kwargs():
    """Check that AverageRegressor passes sample_weight as kwargs"""
    class MockRegressor(BaseEstimator, RegressorMixin):
        """Mock Regressor to check that sample_weight is received as kwargs"""
        def fit(self, X, y, *args, **sample_weight):
            assert_true('sample_weight' in sample_weight)

    clf = MockRegressor()
    eclf = AverageRegressor(estimators=[('mock', clf)])

    # Should not raise an error.
    eclf.fit(X_train, y_train, sample_weight=np.ones((len(y_train),)))


def test_sample_weight():
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


def test_set_params():
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
