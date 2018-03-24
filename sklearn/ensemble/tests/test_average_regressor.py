"""
Testing for the AverageRegressor module (sklearn.ensemble.AverageRegressor).
"""

# Author: Mohamed Ali Jamaoui m.ali.jamaoui@gmail.com
#
# License: BSD 3 clause

from sklearn import datasets
from sklearn.ensemble import AverageRegressor
from sklearn.utils.validation import check_random_state
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.exceptions import NotFittedError

from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import make_scorer


rng = check_random_state(0)
# load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_estimator_init():
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)

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
    msg = ('Number of classifiers and weights must be equal'
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
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)

    ensemble = AverageRegressor(estimators=[('lr1', LinearRegression()),
                                        ('lr2', LinearRegression())])
    msg = ("This AverageRegressor instance is not fitted yet. Call \'fit\'"
           " with appropriate arguments before using this method.")
    assert_raise_message(NotFittedError, msg, ensemble.predict, X_train)


def test_single_estimator():
    # Check singleton ensembles.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    clf1 = AverageRegressor(estimators=[
        ('knn', KNeighborsRegressor())]).fit(X_train, y_train)


    clf2 = KNeighborsRegressor().fit(X_train, y_train)

    assert_array_almost_equal(clf1.predict(X_test), clf2.predict(X_test))


def test_average_prediction():
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

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
    
    assert_almost_equal(scores.mean(), 0.85, decimal=2)


def test_parallel_regression():
    # Check parallel regression.
    rng = check_random_state(0)

    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

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
