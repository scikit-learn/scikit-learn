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
from sklearn.exceptions import NotFittedError

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression


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
