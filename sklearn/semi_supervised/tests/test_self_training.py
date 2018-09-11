import numpy as np

from sklearn.base import BaseEstimator

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_raise_message

from sklearn.semi_supervised import SelfTraining
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import ParameterGrid
from sklearn.linear_model import Perceptron
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_iris
from sklearn.utils import check_random_state

# Author: Oliver Rausch
# License: BSD 3 clause


rng = check_random_state(0)

# load the iris dataset and randomly permute it
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                    iris.target,
                                                    random_state=rng)

limit = 50

y_train_missing_labels = y_train.copy()
y_train_missing_labels[limit:] = -1

def test_classification():
    # Check classification for various parameter settings.
    grid = ParameterGrid({"max_iter": [1, 50, 100],
                          "threshold": [0.0, 0.5, 1.0]})

    for base_estimator in [DummyClassifier(),
                           DecisionTreeClassifier(),
                           KNeighborsClassifier(),
                           SVC(gamma="scale", probability=True)]:
        for params in grid:
            SelfTraining(base_estimator,
                         **params).fit(X_train, y_train_missing_labels).predict(X_test)


def test_missing_predict_proba():
    # Check that an error is thrown if predict_proba is not implemented
    base_estimator = SVC(gamma="scale")
    self_training = SelfTraining(base_estimator)
    message = "The base_estimator should implement predict_proba!"
    assert_raise_message(ValueError, message, self_training.fit, X_train,
                         y_train)

def test_single_estimator():
    # Check classification for single iteration.
    # Fitting a SelfTraining estimator with one iteration and 100 unlabeled 
    # datapoints should give the same results as fitting a normal classifier 
    # with only 50 labeled datapoints.

    clf1 = SelfTraining(base_estimator=KNeighborsClassifier(),
                            max_iter=0).fit(X_train, y_train_missing_labels)

    clf2 = KNeighborsClassifier().fit(X_train[:limit], y_train[:limit])

    assert_array_equal(clf1.predict(X_test), clf2.predict(X_test))
