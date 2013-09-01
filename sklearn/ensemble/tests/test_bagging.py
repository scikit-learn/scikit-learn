"""
Testing for the ensemble module (sklearn.ensemble.ensemble).
"""

# Authors: Maheshakya Wijewardena
# License: BSD 3 clause

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_true
#from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import SkipTest

from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.grid_search import GridSearchCV, ParameterGrid
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.ensemble import _partition_estimators
from sklearn.linear_model import Perceptron, LogisticRegression
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC, SVR

from sklearn.cross_validation import train_test_split
from sklearn.datasets import load_boston, load_iris
from sklearn.utils import check_random_state

rng = check_random_state(0)

# also load the iris dataset
# and randomly permute it
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

def test__partition_estimators():
    """test for parttition_estimators() function"""
    ensemble = RandomForestClassifier(n_estimators=13, n_jobs= 5)
    
    assert_equal((5,[3,3,3,2,2],[0,3,6,9,11,13]),_partition_estimators(ensemble))
    



if __name__ == "__main__":
    import nose
    nose.runmodule()
