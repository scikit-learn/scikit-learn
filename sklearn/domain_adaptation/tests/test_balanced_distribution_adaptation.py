from math import ceil

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from sklearn.ensemble import StackingClassifier
from sklearn.exceptions import NotFittedError
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_iris, make_blobs
from sklearn.metrics import accuracy_score

from sklearn.domain_adaptation import BalancedDistributionAdaptation

# Author: Oliver Rausch <rauscho@ethz.ch>
# License: BSD 3 clause

# load the iris dataset and randomly permute it
iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                    iris.target,
                                                    random_state=0)


@pytest.mark.parametrize("dim, lamb, mu",
                         [(-1, 1.0), (-100, -2), (-10, 10)])
def test_invalid_params(dim, lamb, mu):
    # Test negative iterations
    tca = BalancedDistributionAdaptation(dim=dim)
    with pytest.raises(ValueError, match="dim must be >= 0 or None"):
        tca.fit(X_train, X_test)

    tca = BalancedDistributionAdaptation(lamb=lamb)
    with pytest.raises(ValueError, match="lamb must be >= 0"):
        tca.fit(X_train, X_test)

    tca = BalancedDistributionAdaptation(mu=mu)
    with pytest.raises(ValueError, match="mu must be >= 0 and <= 1"):
        tca.fit(X_train, X_test)
