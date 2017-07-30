import pickle
import numpy as np
from sklearn import datasets
from sklearn.freeze import FreezeWrap
from sklearn.feature_selection import SelectKBest
from sklearn.tree import DecisionTreeClassifier
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false


def test_freeze():
    X, y = datasets.load_iris(return_X_y=True)

    est = SelectKBest(k=1).fit(X, y)

    frozen_est = FreezeWrap(est)

    dumped = pickle.dumps(frozen_est)
    frozen_est2 = pickle.loads(dumped)
    assert_false(frozen_est is frozen_est2)

    # Test fit_transform where expected
    assert_true(hasattr(est, 'fit_transform'))
    assert_true(hasattr(frozen_est, 'fit_transform'))
    assert_false(est.fit_transform is frozen_est.fit_transform)
    frozen_est.fit_transform([np.arange(X.shape[1])], [0])

    # Test fit_transform not available when not on base
    est = DecisionTreeClassifier().fit(X, y)
    frozen_est = FreezeWrap(est)
    assert_false(hasattr(est, 'fit_transform'))
    assert_false(hasattr(frozen_est, 'fit_transform'))

    # TODO: much more
