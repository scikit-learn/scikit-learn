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
    assert_array_equal(est.scores_, frozen_est2.scores_)

    # scores should be unaffected by new fit
    assert_true(frozen_est2.fit() is frozen_est2)
    assert_array_equal(est.scores_, frozen_est2.scores_)

    # Test fit_transform where expected
    assert_true(hasattr(est, 'fit_transform'))
    assert_true(hasattr(frozen_est, 'fit_transform'))
    assert_false(est.fit_transform is frozen_est.fit_transform)
    frozen_est.fit_transform([np.arange(X.shape[1])], [0])
    # scores should be unaffected by new fit_transform
    assert_array_equal(est.scores_, frozen_est.scores_)

    # Test fit_transform not set when not needed
    est = DecisionTreeClassifier().fit(X, y)
    frozen_est = FreezeWrap(est)
    assert_false(hasattr(est, 'fit_transform'))
    assert_false(hasattr(frozen_est, 'fit_transform'))
