import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false

from sklearn.decomposition.mf import MatrixFactorization


def test_fit_transform():
    random_state = np.random.mtrand.RandomState(0)

    n_samples = 100
    n_features = 25
    rank = 15

    L = random_state.rand(n_samples, rank)
    R = random_state.rand(rank, n_features)

    # Keep 20% of X to train and test
    X = np.dot(L, R)
    X_train = X.copy()
    X_train[random_state.rand(n_samples, n_features) > 0.8] = np.nan
    X_train_bak = X_train.copy()

    # Test with different numbers of components
    for n_components in [None, 10, 30]:
        X_train = X_train_bak.copy()
        mf = MatrixFactorization(n_components=n_components, random_state=0)
        L = mf.fit_transform(X_train)

        # Check the format of the output
        assert_equal(L.shape[0], n_samples)
        assert_equal(mf.components_.shape[1], n_features)
        if n_components is not None:
            assert_true(L.shape[1], n_components + 2)
            assert_true(mf.components_.shape[0], n_components + 2)

        # Test with dense and sparse input
        L_all = []
        R_all = []

        mf = MatrixFactorization(n_components=n_components, random_state=0)
        L_all.append(mf.fit_transform(X_train))
        R_all.append(mf.components_)

        X_train[np.isnan(X_train)] = 0
        mf = MatrixFactorization(n_components=n_components, random_state=0,
                                 missing_values=0)
        L_all.append(mf.fit_transform(X_train))
        R_all.append(mf.components_)

        assert_array_equal(L_all[0], L_all[1])
        assert_array_equal(R_all[0], R_all[1])

        mf = MatrixFactorization(n_components=n_components, random_state=0)
        L_all.append(mf.fit_transform(sp.coo_matrix(X_train)))
        R_all.append(mf.components_)

        assert_array_equal(L_all[1], L_all[2])
        assert_array_equal(R_all[0], R_all[1])

        # Test that no exception is raised with negative parameters
        mf = MatrixFactorization(
            n_components=n_components,
            learning_rate=-1,
            regularization=-1)
        L = mf.fit_transform(sp.coo_matrix(X_train))