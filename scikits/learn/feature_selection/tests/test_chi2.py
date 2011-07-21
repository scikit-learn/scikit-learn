import numpy as np
from numpy.testing import assert_, assert_equal
from scipy.sparse import csr_matrix

from .. import Chi2

# Feature 0 is highly informative for class 1;
# feature 1 is the same everywhere;
# feature 2 is a bit informative for class 2.
X = ([[2, 1, 2],
      [9, 1, 1],
      [6, 1, 2],
      [0, 1, 3]])
y = [0, 1, 2, 2]


def test_chi2():
    """Test Chi2 feature extraction"""

    chi2 = Chi2(n_features=1).fit(X, y)
    assert_equal(chi2.top_features_, [0])
    assert_equal(chi2.transform(X), np.array(X)[:, [0]])

    chi2 = Chi2(n_features=2).fit(X, y)
    assert_equal(sorted(chi2.top_features_), [0, 2])

    Xsp = csr_matrix(X, dtype=np.float)
    chi2 = Chi2(n_features=2).fit(Xsp, y)
    assert_equal(sorted(chi2.top_features_), [0, 2])
    Xtrans = chi2.transform(Xsp)
    assert_equal(Xtrans.shape, [Xsp.shape[0], 2])

    # == doesn't work on scipy.sparse matrices
    Xtrans = Xtrans.toarray()
    Xtrans2 = Chi2(n_features=2).fit_transform(Xsp, y).toarray()
    assert_equal(Xtrans, Xtrans2)


def test_chi2_boolean():
    """Test Chi2 on boolean data
    
    Or, assert that X.sum() on a boolean array yields an int array."""

    #Xb = np.asanyarray(X) > 2
    Xb = np.array([[0, 0, 0],
                   [0, 0, 1],
                   [1, 0, 1],
                   [1, 0, 0]], dtype=bool)
    chi2 = Chi2(n_features=1).fit(Xb, y)
    assert_equal(chi2.top_features_, [0])
