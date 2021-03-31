import numpy as np
from sklearn.neighbors import KNeighborsClassifier
import pytest
from sklearn.domain_adaptation import BalancedDistributionAdaptation

# Author: Jindong Wang <jindongwang@outlook.com>
# License: BSD 3 clause

Xs, Xt = np.random.randn(100, 50), np.random.randn(90, 50)
Ys, Yt = np.random.randint(0, 3, 100), np.random.randint(0, 3, 90)


def test_basic_output():
    # expected_output = {
    #     0.2227,
    #     'ypre should be 90 x 1 array'
    # }
    clf = KNeighborsClassifier()
    tca = BalancedDistributionAdaptation()
    acc, ypre, _, _ = tca.fit_predict(Xs, Ys, Xt, Yt, clf)
    print(acc, ypre)


@pytest.mark.parametrize("dim, lamb",
                         [(-1, 1.0), (-100, -2), (-10, 10)])
def test_invalid_params(dim, lamb):
    # Test negative iterations
    bda = BalancedDistributionAdaptation(dim=dim)
    with pytest.raises(ValueError, match="dim must be >= 0 or None"):
        bda.fit(Xs, Xt)

    bda = BalancedDistributionAdaptation(lamb=lamb)
    with pytest.raises(ValueError, match="lamb must be >= 0 or None"):
        bda.fit(Xs, Xt)
