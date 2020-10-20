# Author: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

import pytest
import numpy as np

from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture


@pytest.mark.parametrize(
    "estimator",
    [GaussianMixture(),
     BayesianGaussianMixture()]
)
def test_gaussian_mixture_n_iter(estimator):
    # check that n_iter is the number of iteration performed.
    rng = np.random.RandomState(0)
    X = rng.rand(10, 5)
    max_iter = 1
    estimator.set_params(max_iter=max_iter)
    estimator.fit(X)
    assert estimator.n_iter_ == max_iter
