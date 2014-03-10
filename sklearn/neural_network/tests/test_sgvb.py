import sys
import re

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from sklearn.utils.testing import (assert_almost_equal, assert_array_equal,
                                   assert_true)

from sklearn.datasets import load_digits
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.neural_network import SGVB
from sklearn.utils.validation import assert_all_finite

np.seterr(all='warn')

Xdigits = load_digits().data
Xdigits -= Xdigits.min()
Xdigits /= Xdigits.max()

# To do:
# Do same tests with continuous data. Olivetti faces?

def test_fit():
    X = Xdigits.copy()

    sgvb = SGVB(n_components_decoder=64, n_components_encoder=64, n_hidden_variables=10, learning_rate=0.1,
                   batch_size=10, n_iter=20, sampling_rounds=2, random_state=9)
    sgvb.fit(X)

    assert_almost_equal(sgvb.score(X), -25., decimal=0)

    # in-place tricks shouldn't have modified X
    assert_array_equal(X, Xdigits)


def test_transform():
    X = Xdigits[:100]
    sgvb = SGVB(n_components_decoder=16, n_components_encoder=16, batch_size=5,
                        n_iter=5, random_state=42)
    sgvb.fit(X)

    Xt1 = sgvb.transform(X)
    Xt2 = np.tanh(X.dot(sgvb.params["W1"].T) + sgvb.params["b1"].T)

    assert_array_equal(Xt1, Xt2)

def test_fit_transform():
    X = Xdigits[:100]
    sgvb = SGVB(n_components_decoder=16, n_components_encoder=16, batch_size=5,
                        n_iter=5, random_state=42)

    Xt1 = sgvb.fit_transform(X)
    Xt2 = np.tanh(X.dot(sgvb.params["W1"].T) + sgvb.params["b1"].T)

    assert_array_equal(Xt1, Xt2)


def test_sgvb_verbose():
    sgvb = SGVB(n_iter=2, verbose=10)
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        sgvb.fit(Xdigits)
    finally:
        sys.stdout = old_stdout



