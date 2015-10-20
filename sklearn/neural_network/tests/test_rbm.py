import sys
import re

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix
from sklearn.utils.testing import (assert_almost_equal, assert_array_equal,
                                   assert_true)

from sklearn.datasets import load_digits
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.neural_network import BernoulliRBM
from sklearn.utils.validation import assert_all_finite
np.seterr(all='warn')

Xdigits = load_digits().data
Xdigits -= Xdigits.min()
Xdigits /= Xdigits.max()


def test_fit():
    X = Xdigits.copy()

    rbm = BernoulliRBM(n_components=64, learning_rate=0.1,
                       batch_size=10, n_iter=7, random_state=9)
    rbm.fit(X)

    assert_almost_equal(rbm.score_samples(X).mean(), -21., decimal=0)

    # in-place tricks shouldn't have modified X
    assert_array_equal(X, Xdigits)


def test_partial_fit():
    X = Xdigits.copy()
    rbm = BernoulliRBM(n_components=64, learning_rate=0.1,
                       batch_size=20, random_state=9)
    n_samples = X.shape[0]
    n_batches = int(np.ceil(float(n_samples) / rbm.batch_size))
    batch_slices = np.array_split(X, n_batches)

    for i in range(7):
        for batch in batch_slices:
            rbm.partial_fit(batch)

    assert_almost_equal(rbm.score_samples(X).mean(), -21., decimal=0)
    assert_array_equal(X, Xdigits)


def test_transform():
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=16, batch_size=5,
                        n_iter=5, random_state=42)
    rbm1.fit(X)

    Xt1 = rbm1.transform(X)
    Xt2 = rbm1._mean_hiddens(X)

    assert_array_equal(Xt1, Xt2)


def test_small_sparse():
    # BernoulliRBM should work on small sparse matrices.
    X = csr_matrix(Xdigits[:4])
    BernoulliRBM().fit(X)       # no exception


def test_small_sparse_partial_fit():
    for sparse in [csc_matrix, csr_matrix]:
        X_sparse = sparse(Xdigits[:100])
        X = Xdigits[:100].copy()

        rbm1 = BernoulliRBM(n_components=64, learning_rate=0.1,
                            batch_size=10, random_state=9)
        rbm2 = BernoulliRBM(n_components=64, learning_rate=0.1,
                            batch_size=10, random_state=9)

        rbm1.partial_fit(X_sparse)
        rbm2.partial_fit(X)

        assert_almost_equal(rbm1.score_samples(X).mean(),
                            rbm2.score_samples(X).mean(),
                            decimal=0)


def test_sample_hiddens():
    rng = np.random.RandomState(0)
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=2, batch_size=5,
                        n_iter=5, random_state=42)
    rbm1.fit(X)

    h = rbm1._mean_hiddens(X[0])
    hs = np.mean([rbm1._sample_hiddens(X[0], rng) for i in range(100)], 0)

    assert_almost_equal(h, hs, decimal=1)


def test_fit_gibbs():
    # Gibbs on the RBM hidden layer should be able to recreate [[0], [1]]
    # from the same input
    rng = np.random.RandomState(42)
    X = np.array([[0.], [1.]])
    rbm1 = BernoulliRBM(n_components=2, batch_size=2,
                        n_iter=42, random_state=rng)
    # you need that much iters
    rbm1.fit(X)
    assert_almost_equal(rbm1.components_,
                        np.array([[0.02649814], [0.02009084]]), decimal=4)
    assert_almost_equal(rbm1.gibbs(X), X)
    return rbm1


def test_fit_gibbs_sparse():
    # Gibbs on the RBM hidden layer should be able to recreate [[0], [1]] from
    # the same input even when the input is sparse, and test against non-sparse
    rbm1 = test_fit_gibbs()
    rng = np.random.RandomState(42)
    from scipy.sparse import csc_matrix
    X = csc_matrix([[0.], [1.]])
    rbm2 = BernoulliRBM(n_components=2, batch_size=2,
                        n_iter=42, random_state=rng)
    rbm2.fit(X)
    assert_almost_equal(rbm2.components_,
                        np.array([[0.02649814], [0.02009084]]), decimal=4)
    assert_almost_equal(rbm2.gibbs(X), X.toarray())
    assert_almost_equal(rbm1.components_, rbm2.components_)


def test_gibbs_smoke():
    # Check if we don't get NaNs sampling the full digits dataset.
    # Also check that sampling again will yield different results.
    X = Xdigits
    rbm1 = BernoulliRBM(n_components=42, batch_size=40,
                        n_iter=20, random_state=42)
    rbm1.fit(X)
    X_sampled = rbm1.gibbs(X)
    assert_all_finite(X_sampled)
    X_sampled2 = rbm1.gibbs(X)
    assert_true(np.all((X_sampled != X_sampled2).max(axis=1)))


def test_score_samples():
    # Test score_samples (pseudo-likelihood) method.
    # Assert that pseudo-likelihood is computed without clipping.
    # See Fabian's blog, http://bit.ly/1iYefRk
    rng = np.random.RandomState(42)
    X = np.vstack([np.zeros(1000), np.ones(1000)])
    rbm1 = BernoulliRBM(n_components=10, batch_size=2,
                        n_iter=10, random_state=rng)
    rbm1.fit(X)
    assert_true((rbm1.score_samples(X) < -300).all())

    # Sparse vs. dense should not affect the output. Also test sparse input
    # validation.
    rbm1.random_state = 42
    d_score = rbm1.score_samples(X)
    rbm1.random_state = 42
    s_score = rbm1.score_samples(lil_matrix(X))
    assert_almost_equal(d_score, s_score)

    # Test numerical stability (#2785): would previously generate infinities
    # and crash with an exception.
    with np.errstate(under='ignore'):
        rbm1.score_samples([np.arange(1000) * 100])


def test_rbm_verbose():
    rbm = BernoulliRBM(n_iter=2, verbose=10)
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        rbm.fit(Xdigits)
    finally:
        sys.stdout = old_stdout


def test_sparse_and_verbose():
    # Make sure RBM works with sparse input when verbose=True
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    from scipy.sparse import csc_matrix
    X = csc_matrix([[0.], [1.]])
    rbm = BernoulliRBM(n_components=2, batch_size=2, n_iter=1,
                       random_state=42, verbose=True)
    try:
        rbm.fit(X)
        s = sys.stdout.getvalue()
        # make sure output is sound
        assert_true(re.match(r"\[BernoulliRBM\] Iteration 1,"
                             r" pseudo-likelihood = -?(\d)+(\.\d+)?,"
                             r" time = (\d|\.)+s",
                             s))
    finally:
        sys.stdout = old_stdout
