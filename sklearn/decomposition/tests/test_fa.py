import numpy as np

from .. import FA


def test_fa_generativ():
    """Factor Analysis generative
    story: how well does it model
    some input?
    """

    # Some random settings for the generative model
    W = np.random.randn(3, 5)
    # latent variable of dim 3, 20 of it
    h = np.random.randn(20, 3)
    # using gamma to model different noise variance
    # per component
    noise = np.random.gamma(1, size=5) * np.random.randn(20, 5)

    # generate observations
    # wlog, mean is 0
    X = np.dot(h, W) + noise

    fa = FA(n_components=3)
    fa.fit(X)

    data = X
    W = fa.W
    psi = fa.psi
    # Sample Covariance
    scov = np.cov(data, rowvar=0, bias=1)
    # Model Covariance
    mcov = np.dot(W.T, W) + np.diag(psi)
    diff =  np.sum(np.abs(scov-mcov))/W.size
    assert diff < 0.1, "Mean absolute difference is %f" % diff
