import numpy as np
from numpy.testing import assert_allclose
from sklearn.neighbors.kde import KernelDensity
from sklearn.neighbors.ball_tree import kernel_norm


def compute_kernel_slow(Y, X, kernel, h):
    d = np.sqrt(((Y[:, None, :] - X) ** 2).sum(-1))
    norm = kernel_norm(h, X.shape[1], kernel) / X.shape[0]

    if kernel == 'gaussian':
        return norm * np.exp(-0.5 * (d * d) / (h * h)).sum(-1)
    elif kernel == 'tophat':
        return norm * (d < h).sum(-1)
    elif kernel == 'epanechnikov':
        return norm * ((1.0 - (d * d) / (h * h)) * (d < h)).sum(-1)
    elif kernel == 'exponential':
        return norm * (np.exp(-d / h)).sum(-1)
    elif kernel == 'linear':
        return norm * ((1 - d / h) * (d < h)).sum(-1)
    elif kernel == 'cosine':
        return norm * (np.cos(0.5 * np.pi * d / h) * (d < h)).sum(-1)
    else:
        raise ValueError('kernel not recognized')


def test_KernelDensity(n_samples=100, n_features=3):
    np.random.seed(0)
    X = np.random.random((n_samples, n_features))
    Y = np.random.random((n_samples, n_features))

    for kernel in ['gaussian', 'tophat', 'epanechnikov',
                   'exponential', 'linear', 'cosine']:
        for bandwidth in [0.01, 0.1, 1]:
            dens_true = compute_kernel_slow(Y, X, kernel, bandwidth)
            def check_results(kernel, bandwidth, atol, rtol):
                kde = KernelDensity(kernel=kernel, bandwidth=bandwidth,
                                    atol=atol, rtol=rtol)
                log_dens = kde.fit(X).eval(Y)
                assert_allclose(np.exp(log_dens), dens_true, atol=atol,
                                rtol=max(1E-7, rtol))

            for rtol in [0, 1E-5]:
                for atol in [1E-6, 1E-2]:
                    for breadth_first in (True, False):
                        yield (check_results, kernel, bandwidth, atol, rtol)


if __name__ == '__main__':
    import nose
    nose.runmodule()
