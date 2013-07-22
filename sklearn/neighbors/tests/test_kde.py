import numpy as np
from numpy.testing import assert_allclose, assert_raises, assert_equal
from sklearn.neighbors import KernelDensity, BallTree, KDTree, NearestNeighbors
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
                assert_allclose(np.exp(kde.score(Y)),
                                np.prod(dens_true),
                                atol=atol, rtol=max(1E-7, rtol))

            for rtol in [0, 1E-5]:
                for atol in [1E-6, 1E-2]:
                    for breadth_first in (True, False):
                        yield (check_results, kernel, bandwidth, atol, rtol)


def test_KernelDensity_sampling(n_samples=100, n_features=3):
    np.random.seed(0)
    X = np.random.random((n_samples, n_features))

    bandwidth = 0.2

    for kernel in ['gaussian', 'tophat']:
        # draw a tophat sample
        kde = KernelDensity(bandwidth, kernel=kernel).fit(X)
        samp = kde.sample(100)
        assert_equal(X.shape, samp.shape)

        # check that samples are in the right range
        nbrs = NearestNeighbors(n_neighbors=1).fit(X)
        dist, ind = nbrs.kneighbors(X, return_distance=True)
        
        if kernel == 'tophat':
            assert np.all(dist < bandwidth)
        elif kernel == 'gaussian':
            # 5 standard deviations is safe for 100 samples, but there's a
            # very small chance this test could fail.
            assert np.all(dist < 5 * bandwidth)

    # check unsupported kernels
    for kernel in ['epanechnikov', 'exponential', 'linear', 'cosine']:
        kde = KernelDensity(bandwidth, kernel=kernel).fit(X)
        assert_raises(NotImplementedError, kde.sample, 100)
    
    


def test_kde_algorithm_metric_choice():
    """Smoke test for various metrics and algorithms"""
    np.random.seed(0)
    X = np.random.random((10, 2))  # 2 features required for haversine dist.
    Y = np.random.random((10, 2))

    for algorithm in ['auto', 'ball_tree', 'kd_tree']:
        for metric in ['euclidean', 'minkowski', 'manhattan',
                       'chebyshev', 'haversine']:
            if algorithm == 'kd_tree' and metric not in KDTree.valid_metrics:
                assert_raises(ValueError, KernelDensity,
                              algorithm=algorithm, metric=metric)
            else:
                kde = KernelDensity(algorithm=algorithm, metric=metric)
                kde.fit(X)
                y_dens = kde.eval(Y)
                assert_equal(y_dens.shape, Y.shape[:1])

def test_kde_score(n_samples=100, n_features=3):
    np.random.seed(0)
    X = np.random.random((n_samples, n_features))
    Y = np.random.random((n_samples, n_features))

def test_kde_badargs():
    np.random.seed(0)
    X = np.random.random((10, 2))
    Y = np.random.random((10, 2))

    assert_raises(ValueError, KernelDensity,
                  algorithm='blah')
    assert_raises(ValueError, KernelDensity,
                  bandwidth=0)
    assert_raises(ValueError, KernelDensity,
                  kernel='blah')
    assert_raises(ValueError, KernelDensity,
                  metric='blah')
    assert_raises(ValueError, KernelDensity,
                  algorithm='kd_tree', metric='blah')

if __name__ == '__main__':
    import nose
    nose.runmodule()
