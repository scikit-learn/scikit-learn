"""
===========================================
Sparse coding with a precomputed dictionary
===========================================

Transform a signal as a sparse combination of Ricker wavelets. This example
visually compares different sparse coding methods using the
:class:`sklearn.decomposition.SparseCoder` estimator.
"""
print __doc__

import numpy as np
import matplotlib.pylab as pl

from sklearn.decomposition import SparseCoder


def ricker_function(resolution, center, width):
    """Discrete sub-sampled Ricker (mexican hat) wavelet"""
    x = np.linspace(0, resolution - 1, resolution)
    x = (2 / ((np.sqrt(3 * width) * np.pi ** 1 / 4))) * (
         1 - ((x - center) ** 2 / width ** 2)) * np.exp(
         (-(x - center) ** 2) / (2 * width ** 2))
    return x


def ricker_matrix(width, resolution, n_atoms):
    """Dictionary of Ricker (mexican hat) wavelets"""
    centers = np.linspace(0, resolution - 1, n_atoms)
    D = np.empty((n_atoms, resolution))
    for i, center in enumerate(centers):
        D[i] = ricker_function(resolution, center, width)
    D /= np.sqrt(np.sum(D ** 2, axis=1))[:, np.newaxis]
    return D


resolution = 1024
subsampling = 3  # subsampling factor
width = 10
n_atoms = resolution / subsampling

# Compute a wavelet dictionary
D = ricker_matrix(width=width, resolution=resolution, n_atoms=n_atoms)

# Generate a signal
y = np.linspace(0, resolution - 1, resolution)
y = np.sin(y / 30) + np.sin(y / 15)

# List the different sparse coding methods in the following format:
# (title, transform_algorithm, transform_alpha, transform_n_nozero_coefs)
estimators = [('Orthogonal Matching Pursuit', 'omp', None, 32),
              ('Least-angle Regression', 'lars', None, 32),
              ('Lasso', 'lasso_cd', 2., None),
              ('Soft thresholding', 'threshold', 2.5, None)]

# Do a wavelet approximation
for title, algo, alpha, n_nonzero in estimators:
    coder = SparseCoder(dictionary=D, transform_n_nonzero_coefs=n_nonzero,
                        transform_alpha=alpha, transform_algorithm=algo)
    x = coder.transform(y)
    pl.figure()
    density = len(np.flatnonzero(x))
    pl.title(title + ': %s nonzero coefs' % density)
    pl.plot(y, ls='dotted')
    pl.plot(np.ravel(np.dot(x, D)))
    pl.axis('tight')

# Soft thresholding debiasing
_, idx = np.where(x != 0)
x[0, idx], _, _, _ = np.linalg.lstsq(D[idx, :].T, y)
pl.figure()
pl.title('Soft thresholding w/ debiasing')
pl.plot(y, ls='dotted')
pl.plot(np.ravel(np.dot(x, D)))
pl.axis('tight')

pl.show()
