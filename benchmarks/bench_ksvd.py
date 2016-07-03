"""Benchmark for the K-SVD function."""


from time import time

import numpy as np
import scipy.misc

from sklearn.decomposition import dict_learning
from sklearn.feature_extraction.image import extract_patches_2d
from sklearn.utils.random import choice


RANDOM_STATE = np.random.RandomState(0)


def compute_ksvd(patch_size, samples_count, n_components, n_nonzero_coefs,
                 max_iter, approximate_svd):
    image = scipy.misc.lena() / 256.0
    image_distorted = image + 0.1 * RANDOM_STATE.randn(*image.shape)

    patch_size = (patch_size, patch_size)
    data_full = extract_patches_2d(image_distorted, patch_size)
    data_full = data_full.reshape(data_full.shape[0], -1)
    data_full -= np.mean(data_full, axis=0)
    data_full /= np.std(data_full, axis=0)
    samples_indexes = choice(range(data_full.shape[0]), replace=False,
                             size=samples_count, random_state=RANDOM_STATE)
    data = data_full[samples_indexes, :]

    t0 = time()
    if approximate_svd:
        method = 'ksvd'
    else:
        method = 'exact_ksvd'

    dict_learning(data, n_components, n_nonzero_coefs=n_nonzero_coefs,
                  max_iter=max_iter, method=method)
    dt = time() - t0
    return dt


BENCHMARKS_PARAMETERS = [
    dict(patch_size=4, samples_count=5000, n_components=50,
         n_nonzero_coefs=2, max_iter=25, approximate_svd=True),
    dict(patch_size=4, samples_count=5000, n_components=50,
         n_nonzero_coefs=2, max_iter=25, approximate_svd=False),
    dict(patch_size=4, samples_count=10000, n_components=100,
         n_nonzero_coefs=2, max_iter=25, approximate_svd=True),
    dict(patch_size=6, samples_count=10000, n_components=100,
         n_nonzero_coefs=2, max_iter=25, approximate_svd=True),
]


def main():
    for parameters in BENCHMARKS_PARAMETERS:
        runtime = compute_ksvd(**parameters)
        print("Running ", parameters)
        print("Time: {:.2f}s".format(runtime))
        print()


if __name__ == '__main__':
    main()
