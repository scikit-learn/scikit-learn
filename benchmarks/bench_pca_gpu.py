"""
A comparison of PCA with/wthout GPU

Obtained with NVIDIA Tesla V100

With svd_solver='full' and iterated_power=2:
        Without GPU : 1.220
        With GPU : 0.727

With svd_solver='full' and iterated_power=10:
        Without GPU : 1.204
        With GPU : 0.013

With svd_solver='randomized' and iterated_power=2:
        Without GPU : 0.106
        With GPU : 0.641

With svd_solver='randomized' and iterated_power=10:
        Without GPU : 0.261
        With GPU : 0.276
"""

import sklearn
from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
import numpy as np
import cupy as cp
import time


if __name__ == '__main__':
    import cupy
    X_np, y_np = make_classification(n_samples=10000, n_features=100)
    X_np = X_np.astype(np.float32)
    X_cp = cp.asarray(X_np)

    for svd_solver in ["full", "randomized"]:
        for iterated_power in [2, 10]:
            pca_np = PCA(n_components=3, svd_solver=svd_solver, copy=True,
                         random_state=0, iterated_power=iterated_power)
            t0 = time.time()
            pca_np.fit_transform(X_np)
            without_gpu_time = time.time() - t0

            with_gpu_time = 0
            with sklearn.config_context(enable_duck_array=True):
                pca_cp = PCA(**pca_np.get_params())
                t0 = time.time()
                pca_cp.fit_transform(X_cp)
                with_gpu_time = time.time() - t0

            msg = 'With svd_solver=\'{}\' and iterated_power={}:'
            print(msg.format(svd_solver, iterated_power))
            print('\tWithout GPU : {:.3f}'.format(without_gpu_time))
            print('\tWith GPU : {:.3f}\n'.format(with_gpu_time))
