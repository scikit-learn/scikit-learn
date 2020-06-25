"""
A comparison of PCA with/wthout GPU

Obtained with NVIDIA Tesla V100

With svd_solver='full' and iterated_power=2:
    Without GPU : runtime: 1.194s, explained variance: 5.658
    With GPU : runtime: 0.012s, explained variance: 5.658

With svd_solver='full' and iterated_power=10:
    Without GPU : runtime: 1.195s, explained variance: 5.658
    With GPU : runtime: 0.012s, explained variance: 5.658

With svd_solver='randomized' and iterated_power=2:
    Without GPU : runtime: 0.033s, explained variance: 5.145
    With GPU : runtime: 0.006s, explained variance: 5.144

With svd_solver='randomized' and iterated_power=10:
    Without GPU : runtime: 0.263s, explained variance: 5.624
    With GPU : runtime: 0.016s, explained variance: 5.636
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
            pca_np.fit_transform(X_np)  # Warm-up
            pca_np = PCA(**pca_np.get_params())  # Overwrite to clean model
            t0 = time.time()
            pca_np.fit_transform(X_np)
            without_gpu_time = time.time() - t0
            exp_var_np = pca_np.explained_variance_.sum()

            with_gpu_time = 0
            exp_var_cp = 0
            with sklearn.config_context(enable_duck_array=True):
                pca_cp = PCA(**pca_np.get_params())
                pca_cp.fit_transform(X_cp)  # Warm-up
                pca_cp = PCA(**pca_np.get_params())  # Overwrite to clean model
                t0 = time.time()
                pca_cp.fit_transform(X_cp)
                with_gpu_time = time.time() - t0
                exp_var_cp = cp.asnumpy(pca_cp.explained_variance_.sum())

            msg = 'With svd_solver=\'{}\' and iterated_power={}:'
            print(msg.format(svd_solver, iterated_power))
            m1 = '\tWithout GPU : runtime: {:.3f}s, explained variance: {:.3f}'
            print(m1.format(without_gpu_time, exp_var_np))
            m2 = '\tWith GPU : runtime: {:.3f}s, explained variance: {:.3f}\n'
            print(m2.format(with_gpu_time, exp_var_cp))
