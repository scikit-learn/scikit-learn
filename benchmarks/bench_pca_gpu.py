"""
A comparison of PCA with/wthout GPU

Obtained with NVIDIA Tesla V100

With svd_solver='full' and iterated_power=2:
        Without GPU : runtime: 1.197s, explained variance: 6.967
        With CuPy : runtime: 0.016s, explained variance: 6.967
        With cuML : runtime: 0.007s, explained variance: 6.967

With svd_solver='full' and iterated_power=10:
        Without GPU : runtime: 1.210s, explained variance: 6.967
        With CuPy : runtime: 0.016s, explained variance: 6.967
        With cuML : runtime: 0.007s, explained variance: 6.967

With svd_solver='randomized' and iterated_power=2:
        Without GPU : runtime: 0.032s, explained variance: 6.745
        With CuPy : runtime: 0.012s, explained variance: 6.693

With svd_solver='randomized' and iterated_power=10:
        Without GPU : runtime: 0.204s, explained variance: 6.945
        With CuPy : runtime: 0.040s, explained variance: 6.949
"""

import sklearn
from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
from cuml.decomposition import PCA as cumlPCA
import numpy as np
import cupy as cp
import time


if __name__ == '__main__':
    import cupy
    X_np, y_np = make_classification(n_samples=10000, n_features=100)
    X_np = X_np.astype(np.float32)
    X_cp = cp.asfortranarray(cp.asarray(X_np))

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

            cupy_time = 0
            exp_var_cp = 0
            with sklearn.config_context(enable_duck_array=True):
                pca_cp = PCA(**pca_np.get_params())
                pca_cp.fit_transform(X_cp)  # Warm-up
                pca_cp = PCA(**pca_np.get_params())  # Overwrite to clean model
                t0 = time.time()
                pca_cp.fit_transform(X_cp)
                cupy_time = time.time() - t0
                exp_var_cp = cp.asnumpy(pca_cp.explained_variance_.sum())

            if svd_solver == 'full':
                pca_cuml = cumlPCA(**pca_np.get_params())
                pca_cuml.fit_transform(X_cp)  # Warm-up
                pca_cuml = cumlPCA(**pca_np.get_params())  # Overwrite to clean model
                t0 = time.time()
                pca_cuml.fit_transform(X_cp)
                cuml_time = time.time() - t0
                exp_var_cuml = pca_np.explained_variance_.sum()

            msg = 'With svd_solver=\'{}\' and iterated_power={}:'
            print(msg.format(svd_solver, iterated_power))
            m1 = '\tWithout GPU : runtime: {:.3f}s, explained variance: {:.3f}'
            print(m1.format(without_gpu_time, exp_var_np))
            m2 = '\tWith CuPy : runtime: {:.3f}s, explained variance: {:.3f}'
            print(m2.format(cupy_time, exp_var_cp))

            m3 = '\tWith cuML : runtime: {:.3f}s, explained variance: {:.3f}\n'
            if svd_solver == 'full':
                print(m3.format(cuml_time, exp_var_cuml))
            else:
                print()