"""
=======================================
Kernel PCA Solvers comparison benchmark
=======================================

This example shows that the various solvers provided in Kernel PCA can help
drastically improve its execution speed when an approximate solution is
sufficient.
"""
print(__doc__)

# Authors: Sylvain MARIE
# License: BSD 3 clause

from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.decomposition import KernelPCA
from sklearn.datasets import make_circles


# 1- Generate random data
# -----------------------
np.random.seed(0)
n_train, n_test, n_features = 2000, 1000, 2
X, y = make_circles(n_samples=(n_train + n_test), factor=.3, noise=.05)
X_train, X_test = X[:n_train, :], X[n_train:, :]


# 2- Design the Experiment
# ------------------------
n_compo_to_try = 10
n_compo_range = [np.floor(np.exp((x / n_compo_to_try) * np.log(n_train)))
                 for x in range(0, n_compo_to_try + 1)]
n_iter = 3
arpack_all = False  # set this flag to True if you wish to run arpack for all


# 3- Benchmark
# ------------
# init
ref_time = np.empty((len(n_compo_range), n_iter)) * np.nan
a_time = np.empty((len(n_compo_range), n_iter)) * np.nan
r_time = np.empty((len(n_compo_range), n_iter)) * np.nan
# loop
for j, n_components in enumerate(n_compo_range):

    n_components = int(n_components)
    print("Performing kPCA with n_components = %i" % n_components)

    # A- reference (dense)
    for i in range(n_iter):
        start_time = datetime.now()
        ref_pred = KernelPCA(n_components, eigen_solver="dense") \
            .fit(X_train).transform(X_test)
        ref_time[j, i] = (datetime.now() - start_time).total_seconds()

    # B- arpack
    if arpack_all or n_components < 100:
        for i in range(n_iter):
            start_time = datetime.now()
            a_pred = KernelPCA(n_components, eigen_solver="arpack") \
                .fit(X_train).transform(X_test)
            # check that the result is still correct despite the approx
            assert_array_almost_equal(np.abs(a_pred), np.abs(ref_pred))
            a_time[j, i] = (datetime.now() - start_time).total_seconds()

    # C- randomized
    for i in range(n_iter):
        start_time = datetime.now()
        r_pred = KernelPCA(n_components, eigen_solver="randomized") \
            .fit(X_train).transform(X_test)
        # check that the result is still correct despite the approximation
        assert_array_almost_equal(np.abs(r_pred), np.abs(ref_pred))
        r_time[j, i] = (datetime.now() - start_time).total_seconds()

# Compute statistics for the 3 methods
avg_ref_time = ref_time.mean(axis=1)
std_ref_time = ref_time.std(axis=1)
avg_a_time = a_time.mean(axis=1)
std_a_time = a_time.std(axis=1)
avg_r_time = r_time.mean(axis=1)
std_r_time = r_time.std(axis=1)


# 4- Plots
# --------
plt.figure(figsize=(15, 20))

# Display 1 plot with error bars per method
plt.errorbar(n_compo_range, avg_ref_time, yerr=std_ref_time,
             marker='x', linestyle='', color='r', label='full')
plt.errorbar(n_compo_range, avg_a_time, yerr=std_a_time, marker='x',
             linestyle='', color='g', label='arpack')
plt.errorbar(n_compo_range, avg_r_time, yerr=std_r_time, marker='x',
             linestyle='', color='b', label='randomized')
plt.legend(loc='upper left')

# customize axes
ax = plt.gca()
ax.set_xscale('log')
ax.set_xlim(0, max(n_compo_range) * 1.1)
ax.set_ylabel("Execution time (s)")
ax.set_xlabel("n_components")

plt.title("Execution time comparison of kPCA on %i samples with %i "
          "features, according to the choice of `eigen_solver`"
          "" % (n_train, n_features))

plt.show()
