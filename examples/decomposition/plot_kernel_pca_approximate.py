"""
========================
Kernel PCA Approximation
========================

This example shows that using a limited number of components in Kernel PCA can
lead to a good approximation while being much faster to execute for large
datasets.

Description:
------------
2000 Random data samples are generated similarly to the "Kernel PCA" example
(the circles). A KernelPCA model is fit on that data with a reduced number of
principal components, varying from 4 to 1999. For each value, data is
transformed and inverse-transformed.

Original data is displayed in top left corner. Then subsequent plots in the
top row show approximations of this data made by KernelPCA transform + inverse
transform, for various number of components. The bottom row shows the data
samples projected along the first two principal components for each case.

What you can observe:
---------------------
With quite small `n_components`, the approximation is almost identical to the
full KernelPCA using 2000 components, while execution is much faster
(execution times are included in the plot titles)

Going further:
--------------
You can have a look at this `Kernel PCA Solvers comparison benchmark
<https://github.com/scikit-learn/scikit-learn/blob/master/benchmarks/bench_kernel_pca_solvers_time_vs_n_samples.py>`_
to compare execution times in more details.
"""
print(__doc__)  # noqa

# Authors: Sylvain MARIE, Schneider Electric
# License: BSD 3 clause

from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import KernelPCA
from sklearn.datasets import make_circles


# 1- Generate random data
# -----------------------
np.random.seed(0)
n_samples = 2000
X, y = make_circles(n_samples=n_samples, factor=.3, noise=.05)


# 2- Design experiment / Init plots
# ---------------------------------
plt.figure(figsize=(11, 5.5))

# Set this to true to create a logarithmic grid
use_log_grid = False
if use_log_grid:
    grid_size = 4
    n_components_range = [np.round(np.exp((x / grid_size) * np.log(n_samples)))
                          for x in range(1, grid_size + 1)]
else:
    # hardcoded for the example gallery
    n_components_range = [4, 40, 2000]

# top left: original
nb_cols = 1 + len(n_components_range)
plt.subplot(2, nb_cols, 1, aspect='equal')
plt.title("Original space")
reds = y == 0
blues = y == 1

plt.scatter(X[reds, 0], X[reds, 1], c="red",
            s=20, edgecolor='k')
plt.scatter(X[blues, 0], X[blues, 1], c="blue",
            s=20, edgecolor='k')
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")

# 2- Run experiments
# -----------------------
for i, n_compo in enumerate(reversed(n_components_range)):
    n_compo = int(n_compo)
    print("Fitting kPCA for n_components=%i..." % n_compo)
    start_time = datetime.now()

    # fit and transform
    kpca = KernelPCA(n_components=n_compo, kernel="rbf",
                     fit_inverse_transform=True, gamma=10)
    X_kpca = kpca.fit_transform(X)
    X_back = kpca.inverse_transform(X_kpca)
    elapsed = (datetime.now() - start_time).total_seconds()
    print("...done. Elapsed time: %.2fs" % elapsed)

    # original space after inverse transform
    plt.subplot(2, nb_cols, i + 2, aspect='equal')
    plt.scatter(X_back[reds, 0], X_back[reds, 1], c="red",
                s=20, edgecolor='k')
    plt.scatter(X_back[blues, 0], X_back[blues, 1], c="blue",
                s=20, edgecolor='k')
    plt.title("kPCA approx (n=%i)(%.2fs)" % (n_compo, elapsed))
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")

    # plot the first two principal components
    plt.subplot(2, nb_cols, i + 2 + nb_cols, aspect='equal')
    plt.scatter(X_kpca[reds, 0], X_kpca[reds, 1], c="red",
                s=20, edgecolor='k')
    plt.scatter(X_kpca[blues, 0], X_kpca[blues, 1], c="blue",
                s=20, edgecolor='k')
    plt.title("First 2PC (n=%i)" % n_compo)
    plt.xlabel(r"1st PC in space induced by $\phi$")
    plt.ylabel("2nd component")

plt.tight_layout()
plt.show()
