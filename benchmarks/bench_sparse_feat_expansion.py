#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from sklearn.preprocessing import PolynomialFeatures
from time import time

degree = 2
trials = 3
num_rows = 1000
dimensionalities = np.array([2 ** i for i in range(1, 11)])
degrees = [2, 3]
csr_times = {d: np.zeros(len(dimensionalities)) for d in degrees}
density = 0.01

for trial in range(trials):
    for degree in degrees:
        transform = PolynomialFeatures(
            degree=degree, include_bias=False, interaction_only=False
        )

        for dim_index, dim in enumerate(dimensionalities):
            print(trial, density, dim)
            X_csr = sparse.random(num_rows, dim, density).tocsr()
            X_dense = X_csr.toarray()
            # CSR
            t0 = time()
            transform.fit_transform(X_csr)
            csr_times[degree][dim_index] += time() - t0

csr_linestyle = (0, (3, 1, 1, 1, 1, 1))  # densely dashdotdotted
dense_linestyle = (0, ())  # solid

fig, axes = plt.subplots(nrows=len(degrees), ncols=1, figsize=(8, 10))
for degree, ax in zip(degrees, axes):
    ax.plot(
        dimensionalities,
        csr_times[degree] / trials,
        label="csr",
        linestyle=csr_linestyle,
    )
    ax.set_title(
        "density %0.2f, degree=%d, n_samples=%d" % (density, degree, num_rows)
    )
    ax.legend()
    ax.set_xlabel("Dimensionality")
    ax.set_ylabel("Time (seconds)")

plt.tight_layout()
plt.show()
