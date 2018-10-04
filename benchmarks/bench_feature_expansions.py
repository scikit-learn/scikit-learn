import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from sklearn.preprocessing import PolynomialFeatures
from time import time

degree = 2
trials = 3
num_rows = 1000
dimensionalities = np.array([6, 12, 25, 50])
densities = np.array([0.01, 0.1, 1.0])
csr_times = {d: np.zeros(len(dimensionalities)) for d in densities}
csc_times = {d: np.zeros(len(dimensionalities)) for d in densities}
dense_times = {d: np.zeros(len(dimensionalities)) for d in densities}
transform = PolynomialFeatures(degree=degree, include_bias=False,
                               interaction_only=False)

for trial in range(trials):
    for density in densities:
        for dim_index, dim in enumerate(dimensionalities):
            print(trial, density, dim)
            X_csr = sparse.random(num_rows, dim, density).tocsr()
            X_csc = X_csr.tocsc()
            X_dense = X_csr.toarray()
            # CSR
            t0 = time()
            transform.fit_transform(X_csr)
            csr_times[density][dim_index] += time() - t0
            # CSC
            t0 = time()
            transform.fit_transform(X_csc)
            csc_times[density][dim_index] += time() - t0
            # Dense
            t0 = time()
            transform.fit_transform(X_dense)
            dense_times[density][dim_index] += time() - t0

csc_linestyle = (0, (5, 10))  # loosely dashed
csr_linestyle = (0, (3, 1, 1, 1, 1, 1))  # densely dashdotdotted
dense_linestyle = (0, ())  # solid

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 10))
for density, ax in zip(densities, axes):

    ax.loglog(dimensionalities, csr_times[density] / trials,
              label='csr', linestyle=csr_linestyle)
    ax.loglog(dimensionalities, dense_times[density] / trials,
              label='dense', linestyle=dense_linestyle)
    ax.loglog(dimensionalities, csc_times[density] / trials,
              label='csc', linestyle=csc_linestyle)
    ax.set_title("density %0.2f, degree=%d, n_samples=%d" %
                 (density, degree, num_rows))
    ax.legend()
    ax.set_xlabel('Dimensionality')
    ax.set_ylabel('Time (seconds)')

plt.tight_layout()
plt.show()
