import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from sklearn.preprocessing import PolynomialFeatures
from time import time

degree = 2
trials = 1
num_rows = 100
dimensionalities = np.array([25, 50, 100, 200])
densities = np.array([0.01, 0.1, 1.0])
colors = ['#d7191c', '#abdda4', '#2b83ba']
assert(len(colors) == len(densities))
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

# The CSC algorithm takes so long that it drastically overshadows dense and
# CSR algorithm times, so don't include it in the dense / CSR plot. Instead,
# compare it with the CSR method in the bottom plot.
csc_linestyle = (0, (5, 10))  # loosely dashed
csr_linestyle = (0, (3, 1, 1, 1, 1, 1))  # densely dashdotdotted
dense_linestyle = (0, ())  # solid
for color, density in zip(colors, densities):
    plt.loglog(dimensionalities, csr_times[density] / trials,
               label='csr, density=%s' % (density,),
               linestyle=csr_linestyle, color=color, alpha=0.7)
    plt.loglog(dimensionalities, dense_times[density] / trials,
               label='dense, density=%s' % (density,),
               linestyle=dense_linestyle, color=color, alpha=0.7)
    plt.loglog(dimensionalities, csc_times[density] / trials,
               label='csc, density=%s' % (density,),
               linestyle=csc_linestyle, color=color, alpha=0.7)
plt.legend()
plt.xlabel('Dimensionality')
plt.ylabel('Average time (seconds) over %s trials' % (trials,), labelpad=20)
plt.title('Logscale time to compute degree=%s polynomial features of a %s row '
          'matrix' % (degree, num_rows,))
plt.show()
