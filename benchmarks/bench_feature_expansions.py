import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from sklearn.preprocessing import PolynomialFeatures
from sys import argv
from time import time

degree = 2
trials = 5
num_rows = 100
dimensionalities = [50, 100, 200, 400, 800]
densities = [0.1, 0.4, 0.7, 1.0]
colors = ['r', 'b', 'g', 'k']
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
csc_linestyle = (0, (5, 10)) # loosely dashed
csr_linestyle = (0, (3, 1, 1, 1, 1, 1)) # densely dashdotdotted
dense_linestyle = (0, ()) # solid
fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig.add_subplot(111, frameon=False)
for color, density in zip(colors, densities):
    ax0.plot(dimensionalities, csr_times[density] / trials,
             label='csr, density=%s' % (density,),
             linestyle=csr_linestyle)
    ax0.plot(dimensionalities, dense_times[density] / trials,
             label='dense, density=%s' % (density,),
             linestyle=dense_linestyle)
    ax0.legend()
for color, density in zip(colors, densities):
    ax1.plot(dimensionalities, csr_times[density] / trials,
             label='csr, density=%s' % (density,),
             linestyle=csr_linestyle)
    ax1.plot(dimensionalities, csc_times[density] / trials,
             label='csc, density=%s' % (density,),
             linestyle=csc_linestyle)
    ax1.legend()
plt.xlabel('Dimensionality')
plt.ylabel('Average time (seconds) over %s trials' % (trials,), labelpad=20)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off',
                right='off')
plt.title('Time to compute degree=%s polynomial features of a %s row matrix'
          % (degree, num_rows,))
plt.show()
