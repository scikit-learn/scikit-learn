import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from sklearn.preprocessing import PolynomialFeatures
from sys import argv
from time import time

if len(argv) < 2:
    print("Please specify where to save the resulting plot "
          "via the a command line param.")
    exit()

degree = 2
trials = 5
num_rows = 100
densities = [0.1, 0.4, 0.7, 1.0]
dimensionalities = [100, 200, 400, 800, 1600]
dense_times = {d: np.zeros(len(dimensionalities)) for d in densities}
sparse_times = {d: np.zeros(len(dimensionalities)) for d in densities}
transform = PolynomialFeatures(degree=degree, include_bias=False,
                               interaction_only=False)

for trial in range(trials):
    for density in densities:
        for dim_index, dim in enumerate(dimensionalities):
            print(trial, density, dim)
            X_csr = sparse.random(num_rows, dim, density).tocsr()
            X_dense = X_csr.toarray()
            # Sparse
            t0 = time()
            transform.fit_transform(X_csr)
            sparse_times[density][dim_index] += time() - t0
            # Dense
            t0 = time()
            transform.fit_transform(X_dense)
            dense_times[density][dim_index] += time() - t0

colors = ['r', 'b', 'g', 'k']
assert(len(colors) == len(densities))
for color, density in zip(colors, densities):
    plt.plot(dimensionalities, sparse_times[density] / trials,
             label='csr, density=%s' % (density,),
             linestyle='--', color=color)
    plt.plot(dimensionalities, dense_times[density] / trials,
             label='dense, density=%s' % (density,),
             linestyle='-', color=color)
plt.xlabel('Dimensionality')
plt.ylabel('Average time (seconds) over %s trials' % (trials,))
plt.title('Time to compute degree=%s polynomial features of a %s row matrix'
          % (degree, num_rows,))
plt.legend()
plt.savefig(argv[1])
