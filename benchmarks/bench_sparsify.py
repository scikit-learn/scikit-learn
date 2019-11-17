"""
Benchmark SGD prediction time with dense/sparse coefficients.

Invoke with
-----------

$ kernprof.py -l sparsity_benchmark.py
$ python -m line_profiler sparsity_benchmark.py.lprof

Typical output
--------------

input data sparsity: 0.050000
true coef sparsity: 0.000100
test data sparsity: 0.027400
model sparsity: 0.000024
r^2 on test data (dense model) : 0.233651
r^2 on test data (sparse model) : 0.233651
Wrote profile results to sparsity_benchmark.py.lprof
Timer unit: 1e-06 s

File: sparsity_benchmark.py
Function: benchmark_dense_predict at line 51
Total time: 0.532979 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    51                                           @profile
    52                                           def benchmark_dense_predict():
    53       301          640      2.1      0.1      for _ in range(300):
    54       300       532339   1774.5     99.9          clf.predict(X_test)

File: sparsity_benchmark.py
Function: benchmark_sparse_predict at line 56
Total time: 0.39274 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    56                                           @profile
    57                                           def benchmark_sparse_predict():
    58         1        10854  10854.0      2.8      X_test_sparse = csr_matrix(X_test)
    59       301          477      1.6      0.1      for _ in range(300):
    60       300       381409   1271.4     97.1          clf.predict(X_test_sparse)
"""

from scipy.sparse.csr import csr_matrix
import numpy as np
from sklearn.linear_model import SGDRegressor
from sklearn.metrics import r2_score

np.random.seed(42)


def sparsity_ratio(X):
    return np.count_nonzero(X) / float(n_samples * n_features)

n_samples, n_features = 5000, 300
X = np.random.randn(n_samples, n_features)
inds = np.arange(n_samples)
np.random.shuffle(inds)
X[inds[int(n_features / 1.2):]] = 0  # sparsify input
print("input data sparsity: %f" % sparsity_ratio(X))
coef = 3 * np.random.randn(n_features)
inds = np.arange(n_features)
np.random.shuffle(inds)
coef[inds[n_features // 2:]] = 0  # sparsify coef
print("true coef sparsity: %f" % sparsity_ratio(coef))
y = np.dot(X, coef)

# add noise
y += 0.01 * np.random.normal((n_samples,))

# Split data in train set and test set
n_samples = X.shape[0]
X_train, y_train = X[:n_samples // 2], y[:n_samples // 2]
X_test, y_test = X[n_samples // 2:], y[n_samples // 2:]
print("test data sparsity: %f" % sparsity_ratio(X_test))

###############################################################################
clf = SGDRegressor(penalty='l1', alpha=.2, max_iter=2000,
                   tol=None)
clf.fit(X_train, y_train)
print("model sparsity: %f" % sparsity_ratio(clf.coef_))


def benchmark_dense_predict():
    for _ in range(300):
        clf.predict(X_test)


def benchmark_sparse_predict():
    X_test_sparse = csr_matrix(X_test)
    for _ in range(300):
        clf.predict(X_test_sparse)


def score(y_test, y_pred, case):
    r2 = r2_score(y_test, y_pred)
    print("r^2 on test data (%s) : %f" % (case, r2))

score(y_test, clf.predict(X_test), 'dense model')
benchmark_dense_predict()
clf.sparsify()
score(y_test, clf.predict(X_test), 'sparse model')
benchmark_sparse_predict()
