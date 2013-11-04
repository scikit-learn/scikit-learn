import numpy as np
from scipy.sparse.base import issparse
from scipy.sparse.csr import csr_matrix
from sklearn.linear_model.stochastic_gradient import SGDRegressor
from sklearn.metrics import r2_score

np.random.seed(42)

n_samples, n_features = 1000, 200
X = np.random.randn(n_samples, n_features)
coef = 3 * np.random.randn(n_features)
inds = np.arange(n_features)
np.random.shuffle(inds)
coef[inds[10:]] = 0  # sparsify coef
y = np.dot(X, coef)

# add noise
y += 0.01 * np.random.normal((n_samples,))

# Split data in train set and test set
n_samples = X.shape[0]
X_train, y_train = X[:n_samples / 2], y[:n_samples / 2]
X_test, y_test = X[n_samples / 2:], y[n_samples / 2:]

print("X.shape:", X.shape, "y.shape:", y.shape)

###############################################################################
clf = SGDRegressor(penalty='l1', alpha=.2, fit_intercept=False, n_iter=2000)

clf.fit(X_train, y_train)


@profile
def benchmark_dense_predict():
    print "coeffs sparse ?", issparse(clf.coef_)
    for _ in range(1000):
        clf.predict(X_test)

@profile
def benchmark_sparse_predict():
    print "coeffs sparse ?", issparse(clf.coef_)
    for _ in range(1000):
        clf.predict(X_test)

def score(y_test, y_pred, case):
    r2 = r2_score(y_test, y_pred)
    print("r^2 on test data (%s) : %f" % (case, r2))

score(y_test, clf.predict(X_test), 'dense model')
#benchmark_dense_predict()
clf.sparsify()
score(y_test, clf.predict(X_test), 'sparse model')
#benchmark_sparse_predict()
