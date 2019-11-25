"""
====================================================================
Ledoit Wolf and OAS Linear Discriminant Analysis for classification
====================================================================

Shows how shrinkage improves classification.
"""
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import toeplitz, cholesky
from sklearn.covariance import OAS
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

np.random.seed(0)

n_features = 100
# simulation covariance matrix (AR(1) process)
r = 0.1
real_cov = toeplitz(r**np.arange(n_features))
coloring_matrix = cholesky(real_cov)
loc = 0.3
n_repeat = 100

n_samples_range = np.arange(6, 32, 2)
lw_accuracy = np.zeros((n_samples_range.size, n_repeat))
oa_accuracy = np.zeros((n_samples_range.size, n_repeat))
for i, n_samples in enumerate(n_samples_range):
    n_train = n_test = int(n_samples / 2)
    for j in range(n_repeat):
        X1 = np.dot(np.random.normal(size=(n_samples, n_features)),
                    coloring_matrix.T)
        X2 = X1 + loc
        X_train = np.concatenate([X1[:n_train], X2[:n_train]], axis=0)
        Y_train = np.array([0] * n_train + [1] * n_train)
        X_test = np.concatenate([X1[n_train:], X2[n_train:]], axis=0)
        Y_test = np.array([0] * n_train + [1] * n_train)

        oa = OAS(store_precision=False, assume_centered=False)
        lda_oa = LinearDiscriminantAnalysis(solver="lsqr",
                                            covariance_estimator=oa)
        lda_lw = LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")

        lda_oa.fit(X_train, Y_train)
        oa_accuracy[i, j] = np.mean(lda_oa.predict(X_test) == Y_test)

        lda_lw.fit(X_train, Y_train)
        lw_accuracy[i, j] = np.mean(lda_lw.predict(X_test) == Y_test)

# plot MSE
plt.figure()
plt.errorbar(n_samples_range,
             lw_accuracy.mean(1),
             yerr=lw_accuracy.std(1),
             label='Ledoit-Wolf',
             color='navy')
plt.errorbar(n_samples_range,
             oa_accuracy.mean(1),
             yerr=oa_accuracy.std(1),
             label='OAS',
             color='darkorange')
plt.ylabel("Test accuracy")
plt.legend(loc="lower right")
plt.title("Discriminant analysis with different covariance estimators")
plt.xlim(5, 31)
plt.show()
