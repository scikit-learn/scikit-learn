"""
=====================
BoxCox transformation
=====================
This example shows the effect of boxcox transform.
The transformation evaluates lambda to maximise the log-likelihood
and tranforms exponential distribution to approximate normal distribution

"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import BoxCoxTransformer
rng = np.random.RandomState(42)

n_samples = 3000
X = rng.randn(n_samples, 1)
X[:, 0] = np.exp(X[:, 0])

X_train = X[:n_samples // 2]
X_test = X[n_samples // 2:]

bct = BoxCoxTransformer(feature_indices=None)
bct.fit(X_train)
X_tr = bct.transform(X_test)

num_bins = 50

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.hist(X_test[:, 0], num_bins, normed=1, facecolor='green', alpha=0.5)
ax1.set_title('Probplot before Box-Cox transformation')
ax1.set_xlabel('')
ax1.set_ylabel('Probability')

ax2.hist(X_tr[:, 0], num_bins, normed=1, facecolor='green', alpha=0.5)
ax2.set_title('Probplot after Box-Cox transformation')
ax2.set_xlabel('')
ax2.set_ylabel('Probability')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)

plt.show()
