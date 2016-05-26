"""
============================================
BoxCox transformation on Boston housing data
============================================

This example shows the effect of boxcox transform on the boston
housing dataset.

"""

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_boston
from sklearn.preprocessing import BoxCoxTransformer
rng = np.random.RandomState(0)

dataset = load_boston()
X = dataset.data
y = dataset.target
n_samples = X.shape[0]

X = X[:, np.all(X > 0, axis=0)]  # restrict X to positive features
X_train = X[:n_samples // 2]
X_test = X[n_samples // 2:]

bct = BoxCoxTransformer(feature_indices=None)
bct.fit(X_train)
X_tr = bct.transform(X_test)

num_bins = 50

for i in [0, 10]:
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.hist(X_test[:, i], num_bins, normed=1, facecolor='green', alpha=0.5)
    ax1.set_title('Probplot before Box-Cox transformation %d' % i)
    ax1.set_xlabel('')
    ax1.set_ylabel('Probability')

    ax2.hist(X_tr[:, i], num_bins, normed=1, facecolor='green', alpha=0.5)
    ax2.set_title('Probplot after Box-Cox transformation %d' % i)
    ax2.set_xlabel('')
    ax2.set_ylabel('Probability')

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

plt.show()
