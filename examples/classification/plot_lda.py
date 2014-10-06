"""
====================================================================
Normal and Shrinkage Linear Discriminant Analysis for classification
====================================================================

Shows how shrinkage improves classification.
"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_blobs
from sklearn.lda import LDA


n_train = 20  # samples for training
n_test = 200  # samples for testing
n_averages = 50  # how often to repeat classification
n_features_max = 75  # maximum number of features


def generate_data(n_samples, n_features):
    """Generate `n_samples` samples of data with 1 discriminative features
    and n_`features` non-discriminative features."""
    X, y = make_blobs(n_samples=n_samples, n_features=1, centers=[[-2], [2]])

    # add non-discriminative features
    if n_features > 1:
        X = np.hstack([X, np.random.randn(n_samples, n_features - 1)])
    return X, y

acc_clf1, acc_clf2 = [], []
m_range = range(1, n_features_max + 1)
for m in m_range:
    score_clf1, score_clf2 = 0, 0
    for i in range(n_averages):
        X, y = generate_data(n_train, m)

        clf1 = LDA(solver='lsqr', alpha='ledoit_wolf').fit(X, y)
        clf2 = LDA(solver='lsqr', alpha=None).fit(X, y)

        X, y = generate_data(n_test, m)
        score_clf1 += clf1.score(X, y)
        score_clf2 += clf2.score(X, y)

    acc_clf1.append(score_clf1 / n_averages)
    acc_clf2.append(score_clf2 / n_averages)

m_range = np.array(m_range) / n_train

plt.plot(m_range, acc_clf1, linewidth=2, label="LDA with shrinkage", color='r')
plt.plot(m_range, acc_clf2, linewidth=2, label="LDA", color='g')

plt.xlabel('n_features / n_samples')
plt.ylabel('Classification accuracy')

plt.legend(loc=1, prop={'size': 8})
plt.suptitle('LDA vs. shrinkage LDA (1 discriminative feature)')
plt.show()
