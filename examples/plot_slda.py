"""
====================================================================
Normal and Shrinkage Linear Discriminant Analysis for classification
====================================================================

Shows how shrinkage improves classification.
"""

import numpy as np
import matplotlib.pyplot as plt

from sklearn.lda import LDA
from sklearn.slda import LDA as sLDA

n_train = 10     # samples per class for training
n_test = 100    # samples per class for training
n_averages = 50


def generate_data(n_samples, n_features):
    """ generate `n_samples` samples of data with 1 discriminative features and n_`features` non-discriminative features."""
    X = np.vstack([np.random.randn(n_samples, 1) - 2, np.random.randn(n_samples, 1) + 2])
    y = np.hstack([[-1]*n_samples, [1]*n_samples])
    # add non discriminative features
    X = np.hstack([X, np.random.randn(2*n_samples, n_features)])
    return X, y


acc_lda, acc_slda, acc_nlda = [], [], []
m_range = range(50)
for m in m_range:
    print(m)
    tmp_lda, tmp_slda, tmp_nlda = 0, 0, 0
    for i in range(n_averages):
        X, y = generate_data(n_train, m)

        lda = LDA().fit(X, y)
        slda = sLDA().fit(X, y)
        nlda = sLDA(shrinkage=None).fit(X, y)

        X, y = generate_data(n_test, m)
        tmp_lda += np.mean(lda.predict(X) == y) / n_averages
        tmp_slda += np.mean(slda.predict(X) == y) / n_averages
        tmp_nlda += np.mean(nlda.predict(X) == y) / n_averages

    acc_lda.append(tmp_lda)
    acc_slda.append(tmp_slda)
    acc_nlda.append(tmp_nlda)

m_range = (1+np.array(m_range))/(2*n_train)

plt.plot(m_range, acc_lda, linewidth=2, label='sklearn.lda')
plt.plot(m_range, acc_nlda, linewidth=2, label='empirical cov')
plt.plot(m_range, acc_slda, linewidth=2, label='Ledoit-Wolf')
#plt.semilogx()

plt.xlabel('n_features / n_samples')
plt.ylabel('Classification accuracy')

plt.legend(loc=3)

plt.suptitle('LDA vs sLDA (1 discriminative feature)')

plt.show()