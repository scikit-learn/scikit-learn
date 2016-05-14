"""
============================================
BoxCox transformation on Boston housing data
============================================

This example shows the effect of boxcox transform on the boston
housing dataset.
  
"""

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from sklearn.datasets import load_boston
from sklearn.preprocessing import BoxCoxTransformer
rng = np.random.RandomState(0)

dataset = load_boston()
X = dataset.data
y = dataset.target
n_samples = X.shape[0]

features_new = np.where(np.all(X>0, axis=0))[0]
X_fit = X[:n_samples/2, features_new]
X_transform = X[n_samples/2:, features_new]
BCTr = BoxCoxTransformer(feature_indices=None)
BCTr.fit(X_fit)
X_tr = BCTr.transform(X_transform)

num_bins = 50

for i in range(11):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    x = X_transform[:, i]
    plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
    ax1.set_title('Probplot before Box-Cox transformation '+str(i))

    plt.xlabel('')
    plt.ylabel('Probability')

    ax2 = fig.add_subplot(212)
    xt = X_tr[:, i]
    plt.hist(xt, num_bins, normed=1, facecolor='green', alpha=0.5)
    ax2.set_title('Probplot after Box-Cox transformation '+str(i))

    plt.xlabel('')
    plt.ylabel('Probability')

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    plt.show()
