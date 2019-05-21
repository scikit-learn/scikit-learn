"""
=================================
Bernoulli Mixture Model Means
=================================

Plot the means of a mixture of Bernoulli distributions. This approximates
the analysis of handwritten digits as described in [Bishop2006]_.

.. [Bishop2006] Christopher M. Bishop. 2006.
Pattern Recognition and Machine Learning
(Information Science and Statistics).
Springer-Verlag, Berlin, Heidelberg.
"""

from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn import mixture
from sklearn.preprocessing import Binarizer

n_pixels = 8

digits = datasets.load_digits()

digit_indexes = defaultdict(list)
for digit_index, digit_value in enumerate(digits.target):
    digit_indexes[digit_value] += [digit_index]

# data comprises 0s and 4s
zeros = digits.data[digit_indexes[0]]
fours = digits.data[digit_indexes[4]]
X = np.concatenate((zeros, fours), axis=0)

# binarize the data
threshold = X.max() // 2
X = Binarizer(threshold).fit_transform(X)

# flatten the images to 1d
X = X.reshape((-1, n_pixels * n_pixels))

# Fit a Bernoulli mixture with EM using two components
bmm = mixture.BernoulliMixture(n_components=2).fit(X)

# check which prototype matches an example 0 best
cluster_id_test = np.mean(
    [bmm.predict(X[i].reshape(1, -1)) for i in range(10)]
    )

if cluster_id_test < 0.5:
    zero_mean_fitted = bmm.means_[0].reshape((n_pixels, n_pixels))
    four_mean_fitted = bmm.means_[1].reshape((n_pixels, n_pixels))
else:
    zero_mean_fitted = bmm.means_[1].reshape((n_pixels, n_pixels))
    four_mean_fitted = bmm.means_[0].reshape((n_pixels, n_pixels))

plotters = [np.mean(zeros, 0).reshape((n_pixels, n_pixels)),
            X[0].reshape((n_pixels, n_pixels)),
            zero_mean_fitted,
            np.mean(fours, 0).reshape((n_pixels, n_pixels)),
            X[-1].reshape((n_pixels, n_pixels)),
            four_mean_fitted]

titles = ["true-prototype", "example data", "fitted prototype"]
for plot_index, plot_data in enumerate(plotters):
    plt.subplot(2, 3, plot_index+1)
    plt.imshow(plot_data)
    plt.xticks(())
    plt.yticks(())
    if plot_index // 3 == 0:
        plt.title(titles[plot_index % 3])

plt.show()
