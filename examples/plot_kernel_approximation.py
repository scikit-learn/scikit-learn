"""
==================================================
Explicit feature map approximation for RBF kernels
==================================================

.. currentmodule:: sklearn.kernel_approximation

An example shows how to use :class:`RBFSampler` to appoximate the feature map
of an RBF kernel for classification with an SVM on the digits dataset.  Results
using a linear SVM in the original space, a linear SVM using the approximate
mapping and using a kernelized SVM are compared.  Timings and accuracy for
varying amounts of Monte Carlo samplings for the approximate mapping are shown.

Sampling more dimensions clearly leads to better classification results, but
comes at a greater cost. This means there is a tradeoff between runtime and
accuracy, given by the parameter n_components.  Note that solving the Linear
SVM and also the approximate kernel SVM could be greatly accelerated by using
stochastic gradient descent via :class:`sklearn.linear_model.SGDClassifier`.
This is not easily possible for the case of the kernelized SVM.

The second plot visualized the decision surfaces of the RBF kernel SVM and
the linear SVM with approximate kernel map.
The plot shows decision surfaces of the classifiers projected onto
the first two principal components of the data. This visualization should
be taken with a grain of salt since it is just an interesting slice through
the decision surface in 64 dimensions. In particular note that
a datapoint (represented as a dot) does not necessarily be classified
into the region it is lying in, since it will not lie on the plane
that the first two principal components span.

The usage of :class:`RBFSampler` is described in detail in
:ref:`kernel_approximation`.

"""
print __doc__

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
#         modified Andreas Mueller
# License: Simplified BSD

# Standard scientific Python imports
import pylab as pl
import numpy as np
from time import time

# Import datasets, classifiers and performance metrics
from sklearn import datasets, svm, pipeline
from sklearn.kernel_approximation import RBFSampler
from sklearn.decomposition import PCA

# The digits dataset
digits = datasets.load_digits(n_class=9)

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.data)
data = digits.data / 16.
data -= data.mean(axis=0)

# We learn the digits on the first half of the digits
data_train, targets_train = data[:n_samples / 2], digits.target[:n_samples / 2]


# Now predict the value of the digit on the second half:
data_test, targets_test = data[n_samples / 2:], digits.target[n_samples / 2:]
#data_test = scaler.transform(data_test)

# Create a classifier: a support vector classifier
kernel_svm = svm.SVC(gamma=.2)
linear_svm = svm.LinearSVC()

# create pipeline from kernel approximation
# and linear svm
feature_map = RBFSampler(gamma=.2, random_state=1)
approx_kernel_svm = pipeline.Pipeline([("feature_map", feature_map),
    ("svm", svm.LinearSVC())])

# fit and predict using linear and kernel svm:

kernel_svm_time = time()
kernel_svm.fit(data_train, targets_train)
kernel_svm_score = kernel_svm.score(data_test, targets_test)
kernel_svm_time = time() - kernel_svm_time

linear_svm_time = time()
linear_svm.fit(data_train, targets_train)
linear_svm_score = linear_svm.score(data_test, targets_test)
linear_svm_time = time() - linear_svm_time

sample_sizes = 50 * np.arange(1, 10)
approx_kernel_scores = []
approx_kernel_times = []

for D in sample_sizes:
    approx_kernel_svm.set_params(feature_map__n_components=D)
    approx_kernel_timing = time()
    approx_kernel_svm.fit(data_train, targets_train)
    approx_kernel_times.append(time() - approx_kernel_timing)
    score = approx_kernel_svm.score(data_test, targets_test)
    approx_kernel_scores.append(score)

# plot the results:
accuracy = pl.subplot(211)
# second y axis for timeings
timescale = pl.subplot(212)

accuracy.plot(sample_sizes, approx_kernel_scores, label="approx. kernel")
timescale.plot(sample_sizes, approx_kernel_times, '--',
        label='approx. kernel')

# horizontal lines for exact rbf and linear kernels:
accuracy.plot([sample_sizes[0], sample_sizes[-1]], [linear_svm_score,
    linear_svm_score], label="linear svm")
timescale.plot([sample_sizes[0], sample_sizes[-1]], [linear_svm_time,
        linear_svm_time], '--', label='linear svm')

accuracy.plot([sample_sizes[0], sample_sizes[-1]], [kernel_svm_score,
    kernel_svm_score], label="rbf svm")
timescale.plot([sample_sizes[0], sample_sizes[-1]], [kernel_svm_time,
        kernel_svm_time], '--', label='rbf svm')

# vertical line for dataset dimensionality = 64
accuracy.plot([64, 64], [0.7, 1], label="n_features")

# legends and labels
accuracy.set_title("Classification accuracy")
timescale.set_title("Training times")
accuracy.set_xlim(sample_sizes[0], sample_sizes[-1])
accuracy.set_xticks(())
accuracy.set_ylim(np.min(approx_kernel_scores), 1)
timescale.set_xlabel("Sampling steps = transformed feature dimension")
accuracy.set_ylabel("Classification accuracy")
timescale.set_ylabel("Training time in seconds")
accuracy.legend(loc='best')
timescale.legend(loc='best')

# visualize the decision surface, projected down to the first
# two principal components of the dataset
pca = PCA(n_components=8).fit(data_train)

X = pca.transform(data_train)

# Gemerate grid along first two principal components
multiples = np.arange(-2, 2, 0.1)
# steps along first component
first = multiples[:, np.newaxis] * pca.components_[0, :]
# steps along second component
second = multiples[:, np.newaxis] * pca.components_[1, :]
# combine
grid = first[np.newaxis, :, :] + second[:, np.newaxis, :]
flat_grid = grid.reshape(-1, data.shape[1])

# title for the plots
titles = ['SVC with rbf kernel',
          'SVC (linear kernel) with rbf feature map\n n_components=100']

pl.figure(figsize=(12, 5))
pl.set_cmap(pl.cm.Paired)


# predict and plot
for i, clf in enumerate((kernel_svm, approx_kernel_svm)):
    # Plot the decision boundary. For that, we will asign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    pl.subplot(1, 2, i + 1)
    Z = clf.predict(flat_grid)

    # Put the result into a color plot
    Z = Z.reshape(grid.shape[:-1])
    pl.set_cmap(pl.cm.Paired)
    pl.contourf(multiples, multiples, Z)
    pl.axis('off')

    # Plot also the training points
    pl.scatter(X[:, 0], X[:, 1], c=targets_train)

    pl.title(titles[i])
pl.show()
