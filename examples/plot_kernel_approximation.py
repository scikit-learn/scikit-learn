"""
==================================================
Explicit feature map approximation for RBF kernels
==================================================

An example shows how to use RBFSampler to appoximate the feature map of an RBF
kernel for classification with an SVM on the digits dataset.
Results using a linear SVM in the original space, a linear SVM using the
approximate mapping and using a kernelized SVM are compared.
Timings and accuracy for varying amounts of Monte Carlo samplings for the
approximate mapping are shown.

Sampling more dimensions clearly leads to better classification results, but
comes at a greater cost. This means there is a tradeoff between runtime and
accuracy, given by the parameter n_components.
"""
print __doc__

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
#         modified Andreas Mueller
# License: Simplified BSD

# Standard scientific Python imports
import pylab as pl
import numpy as np

# Import datasets, classifiers and performance metrics
from sklearn import datasets, svm, pipeline
from sklearn.kernel_approximation import RBFSampler

# The digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
data = digits.images.reshape((n_samples, -1))

# We learn the digits on the first half of the digits
data_train, targets_train = data[:n_samples / 2], digits.target[:n_samples / 2]

# Now predict the value of the digit on the second half:
data_test, targets_test = data[n_samples / 2:], digits.target[n_samples / 2:]

# Create a classifier: a support vector classifier
kernel_svm = svm.SVC(gamma=0.001)
linear_svm = svm.LinearSVC()

# create pipeline from kernel approximation
# and linear svm
feature_map = RBFSampler(gamma=0.001)
approx_kernel_svm = pipeline.Pipeline([("feature_map", feature_map),
    ("svm", svm.LinearSVC())])

# fit and predict using linear and kernel svm:
kernel_svm.fit(data_train, targets_train)
kernel_svm_score = kernel_svm.score(data_test, targets_test)

linear_svm.fit(data_train, targets_train)
linear_svm_score = linear_svm.score(data_test, targets_test)

sample_sizes = 20 * np.arange(1, 15)
approx_kernel_scores = []
for D in sample_sizes:
    approx_kernel_svm.set_params(feature_map__n_components=D)
    approx_kernel_svm.fit(data_train, targets_train)
    score = approx_kernel_svm.score(data_test, targets_test)
    approx_kernel_scores.append(score)

# plot the results:
pl.plot(sample_sizes, approx_kernel_scores, label="approximate kernel")

# horizontal lines for exact rbf and linear kernels:
pl.plot([sample_sizes[0], sample_sizes[-1]], [linear_svm_score,
    linear_svm_score], label="linear svm")
pl.plot([sample_sizes[0], sample_sizes[-1]], [kernel_svm_score,
    kernel_svm_score], label="rbf svm")

# vertical line for dataset dimensionality = 64
pl.plot([64, 64], [0.7, 1], label="original dimensionality")

# legends and labels
pl.xlim(sample_sizes[0], sample_sizes[-1])
pl.ylim(np.min(approx_kernel_scores), 1)
pl.xlabel("Sampling steps = transformed feature dimension")
pl.ylabel("Classification accuracy")
pl.legend(loc='best')
pl.show()
