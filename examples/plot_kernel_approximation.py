"""
==================================================
Explicit feature map approximation for RBF kernels
==================================================

An example shows how to use :class:`RBFSampler` to appoximate the feature map of an RBF
kernel for classification with an SVM on the digits dataset.
Results using a linear SVM in the original space, a linear SVM using the
approximate mapping and using a kernelized SVM are compared.
Timings and accuracy for varying amounts of Monte Carlo samplings for the
approximate mapping are shown.

Sampling more dimensions clearly leads to better classification results, but
comes at a greater cost. This means there is a tradeoff between runtime and
accuracy, given by the parameter n_components.  Note that solving the Linear
SVM and also the approximate kernel SVM could be greatly accelerated by using
stochastic gradient descent via :class:`SGDClassifier`. This is not easily possible for
the case of the kernelized SVM.

The usage of :class:`RBFSampler` is described in detail in :ref:`kernel_approximation`.

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

# The digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.data)
data = digits.data / 16.

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
accuracy = pl.gca()
print(kernel_svm_time)
print(linear_svm_time)
print(approx_kernel_times)
print(kernel_svm_score)
# second y axis for timeings
timescale = accuracy.twinx()

accuracy.plot(sample_sizes, approx_kernel_scores, label="approximate kernel")
timescale.plot(sample_sizes, approx_kernel_times, '--',
        label='runtime approx. kernel')

# horizontal lines for exact rbf and linear kernels:
accuracy.plot([sample_sizes[0], sample_sizes[-1]], [linear_svm_score,
    linear_svm_score], label="linear svm")
timescale.plot([sample_sizes[0], sample_sizes[-1]], [linear_svm_time,
        linear_svm_time], '--', label='runtime linear')

accuracy.plot([sample_sizes[0], sample_sizes[-1]], [kernel_svm_score,
    kernel_svm_score], label="rbf svm")
timescale.plot([sample_sizes[0], sample_sizes[-1]], [kernel_svm_time,
        kernel_svm_time], '--', label='runtime exact kernel')

# vertical line for dataset dimensionality = 64
accuracy.plot([64, 64], [0.7, 1], label="original dimensionality")

# legends and labels
accuracy.set_xlim(sample_sizes[0], sample_sizes[-1])
accuracy.set_ylim(np.min(approx_kernel_scores), 1)
accuracy.set_xlabel("Sampling steps = transformed feature dimension")
accuracy.set_ylabel("Classification accuracy")
timescale.set_ylabel("Training time in seconds (dashed lines)")
accuracy.legend(loc='best')
pl.show()
