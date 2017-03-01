"""
==================================================
Explicit feature map approximation for RBF kernels
==================================================

An example illustrating the approximation of the feature map
of an RBF kernel.

.. currentmodule:: sklearn.kernel_approximation

It shows how to use :class:`RBFSampler` and :class:`Nystroem` to
approximate the feature map of an RBF kernel for classification with an SVM on
the digits dataset. Results using a linear SVM in the original space, a linear
SVM using the approximate mappings and using a kernelized SVM are compared.
Timings and accuracy for varying amounts of Monte Carlo samplings (in the case
of :class:`RBFSampler`, which uses random Fourier features) and different sized
subsets of the training set (for :class:`Nystroem`) for the approximate mapping
are shown.

Please note that the dataset here is not large enough to show the benefits
of kernel approximation, as the exact SVM is still reasonably fast.

Sampling more dimensions clearly leads to better classification results, but
comes at a greater cost. This means there is a tradeoff between runtime and
accuracy, given by the parameter n_components. Note that solving the Linear
SVM and also the approximate kernel SVM could be greatly accelerated by using
stochastic gradient descent via :class:`sklearn.linear_model.SGDClassifier`.
This is not easily possible for the case of the kernelized SVM.

The second plot visualized the decision surfaces of the RBF kernel SVM and
the linear SVM with approximate kernel maps.
The plot shows decision surfaces of the classifiers projected onto
the first two principal components of the data. This visualization should
be taken with a grain of salt since it is just an interesting slice through
the decision surface in 64 dimensions. In particular note that
a datapoint (represented as a dot) does not necessarily be classified
into the region it is lying in, since it will not lie on the plane
that the first two principal components span.

The usage of :class:`RBFSampler` and :class:`Nystroem` is described in detail
in :ref:`kernel_approximation`.

"""
print(__doc__)

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
#         Andreas Mueller <amueller@ais.uni-bonn.de>
# License: BSD 3 clause

# Standard scientific Python imports
import matplotlib.pyplot as plt
import numpy as np
from time import time

# Import datasets, classifiers and performance metrics
from sklearn import datasets, svm, pipeline
from sklearn.kernel_approximation import (RBFSampler,
                                          Nystroem)
from sklearn.decomposition import PCA
from collections import defaultdict


def timing(callable, *args, **kwargs):
    """Time the call of a function and return time and
    result. Passing args and kwargs to function.
    """
    init_time = time()
    result = callable(*args, **kwargs)
    return result, time() - init_time

def fit_score(clf, train, test):
    """Call fit and score on a classifier."""
    clf.fit(*train)
    return clf.score(*test)

def plot_svm_data(ax, x_values, y_vals_dict, **lineargs):
    """Plot several lines in a dict on an axis, using the key of the 
    dict as the label.

    ax: a matplotlib Axes object.
    x_values: iterable of values to plot along x-axis.
    y_values_dict: dict of label:data pairs to plot.
                   plotted as horizontal line if data is scalar.
    lineargs: kwargs to pass to plot call(s)
    """
    artists = []
    for label, y_values in y_vals_dict.items():
        if not hasattr(y_values, '__len__'):
            y_values = [y_values]*len(x_values)
        line = ax.plot(x_values, y_values, label=label, **lineargs)
        artists.append(line)
    return artists

def flat_grid_from_pca(pca, multiples):
    """Generate grid along first two principal components
    steps along first component

    pca: a fitted PCA object.
    multiples: a range of steps to make through components
    """
    first = multiples[:, np.newaxis] * pca.components_[0, :]
    # steps along second component
    second = multiples[:, np.newaxis] * pca.components_[1, :]
    # combine
    grid = first[np.newaxis, :, :] + second[:, np.newaxis, :]
    flat_grid = grid.reshape(-1, grid.shape[-1])
    return flat_grid

def plot_projected_decision_surface(ax, clf, X, multiples, flat_grid):
    """Plot the projected decision surface.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a fitted classifier
    X: an array whose first two columns should be axes of plot.
    multiples: 
    """

    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, x_max]x[y_min, y_max].
    Z = clf.predict(flat_grid)

    # Put the result into a color plot
    s = int(np.sqrt(flat_grid.shape[0]))
    Z = Z.reshape((s, s))
    contour = ax.contourf(multiples, multiples, Z, cmap=plt.cm.Paired)
    return contour


# The digits dataset
digits = datasets.load_digits(n_class=9)

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.data)
data = digits.data / 16.
data -= data.mean(axis=0)

# We learn the digits on the first half of the digits
split_index = int(n_samples / 2)
data_train, targets_train = data[:split_index], digits.target[:split_index]


# We will predict the value of the digit on the second half:
data_test, targets_test = data[split_index:], digits.target[split_index:]

# Create support vector classifier
svms = {'linear svm': svm.LinearSVC(),
        'rbf svm': svm.SVC(gamma=.2)}

# create pipeline from kernel approximation
# and linear svm
approximate_kernel_labels = ("Nystroem approx. kernel",
                             "Fourier approx. kernel")
approximate_kernels = dict(zip(approximate_kernel_labels,
                               (Nystroem(gamma=.2, random_state=1),
                                RBFSampler(gamma=.2, random_state=1))))
approx_svm = {label: pipeline.Pipeline([('feature_map', approx_kernel),
                                        ('svm', svm.LinearSVC())])
              for label, approx_kernel in approximate_kernels.items()}


# fit and predict using linear and kernel svm:
svm_times, svm_scores = {}, {}
for kernel, clf in svms.items():
    score, svm_performance = timing(fit_score, clf, (data_train, targets_train),
                          (data_test, targets_test))
    svm_scores[kernel] = score
    svm_times[kernel] = svm_performance


# create timing, accuracy data for approximate kernel models
sample_sizes = 30 * np.arange(1, 10)
scores, times = defaultdict(list), defaultdict(list)
for D in sample_sizes:
    for feature_map, clf in approx_svm.items():
        clf.set_params(feature_map__n_components=D)
        fitted_clf, _time = timing(clf.fit, data_train, targets_train)
        times[feature_map].append(_time)
        score = fitted_clf.score(data_test, targets_test)
        scores[feature_map].append(score)

# figure layout
fig, subplots = plt.subplots(2, 1, figsize=(8, 8))
# zip together information for creating plots.
plot_layout = zip(subplots,
                  (scores, times),
                  (svm_scores, svm_times),
                  ({'linestyle': 'solid'}, {'linestyle': 'dashed'}))

#Plot performance for svms, and approximate kernel models
for ax, performance, svm_data, lineargs in plot_layout:
    plot_svm_data(ax, sample_sizes, performance, **lineargs)
    plot_svm_data(ax, sample_sizes, svm_data, **lineargs)
        
#Make fine-grained tweaks to plots.
accuracy, timescale = subplots

# plot verticle line.
accuracy.plot([data.shape[1]]*2, [0.7, 1], label="n_features")

# format legends and axes
accuracy.set_title("Classification accuracy")
accuracy.set_xlim(sample_sizes[0], sample_sizes[-1])
accuracy.set_xticks(())
accuracy.set_ylim(np.min(scores[approximate_kernel_labels[1]]), 1)
accuracy.set_ylabel("Classification accuracy")
accuracy.legend(loc='best')

timescale.set_xlabel("Sampling steps = transformed feature dimension")
timescale.set_title("Training times")
timescale.set_ylabel("Training time in seconds")
timescale.legend(loc='best')
plt.tight_layout()

# visualize the decision surface, projected down to the first
# two principal components of the dataset
pca = PCA(n_components=8, random_state=1).fit(data_train)
multiples = np.arange(-2, 2, 0.1)
X = pca.transform(data_train)

flat_grid = flat_grid_from_pca(pca, multiples)

fig, subplots = plt.subplots(1, 3, figsize=(12, 5))

# title for the plots
titles = ['SVC with rbf kernel',
          'SVC (linear kernel)\n with Fourier rbf feature map\n'
          'n_components=100',
          'SVC (linear kernel)\n with Nystroem rbf feature map\n'
          'n_components=100']

# predict and plot
classifiers = [svms['rbf svm'], approx_svm[approximate_kernel_labels[0]],
               approx_svm[approximate_kernel_labels[1]]]
for ax, title, clf in zip(subplots.flatten(), titles, classifiers):
    plot_projected_decision_surface(ax, clf, X, multiples, flat_grid)
    ax.scatter(X[:, 0], X[:, 1], c=targets_train, cmap=plt.cm.Paired, s=20,
               edgecolors='k')
    ax.set_title(title)
    ax.axis('off')
plt.tight_layout()
plt.show()
