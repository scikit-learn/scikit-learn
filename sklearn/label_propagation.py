# encoding=UTF-8
"""Graph based semi-supervised learning with label propagation algorithms

Label propagation in the context of this module refers to a set of
semi-supervised classification algorithms. In the high level, these algorithms
work by forming a fully-connected graph between all points given and solving
for the steady-state distribution of labels at each point. Using these
algorithms assumes that the data can be clustered across a lower dimensional
manifold.

These algorithms perform very well in practice. The cost of running can be very
expensive, at approximately O(N^3) where N is the number of (labeled and
unlabeled) points. The theory (why they perform so well) is motivated by
intuitions from random walk algorithms and geometric relationships in the data.
For more information see the references below.

This algorithm solves a convex optimization problem and will converge to one
global solution. The ordering of input labels will not change the solution.

The algorithms assume maximum entropy priors for unlabeled data in each case
of these algorithms. It may be desired to incorporate prior information in
light of some domain information.

LabelSpreading is recommended for a good general case semi-supervised solution.
LabelPropagation much easier to understand and intuitive, so it may be good
for debugging, feature selection, and graph analysis, but in the most general
case it will be outperformed by LabelSpreading.

Model Features
--------------
Label clamping:
  The algorithm tries to learn distributions of labels over the dataset. In the
  "Hard Clamp" mode, the true ground labels are never allowed to change. They
  are clamped into position. In the "Soft Clamp" mode, they are allowed some
  wiggle room, but some alpha of their original value will always be retained.
  Hard clamp is the same as soft clamping with alpha set to 1.

Kernel:
  A function which projects a vector into some higher dimensional space. See
  the documentation for SVMs for more info on kernels. Only RBF kernels are
  currently supported.

Example
-------
>>> from sklearn import datasets
>>> label_prop_model = LabelPropagation()
>>> iris = datasets.load_iris()
>>> random_unlabeled_points = np.where(np.random.random_integers(0, 1,
...        size=len(iris.target)))
>>> labels = np.copy(iris.target)
>>> labels[random_unlabeled_points] = -1
>>> label_prop_model.fit(iris.data, labels)
... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
LabelPropagation(alpha=1, gamma=20, kernel='rbf', max_iter=30, tol=0.001,
        unlabeled_identifier=-1)

Notes
-----
References:
[1] Yoshua Bengio, Olivier Delalleau, Nicolas Le Roux. In Semi-Supervised
Learning (2006), pp. 193-216
"""
import numpy as np
from .base import BaseEstimator, ClassifierMixin
from .metrics.pairwise import rbf_kernel
from .utils.fixes import dot_out

# Authors: Clay Woolam <clay@woolam.org>
# License: BSD


class BaseLabelPropagation(BaseEstimator, ClassifierMixin):
    """Base class for label propagation module.

    Parameters
    ----------
    kernel : string
      string identifier for kernel function to use
      only 'rbf' kernel is currently supported

    gamma : float
      parameter for rbf kernel

    alpha : float
      clamping factor

    unlabeled_identifier : any object, same class as label objects
      a special identifier label that represents unlabeled examples
      in the training set

    max_iter : float
      change maximum number of iterations allowed

    tol : float
      threshold to consider the system at steady state
    """

    def __init__(self, kernel='rbf', gamma=20, alpha=1,
            unlabeled_identifier=-1, max_iter=30,
            tol=1e-3):

        self.max_iter = max_iter
        self.tol = tol

        # object referring to a point that is unlabeled
        self.unlabeled_identifier = unlabeled_identifier

        # kernel parameters
        self.kernel = kernel
        self.gamma = gamma

        # clamping factor
        self.alpha = alpha

    def _get_kernel(self, X, Y):
        if self.kernel == "rbf":
            return rbf_kernel(X, Y, gamma=self.gamma)
        else:
            raise ValueError("%s is not a valid kernel. Only rbf \
                             supported at this time" % self.kernel)

    def _build_graph(self):
        raise NotImplementedError("Graph construction must be implemented \
                to fit a label propagation model.")

    def predict(self, X):
        """Performs inductive inference across the model.

        Parameters
        ----------
        X : array_like, shape = [n_samples, n_features]

        Returns
        -------
        y : array_like, shape = [n_samples]
            Predictions for input data
        """
        probas = self.predict_proba(X)
        return self.unique_labels_[np.argmax(probas, axis=1)].flatten()

    def predict_proba(self, X):
        """Predict probability for each possible outcome.

        Compute the probability estimates for each single sample in X
        and each possible outcome seen during training (categorical
        distribution).

        Parameters
        ----------
        X : array_like, shape = (n_features, n_features)

        Return
        ------
        inference : array of normalized probability distributions across class
        labels
        """
        X_2d = np.atleast_2d(X)
        inference = np.dot(self._get_kernel(self._X, X_2d).T,
                self.label_distributions)
        normalizer = np.atleast_2d(np.sum(inference, axis=1)).T
        np.divide(inference, normalizer, out=inference)
        return inference

    def fit(self, X, y):
        """Fit a semi-supervised label propagation model on X and y.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_freatures]
          A {n_samples by n_samples} size matrix will be created from this
          (keep dataset fewer than 2000 points)

        y : array_like, shape = [n_samples]
          Signal to predict with unlabeled points marked with a special
          identifier. All unlabeled samples will be transductively assigned
          labels

        Returns
        -------
        Updated LabelPropagation object with a new transduction results
        """
        self._X = np.asanyarray(X)

        # actual graph construction (implementations should override this)
        graph_matrix = self._build_graph()

        # label construction
        # construct a categorical distribution for classification only
        unique_labels = np.unique(y)
        unique_labels = unique_labels[unique_labels !=\
                self.unlabeled_identifier]
        self.unique_labels_ = unique_labels

        n_samples, n_classes = len(y), len(unique_labels)

        y = np.asanyarray(y)
        unlabeled = y == self.unlabeled_identifier
        clamp_weights = np.ones((n_samples, 1))
        clamp_weights[unlabeled, 0] = self.alpha

        # initialize distributions
        self.label_distributions_ = np.zeros((n_samples, n_classes))
        for label in unique_labels:
            self.label_distributions_[y == label, unique_labels == label] = 1

        Y_static = np.copy(self.label_distributions_)
        if self.alpha > 0.:
            Y_static = Y_static * (1 - self.alpha)
        Y_static[unlabeled] = 0

        l_previous = np.zeros((self._X.shape[0], n_classes))
        self.label_distributions_.resize((self._X.shape[0], n_classes))

        remaining_iter = self.max_iter
        while _not_converged(self.label_distributions_, l_previous, self.tol)\
                and remaining_iter > 1:
            l_previous = self.label_distributions_
            self.label_distributions_ = np.dot(graph_matrix,
                    self.label_distributions_)
            # clamp
            self.label_distributions_ = np.multiply(clamp_weights,
                    self.label_distributions_) + Y_static
            remaining_iter -= 1

        normalizer = np.atleast_2d(np.sum(self.label_distributions, axis=1)).T
        np.divide(self.label_distributions, normalizer,
                out=self.label_distributions)
        # set the transduction item
        transduction = self.unique_labels_[np.argmax(self.label_distributions,
                axis=1)]
        self.transduction_ = transduction.flatten()
        return self


class LabelPropagation(BaseLabelPropagation):
    """Semi-supervised learning using Label Spreading strategy.

    Baseline semi-supervised estimator using a stochastic affinity
    matrix and hard clamping.

    Samples with missing label information must be assigned a special
    marker (the -1 integer by default) instead of the usual label value.

    Parameters
    ----------
    kernel : string
      String identifier for kernel function to use.  Only 'rbf' kernel
      is currently supported

    gamma : float
      parameter for rbf kernel

    alpha : float
      clamping factor

    unlabeled_identifier : any object, same class as label objects
      a special identifier label that represents unlabeled examples
      in the training set

    max_iter : float
      change maximum number of iterations allowed

    tol : float
      threshold to consider the system at steady state

    Examples
    --------
    >>> from sklearn import datasets
    >>> label_prop_model = LabelPropagation()
    >>> iris = datasets.load_iris()
    >>> random_unlabeled_points = np.where(np.random.random_integers(0, 1,
    ...    size=len(iris.target)))
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> label_prop_model.fit(iris.data, labels)
    ... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    LabelPropagation(alpha=1, gamma=20, kernel='rbf', max_iter=30, tol=0.001,
            unlabeled_identifier=-1)

    References
    ----------
    Xiaojin Zhu and Zoubin Ghahramani. Learning from labeled and unlabeled data
    with label propagation. Technical Report CMU-CALD-02-107, Carnegie Mellon
    University, 2002 http://pages.cs.wisc.edu/~jerryzhu/pub/CMU-CALD-02-107.pdf

    See Also
    --------
    LabelSpreading : Alternate label proagation strategy more robust to noise
    """
    def _build_graph(self):
        """
        Builds a matrix representing a fully connected graph between each point
        in the dataset.

        This basic implementation creates a non-stochastic affinity matrix, so
        class distributions will exceed 1 (normalization may be desired)
        """
        affinity_matrix = self._get_kernel(self._X, self._X)
        degree_inverse = np.diag(1. / np.sum(affinity_matrix, axis=0))
        dot_out(degree_inverse, affinity_matrix, out=affinity_matrix)
        return affinity_matrix


class LabelSpreading(BaseLabelPropagation):
    """Semi-supervised learning using Label Spreading strategy.

    Similar to the basic Label Propgation algorithm, but uses affinity matrix
    based on the graph laplacian and soft clamping accross the labels. Will be
    more robust to noise & uncertainty in the input labeling.

    Samples with missing label information must be assigned a special
    marker (the -1 integer by default) instead of the usual label value.

    Parameters
    ----------
    kernel : string
      string identifier for kernel function to use
      only 'rbf' kernel is currently supported

    gamma : float
      parameter for rbf kernel

    alpha : float
      clamping factor

    unlabeled_identifier : any object, same class as label objects
      a special identifier label that represents unlabeled examples
      in the training set

    max_iter : float
      change maximum number of iterations allowed

    tol : float
      threshold to consider the system at steady state

    Examples
    --------
    >>> from sklearn import datasets
    >>> label_prop_model = LabelSpreading()
    >>> iris = datasets.load_iris()
    >>> random_unlabeled_points = np.where(np.random.random_integers(0, 1,
    ...    size=len(iris.target)))
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> label_prop_model.fit(iris.data, labels)
    ... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    LabelSpreading(alpha=0.2, gamma=20, kernel='rbf', max_iter=30, tol=0.001,
           unlabeled_identifier=-1)

    References
    ----------
    Dengyong Zhou, Olivier Bousquet, Thomas Navin Lal, Jason Weston,
    Bernhard Schölkopf. Learning with local and global consistency (2004)
    http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.115.3219

    See Also
    --------
    Label Propagation : Unregularized graph based semi-supervised learning
    """

    def __init__(self, kernel='rbf', gamma=20, alpha=0.2,
            unlabeled_identifier=-1, max_iter=30, tol=1e-3):
        # this one has different base parameters
        super(LabelSpreading, self).__init__(kernel=kernel, gamma=gamma,
                alpha=alpha, unlabeled_identifier=unlabeled_identifier,
                max_iter=max_iter, tol=tol)

    def _build_graph(self):
        """Graph matrix for Label Spreading computes the graph laplacian"""
        # compute affinity matrix (or gram matrix)
        n_samples = self._X.shape[0]
        affinity_matrix = self._get_kernel(self._X, self._X)
        affinity_matrix[np.diag_indices(n_samples)] = 0
        degree_matrix = np.diag(1. / np.sum(affinity_matrix, axis=0))
        np.sqrt(degree_matrix, out=degree_matrix)
        dot_out(degree_matrix, affinity_matrix, out=affinity_matrix)

        # final step produces graph laplacian matrix
        dot_out(affinity_matrix, degree_matrix, out=affinity_matrix)
        return affinity_matrix


### Helper functions

def _not_converged(y_truth, y_prediction, tol=1e-3):
    """basic convergence check"""
    return np.sum(np.abs(np.asarray(y_truth - y_prediction))) > tol
