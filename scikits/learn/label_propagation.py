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
>>> from scikits.learn import datasets
>>> label_prop_model = LabelPropagation()
>>> iris = datasets.load_iris()
>>> random_unlabeled_points = np.where(np.random.random_integers(0, 1,
...        size=len(iris.target)))
>>> labels = np.copy(iris.target)
>>> labels[random_unlabeled_points] = -1
>>> label_prop_model.fit(iris.data, labels)
... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
LabelPropagation(...)

Notes
-----
References:
[1] Yoshua Bengio, Olivier Delalleau, Nicolas Le Roux. In Semi-Supervised
Learning (2006), pp. 193-216
"""
import numpy as np
from .base import BaseEstimator, ClassifierMixin
from .metrics.pairwise import rbf_kernel

# Authors: Clay Woolam <clay@woolam.org>


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

    max_iters : float
      change maximum number of iterations allowed
    tol : float
      threshold to consider the system at steady state
    """

    def __init__(self, kernel='rbf', gamma=20, alpha=1,
            unlabeled_identifier=-1, max_iters=30,
            tol=1e-3):

        self.max_iters = max_iters
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
        """
        Performs inductive inference across the model.

        Parameters
        ----------
        X : array_like, shape = [n_points, n_features]

        Returns
        -------
        y : array_like, shape = [n_points]
            Predictions for input data
        """
        ym = self.predict_proba(X)
        return self.unq_labels[np.argmax(ym, axis=1)].flatten()

    def predict_proba(self, X):
        """
        Returns a probability distribution (categorical distribution)
        over labels for a single input point.

        Parameters
        ----------
        X : array_like, shape = (n_features, n_features)
        """
        ary = np.atleast_2d(X)
        return np.dot(self._get_kernel(self._X, ary).T, self.y_)

    def fit(self, X, y):
        """
        Fit a semi-supervised label propagation model based on input data
        matrix X and corresponding label matrix Y.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_freatures]
          A {n_samples by n_samples} size matrix will be created from this
          (keep dataset fewer than 2000 points)
        y : array_like, shape = [n_labeled_samples]
          n_labeled_samples (unlabeled points marked with a special identifier)
          All unlabeled samples will be transductively assigned labels

        Returns
        -------
        updated LabelPropagation object with a new transduction results
        """
        self._X = np.asanyarray(X)

        # actual graph construction (implementations should override this)
        graph_matrix = self._build_graph()

        # label construction
        # construct a categorical distribution for classification only
        unq_labels = np.unique(y)
        unq_labels = unq_labels[unq_labels != self.unlabeled_identifier]
        self.unq_labels = unq_labels

        n_labels, n_classes = len(y), len(unq_labels)

        y_st = np.asanyarray(y)
        self.unlabeled_points = np.where(y_st == self.unlabeled_identifier)
        alpha_ary = np.ones((n_labels, 1))
        alpha_ary[self.unlabeled_points, 0] = self.alpha

        # initialize distributions
        self.y_ = np.zeros((n_labels, n_classes))
        for label in unq_labels:
            self.y_[y_st == label, unq_labels == label] = 1

        Y_alpha = np.copy(self.y_)
        Y_alpha = Y_alpha * (1 - self.alpha)
        Y_alpha[self.unlabeled_points] = 0

        y_p = np.zeros((self._X.shape[0], n_classes))
        self.y_.resize((self._X.shape[0], n_classes))

        max_iters = self.max_iters
        ct = self.tol
        while _not_converged(self.y_, y_p, ct) and max_iters > 1:
            y_p = self.y_
            self.y_ = np.dot(graph_matrix, self.y_)
            # clamp
            self.y_ = np.multiply(alpha_ary, self.y_) + Y_alpha
            max_iters -= 1

        # set the transduction item
        transduction = self.unq_labels[np.argmax(self.y_, axis=1)]
        self.transduction_ = transduction.flatten()
        return self


class LabelPropagation(BaseLabelPropagation):
    """
    Computes a stochastic affinity matrix and uses hard clamping.

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

    max_iters : float
      change maximum number of iterations allowed
    tol : float
      threshold to consider the system at steady state

    Examples
    --------
    >>> from scikits.learn import datasets
    >>> label_prop_model = LabelPropagation()
    >>> iris = datasets.load_iris()
    >>> random_unlabeled_points = np.where(np.random.random_integers(0, 1,
    ...    size=len(iris.target)))
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> label_prop_model.fit(iris.data, labels)
    ... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    LabelPropagation(...)

    References
    ----------
    Xiaojin Zhu and Zoubin Ghahramani. Learning from labeled and unlabeled data
    with label propagation. Technical Report CMU-CALD-02-107, Carnegie Mellon
    University, 2002. http://pages.cs.wisc.edu/~jerryzhu/pub/CMU-CALD-02-107.pdf

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
        degree_inv = np.diag(1. / np.sum(affinity_matrix, axis=0))
        np.dot(degree_inv, affinity_matrix, affinity_matrix)
        return affinity_matrix


class LabelSpreading(BaseLabelPropagation):
    """Semi-supervised learning using Label Spreading strategy.

    Similar to the basic Label Propgation algorithm, but uses affinity matrix
    based on the graph laplacian and soft clamping accross the labels. Will be
    more robust to noise & uncertainty in the input labeling.

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

    max_iters : float
      change maximum number of iterations allowed
    tol : float
      threshold to consider the system at steady state

    Examples
    --------
    >>> from scikits.learn import datasets
    >>> label_prop_model = LabelSpreading()
    >>> iris = datasets.load_iris()
    >>> random_unlabeled_points = np.where(np.random.random_integers(0, 1,
    ...    size=len(iris.target)))
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> label_prop_model.fit(iris.data, labels)
    ... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    LabelSpreading(...)

    References
    ----------
    Dengyong Zhou, Olivier Bousquet, Thomas Navin Lal, Jason Weston,
    Bernhard Schölkopf. Learning with local and global consistency (2004)
    http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.115.3219

    See Also
    --------
    Label Propagation : Unregularized graph based semi-supervised learning
    """

    def __init__(self, alpha=0.2, *args, **kwargs):
        # this one has different base parameters
        super(LabelSpreading, self).__init__(alpha=alpha, *args, **kwargs)

    def _build_graph(self):
        """Graph matrix for Label Spreading computes the graph laplacian"""
        # compute affinity matrix (or gram matrix)
        n_samples = self._X.shape[0]
        affinity_matrix = self._get_kernel(self._X, self._X)
        affinity_matrix[np.diag_indices(n_samples)] = 0
        degree_matrix = np.diag(1. / np.sum(affinity_matrix, axis=0))
        np.sqrt(degree_matrix, degree_matrix)
        np.dot(degree_matrix, affinity_matrix, affinity_matrix)
        # final step produces graph laplacian matrix
        np.dot(affinity_matrix, degree_matrix, affinity_matrix)
        return affinity_matrix

### Helper functions


def _not_converged(y, y_hat, tol=1e-3):
    """basic convergence check"""
    return np.sum(np.abs(np.asarray(y - y_hat))) > tol
