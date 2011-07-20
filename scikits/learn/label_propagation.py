"""
Label propagation in the context of this module refers to a set of 
semisupervised classification algorithms. In the high level, these algorithms 
work by forming a fully-connected graph between all points given and solving 
for the steady-state distribution of labels at each point.

These algorithms perform very well in practice. The cost of running can be very
expensive, at approximately O(N^3) where N is the number of (labeled and
unlabeled) points. The theory (why they perform so well) is motivated by 
intuitions from random walk algorithms and geometric relationships in the data.
For more information see [1].

Model Features
--------------
Label clamping:
  The algorithm tries to learn distributions of labels over the dataset. In the 
  "Hard Clamp" mode, the true ground labels are never allowed to change. They 
  are clamped into position. In the "Soft Clamp" mode, they are allowed some 
  wiggle room, but some alpha of their original value will always be retained.
  Hard clamp is the same as soft clamping with alpha set to 1.

Kernel:
  A function which projects a vector into some higher dimensional space. See the
  documentation for SVMs for more info on kernels.

References
----------
[1] Yoshua Bengio, Olivier Delalleau, Nicolas Le Roux. In Semi-Supervised 
    Learning (2006), pp. 193-216
"""
import numpy as np
from .base import BaseEstimator, ClassifierMixin

# really low epsilon (we don't really want machine eps)
EPSILON = 1e-7
DEFAULT_SIGMA = 0.5

### Main classes
class BaseLabelPropagation(BaseEstimator, ClassifierMixin):
    """
    Base class for label propagation module.

    Parameters
    ----------
    kernel : function (array_1, array_2) -> float
      kernel function to use
    sigma : float
      parameter to initialize the kernel function
    alpha : float
      clamping factor

    max_iters : float
      change maximum number of iterations allowed
    convergence_threshold : float
      threshold to consider the system at steady state
    """

    _default_alpha = 1

    def __init__(self, kernel=None, sigma=None, alpha=None, max_iters=1000, convergence_threshold=1e-3):
        self.max_iters, self.convergence_threshold = max_iters, convergence_threshold
        if sigma is None:
            self.sigma = DEFAULT_SIGMA
        else:
            self.sigma = sigma

        if kernel is None:
            self.kernel = gen_gaussian_kernel(sigma=self.sigma)
        elif hasattr(kernel, '__call__'):
            self.kernel = kernel

        # clamping factor
        if alpha is None:
            self.alpha = self._default_alpha
        else:
            self.alpha = alpha

    def _build_graph(self):
        """
        Builds a matrix representing a fully connected graph between each point
        in the dataset

        This basic implementation creates a non-stochastic affinity matrix, so 
        class probability distributions will exceed 1 (normalization may be desired)
        """
        self._graph_matrix = compute_affinity_matrix(self._X, kernel=self.kernel)

    def predict(self, X):
        return [np.argmax(self.predict_proba(x)) for x in X]

    def predict_proba(self, x):
        s = 0
        ary = np.asanyarray(x)
        if len(ary.shape) == 1:
            ary.resize((ary.shape[0], 1))
        for xj, yj in zip(self._X, self._y):
            Wx = self.kernel(xj, ary)
            s += Wx * yj / (Wx + EPSILON)
        return s

    def fit(self, X, y, **params):
        """
        Fit a semi-supervised label propagation model based on input data 
        matrix X and corresponding label matrix Y. 

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_freatures]
          A {n_samples by n_samples} size matrix will be created from this 
          (keep dataset fewer than 2000 points)
        y : array, shape = [n_labeled_samples]
          n_labeled_samples <= n_samples
          All unlabeled samples will be transductively assigned labels

        Examples
        --------
        >>> samples = [[1,0,0], [0,1,1], [0,1,0], [1,1,0], [0,0,1]]
        >>> labels = [1, -1, -1]
        >>> label_prop_model = BaseLabelPropagation()
        >>> label_prop_model.fit(samples, labels)

        Warning
        -------
        Base class works, but intended for instruction. It will generate a n
        onstochastic affinity matrix.
        """
        self._set_params(**params)
        self._X = np.asanyarray(X)
        self._y = np.asanyarray(y)

        # actual graph construction (implementations should override this)
        self._build_graph()

        # how many labels
        if len(self._y.shape) > 1:
            (num_labels, num_classes) = self._y.shape 
        else:
            # this handles binary class input
            num_labels, num_classes = self._y.shape[0], 1

        y_p = np.zeros((self._X.shape[0], num_classes))
        self._y.resize((self._X.shape[0], num_classes))

        Y_orig = np.asanyarray(y)
        Y_orig.resize((num_labels, num_classes))

        # main loop (minor optimization for hard clamp: put if statement outside of for loop)
        max_iters = self.max_iters
        while not_converged(self._y, y_p, self.convergence_threshold) and max_iters > 0:
            y_p = self._y
            self._y = np.dot(self._graph_matrix, self._y)
            # clamp
            self._y[:num_labels] = self.alpha * Y_orig - (1 - self.alpha) * self._y[:num_labels]
            max_iters -= 1
        return self

class LabelPropagation(BaseLabelPropagation):
    """
    Original label propagation algorithm. Computes a basic affinity matrix and 
    uses hard clamping.
    """
    def _build_graph(self):
        affinity_matrix = compute_affinity_matrix(self._X)
        degree_matrix = np.map(np.sum, affinity_matrix) * np.identity(affinity_matrix.shape[0])
        deg_inv = np.linalg.inv(degree_matrix)
        
        aff_ideg = deg_inv * np.matrix(affinity_matrix)
        self._graph_matrix = aff_ideg

class LabelSpreading(BaseLabelPropagation):
    """
    Similar to the basic Label Propgation algorithm, but uses affinity matrix 
    based on the graph laplacian and uses soft clamping for labels.

    Parameters
    ----------
    alpha : float
      "clamping factor" or how much of the original labels you want to keep
    """
    _default_alpha = 0.8

    def _build_graph(self):
        """
        Graph matrix for Label Spreading uses the Graph Laplacian!
        """
        affinity_matrix = compute_affinity_matrix(self._X, diagonal=0)
        degree_matrix = np.map(np.sum, affinity_matrix) * np.identity(affinity_matrix.shape[0])
        deg_invsq = np.sqrt(np.linalg.inv(degree_matrix))

        laplacian = deg_invsq * np.matrix(affinity_matrix) * deg_invsq
        self._graph_matrix = laplacian

### Helper functions        

def gen_gaussian_kernel(sigma=DEFAULT_SIGMA):
    """ generate a basic gaussian kernel function """
    def gaussian_kernel(x1, x2):
        """
        computes the gaussian kernel between two input vectors
        """
        return np.exp( -np.linalg.norm(x1 - x2) ** 2 / sigma )
    return gaussian_kernel
# default gaussian
gaussian_kernel = gen_gaussian_kernel()

def compute_affinity_matrix(X, kernel=gaussian_kernel, diagonal=1):
    """
    affinity matrix from input matrix (representing a fully connected graph)
    """
    height = X.shape[0]
    aff_mat = np.zeros((height,height)) # square matrix
    for i in xrange(height):
        aff_mat[i,i] = diagonal
        for j in xrange(i+1, height):
            aff = kernel(X[i], X[j])
            aff_mat[i,j] = aff
            aff_mat[j,i] = aff
    return aff_mat

def not_converged(y, y_hat, threshold=1e-3):
    """basic convergence check"""
    return np.sum(np.sum(np.abs(np.asarray(y-y_hat)))) > threshold
