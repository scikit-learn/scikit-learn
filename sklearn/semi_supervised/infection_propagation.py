# coding=utf8
"""
Label propagation in the context of this module refers to a set of
semi-supervised classification algorithms. At a high level, these algorithms
work by forming a fully-connected graph between all points given and solving
for the steady-state distribution of labels at each point.

These algorithms perform very well in practice. The cost of running can be very
expensive, at approximately O(N^3) where N is the number of (labeled and
unlabeled) points. The theory (why they perform so well) is motivated by
intuitions from random walk algorithms and geometric relationships in the data.
For more information see the references below.

Model Features
--------------
Label clamping:
  The algorithm tries to learn distributions of labels over the dataset given
  label assignments over an initial subset. In one variant, the algorithm does
  not allow for any errors in the initial assignment (hard-clamping) while
  in another variant, the algorithm allows for some wiggle room for the initial
  assignments, allowing them to change by a fraction alpha in each iteration
  (soft-clamping).

Kernel:
  A function which projects a vector into some higher dimensional space. This
  implementation supports RBF and KNN kernels. Using the RBF kernel generates
  a dense matrix of size O(N^2). KNN kernel will generate a sparse matrix of
  size O(k*N) which will run much faster. See the documentation for SVMs for
  more info on kernels.

Examples
--------
>>> from sklearn import datasets
>>> from sklearn.semi_supervised import LabelPropagation
>>> label_prop_model = LabelPropagation()
>>> iris = datasets.load_iris()
>>> rng = np.random.RandomState(42)
>>> random_unlabeled_points = rng.rand(len(iris.target)) < 0.3
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

[2] Olivier Delalleau, Yoshua Bengio, Nicolas Le Roux. Efficient
Non-Parametric Function Induction in Semi-Supervised Learning. AISTAT 2005
"""

# Authors: Clay Woolam <clay@woolam.org>
#          Utkarsh Upadhyay <mail@musicallyut.in>
# License: BSD
from abc import ABCMeta, abstractmethod

import warnings
import numpy as np
from scipy import sparse
from copy import deepcopy
from sklearn.semi_supervised import LabelPropagation
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array
from sklearn.exceptions import ConvergenceWarning
class InfectionPropagation(LabelPropagation):
    """InfectionProp model for semi-supervised learning

        This model is based on a paper proposed in AISTAT 2018 based on 
        infection propagation. The main ideas is to use dynamic data models 
        to propagate infection labels from infected nodes to the uninfected ones.
        For more details you can take a look at a blog post on this here: 
        https://dadashkarimi.github.io/rosenfield-2018/

    Parameters
    ----------
    kernel : {'knn', 'rbf', callable}
        String identifier for kernel function to use or the kernel function
        itself. Only 'rbf' and 'knn' strings are valid inputs. The function
        passed should take two inputs, each of shape [n_samples, n_features],
        and return a [n_samples, n_samples] shaped weight matrix

    gamma : float
      parameter for rbf kernel

    n_neighbors : integer > 0
      parameter for knn kernel

    max_iter : integer
      maximum number of iterations allowed

    tol : float
      Convergence tolerance: threshold to consider the system at steady
      state

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

    Attributes
    ----------
    X_ : array, shape = [n_samples, n_features]
        Input array.

    classes_ : array, shape = [n_classes]
        The distinct labels used in classifying instances.

    label_distributions_ : array, shape = [n_samples, n_classes]
        Categorical distribution for each item.

    transduction_ : array, shape = [n_samples]
        Label assigned to each item via the transduction.

    n_iter_ : int
        Number of iterations run.

    Examples
    --------
    >>> from sklearn import datasets
    >>> from sklearn.semi_supervised import LabelSpreading
    >>> label_prop_model = InfectionProp(max_iter=10,kernel='knn')
    >>> iris = datasets.load_iris()
    >>> rng = np.random.RandomState(42)
    >>> random_unlabeled_points = rng.rand(len(iris.target)) < 0.3
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> label_prop_model.fit(iris.data, labels)
    ... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    InfectionProp(...)

    References
    ----------
    @inproceedings{DBLP:conf/aistats/RosenfeldG18,
      author    = {Nir Rosenfeld and
               Amir Globerson},
      title     = {Semi-Supervised Learning with Competitive Infection Models},
      booktitle = {International Conference on Artificial Intelligence and Statistics,
               {AISTATS} 2018, 9-11 April 2018, Playa Blanca, Lanzarote, Canary Islands,
               Spain},
      pages     = {336--346},
      year      = {2018},
      crossref  = {DBLP:conf/aistats/2018},
      url       = {http://proceedings.mlr.press/v84/rosenfeld18a.html},
      timestamp = {Sun, 15 Apr 2018 20:06:04 +0200},
      biburl    = {https://dblp.org/rec/bib/conf/aistats/RosenfeldG18},
      bibsource = {dblp computer science bibliography, https://dblp.org}
    }

    @author: dadashkarimi.github.io
        """

    _variant = 'spreading'

    def __init__(self, kernel='rbf', gamma=20, n_neighbors=7, alpha=0.2,
                 max_iter=30, tol=1e-3, n_jobs=1):

        # this one has different base parameters

        super(InfectionPropagation,self).__init__(kernel=kernel, gamma=gamma,
                                             n_neighbors=n_neighbors,
                                             alpha=alpha, max_iter=max_iter,
                                             tol=tol,
                                             n_jobs=n_jobs)

    def _build_graph(self):
        """Graph matrix for Label Spreading computes the graph laplacian"""
        # compute affinity matrix (or gram matrix)
        if self.kernel == 'knn':
            self.nn_fit = None
        n_samples = self.X_.shape[0]
        affinity_matrix = self._get_kernel(self.X_)
        laplacian = sparse.csgraph.laplacian(affinity_matrix, normed=True)
        laplacian = -laplacian
        if sparse.isspmatrix(laplacian):
            diag_mask = (laplacian.row == laplacian.col)
            laplacian.data[diag_mask] = 0.0
        else:
            laplacian.flat[::n_samples + 1] = 0.0  # set diag to 0.0
        return laplacian


    def fit(self, X, y):
        """Fit a semi-supervised label propagation model based

        All the input data is provided matrix X (labeled and unlabeled)
        and corresponding label matrix y with a dedicated marker value for
        unlabeled samples.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            A {n_samples by n_samples} size matrix will be created from this

        y : array_like, shape = [n_samples]
            n_labeled_samples (unlabeled points are marked as -1)
            All unlabeled samples will be transductively assigned labels

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_X_y(X, y)
        self.X_ = X
        check_classification_targets(y)
        # actual graph construction (implementations should override this)
        graph_matrix = self._build_graph()
        affinity_matrix = graph_matrix# self._get_kernel(self.X_)
        # label construction
        # construct a categorical distribution for classification only
        classes = np.unique(y)
        classes = (classes[classes != -1])
        self.classes_ = classes
        n_samples, n_classes = len(y), len(classes)
      
        alpha = self.alpha
        if self._variant == 'spreading' and \
                (alpha is None or alpha <= 0.0 or alpha >= 1.0):
            raise ValueError('alpha=%s is invalid: it must be inside '
                             'the open interval (0, 1)' % alpha)
        y = np.asarray(y)
        unlabeled = y == -1

        # initialize distributions
        self.label_distributions_ = np.zeros((n_samples, n_classes))
        V = self.X_.shape[0]
        for label in classes:
            self.label_distributions_[y == label, classes == label] = 1
        y_static = np.copy(self.label_distributions_)
        if self._variant == 'propagation':
            # LabelPropagation
            y_static[unlabeled] = 0
        else:
            # LabelSpreading
            y_static *= 1 - alpha

        l_previous = np.zeros((self.X_.shape[0], n_classes))
        np.expand_dims(unlabeled,axis=1)
        if sparse.isspmatrix(graph_matrix):
            graph_matrix = graph_matrix.tocsr()
        if sparse.isspmatrix(affinity_matrix):
            affinity_matrix = affinity_matrix.tocsr()
        import Queue as Q
        def D_theta(p_uv):
            X = np.random.multinomial(1, [p_uv,1-p_uv], size=1)
            return 1 if X[0][0] == 0 else np.inf
        def q_u(y_v):
            return 1.0/(y_v+1)
        for self.n_iter_ in range(self.max_iter):
            l_previous = deepcopy(self.label_distributions_)
            q = Q.PriorityQueue()
            dist = np.full(V, np.inf) # distance
            for j in np.argwhere(unlabeled==False)[:,0]:
                dist[j]=0
                q.put((dist[j],j))
            while not q.empty():
                dist_v, v= q.get()
                neighbors = sparse.find(affinity_matrix[v,:]!=0)[1]
                for u in neighbors:
                        delta_uv = D_theta(affinity_matrix[v,u]/affinity_matrix[v,:].sum())
                        if delta_uv == np.inf: # not infected
                                continue
                        alt = dist [v] + affinity_matrix[u,v] #+ q_u(y[v])
                        if alt < dist[u]:
                                dist[u] = alt
                                y[u] = y[v] # u inherits label from parent v
                                self.label_distributions_[u][classes.tolist().index(y[v])] +=1
                                q.put((dist[u],u))
            if np.abs(self.label_distributions_ - l_previous).sum() < self.tol:
                break
        
        else:
            warnings.warn(
                'max_iter=%d was reached without convergence.' % self.max_iter,
                category=ConvergenceWarning
            )
            self.n_iter_ += 1
        normalizer = np.sum(self.label_distributions_, axis=1)[:, np.newaxis]
        self.label_distributions_ /= normalizer

        # set the transduction item
        transduction = self.classes_[np.argmax(self.label_distributions_,
                                               axis=1)]
        self.transduction_ = transduction.ravel()
        return self


