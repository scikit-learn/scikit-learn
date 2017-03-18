# coding: utf-8
"""
Large Margin Nearest Neighbor Classification
"""

# Author: John Chiotellis <johnyc.code@gmail.com>

# License: BSD 3 clause (C) John Chiotellis

from __future__ import print_function
import warnings


import os
import numpy as np
from scipy import sparse, optimize

from ..neighbors import KNeighborsClassifier
from ..metrics.pairwise import euclidean_distances
from ..utils import gen_batches
from ..utils.fixes import argpartition
from ..utils.multiclass import check_classification_targets
from ..utils.validation import check_is_fitted, check_array, check_X_y, \
    check_random_state
from ..exceptions import DataDimensionalityWarning


class LargeMarginNearestNeighbor(KNeighborsClassifier):
    """Distance Metric Learning for Large Margin Classification

    Large margin nearest neighbor classification (LMNN) is a machine learning
    algorithm for metric learning. It learns a (pseudo-)metric in a
    supervised fashion to improve the classification accuracy of the k-nearest
    neighbor rule.
    The main intuition behind LMNN is to learn a pseudometric under which all
    data instances in the training set are surrounded by at least k instances
    that share the same class label. If this is achieved, the leave-one-out
    error is minimized.
    This implementation follows closely Kilian Weinberger's MATLAB code found
    at <https://bitbucket.org/mlcircus/lmnn> which solves the unconstrained
    problem, finding a linear transformation with L-BFGS instead of solving the
    constrained problem that finds the globally optimal metric.


    Parameters
    ----------
    L : array-like
        Initial transformation in an array with shape (n_features_out,
        n_features_in).  If None `load` will be used to load a transformation
        from a file. (default: None)

    n_neighbors : int
        Number of target neighbors (default: 3)

    max_iter : int
        Maximum number of iterations in the optimization (default: 200)

    use_pca : bool
        Whether to use pca to warm-start the linear transformation.
        If False, the identity will be used. (default: True)

    tol : float
        Tolerance for the optimization  (default: 1e-5)

    n_features_out : int
        Preferred dimensionality of the inputs after the transformation.
        If None it is inferred from `use_pca` and `L`.(default: None)

    max_constr : int
        Maximum number of constraints to enforce per iteration (default: 10
        million).

    use_sparse : bool
        Whether to use a sparse or a dense matrix for the impostor-pairs
        storage. Using a sparse matrix, the distance to impostors is computed
        twice, but it is somewhat faster for larger data sets than using a
        dense matrix. With a dense matrix, the unique impostor pairs have to be
        identified explicitly (default: True).

    load : string
        A file path from which to load a linear transformation. If None, either
        identity or pca will be used based on `use_pca` (default: None).

    save : string
        A file path prefix to save intermediate linear transformations to.
        After every function call, it will be extended with the function call
        number and the `.npy` file extension. If None, nothing will be saved
        (default: None).

    disp : int, optional
        If zero, then no output.  If a positive number, then
        ``0 < disp < 99`` print also f and ``|proj g|`` every disp iterations;
        ``disp = 99``   print details of every iteration except n-vectors;
        ``disp = 100``  print also the changes of active set and final x;
        ``disp > 100``  print details of every iteration including x and g.

    random_state : int
        A seed for reproducibility of random state  (default: None).

    Attributes
    ----------
    L_ : array-like
        The linear transformation used during fitting with shape
        (n_features_out, n_features_in).

    n_neighbors_ : int
        The number of target neighbors (decreased if n_neighbors was not
        realistic for all classes).

    n_features_out_ : int
        The dimensionality of a vector after applying to it the linear
        transformation.

    X_ : array-like
        An array of training samples with shape (n_samples, n_features_in).

    y_ : array-like
        An array of training labels with shape (n_samples,).

    labels_: array-like
        An array of the uniquely appearing class labels with shape (n_classes,)
        and type object.

    classes_: array-like
        An array of the uniquely appearing class labels as integers with shape
        (n_classes,) and type int.

    targets_ : array-like
        An array of target neighbors for each sample with shape (n_samples,
        n_neighbors).

    grad_static_ : array-like
        An array of the gradient component caused by target neighbors, that
        stays fixed throughout the algorithm with shape (n_features_in,
        n_features_in).

    n_iter_ : int
        The number of iterations of the optimizer.

    n_funcalls_ : int
        The number of times the optimizer computes the loss and the gradient.

    details_ : dict
        A dictionary of information created by the L-BFGS optimizer during
        fitting.


    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import LargeMarginNearestNeighbor
    >>> lmnn = LargeMarginNearestNeighbor(n_neighbors=1)
    >>> lmnn.fit(X, y) # doctest: +ELLIPSIS
    LargeMarginNearestNeighbor(...)
    >>> print(lmnn.predict([[1.1]]))
    [0]
    >>> print(lmnn.predict_proba([[0.9]]))
    [[ 1.  0.]]


    References
    ----------
    .. [1] `Weinberger, K. Q.; Saul L. K. (2009). "Distance Metric Learning
    for Large Margin Classification". Journal of Machine Learning Research.
    10: 207–244.
    <http://jmlr.csail.mit.edu/papers/volume10/weinberger09a/weinberger09a.pdf>`_

    .. [2] `Wikipedia entry on Large Margin Nearest Neighbor
           <https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor>`_

    """

    def __init__(self, L=None, n_neighbors=3, n_features_out=None,
                 max_iter=200, tol=1e-5, use_pca=True, max_constr=int(1e7),
                 use_sparse=True, load=None, save=None, disp=0,
                 random_state=None):

        super(LargeMarginNearestNeighbor, self).__init__(
            n_neighbors=n_neighbors)

        # Parameters
        self.L = L
        self.n_features_out = n_features_out
        self.max_iter = max_iter
        self.tol = tol
        self.use_pca = use_pca
        self.max_constr = max_constr
        self.use_sparse = use_sparse
        self.load = load
        self.save = save
        self.disp = disp
        self.random_state = random_state

    def fit(self, X, y):
        """Find a linear transformation by optimization of the unconstrained
        problem, such that the k-nearest neighbor classification accuracy
        improves.

        Parameters
        ----------
        X : array-like
            An array of training samples with shape (n_samples, n_features_in).
        y : array-like
            An array of data labels with shape (n_samples,).

        Returns
        -------
        LargeMarginNearestNeighbor
            self

        """

        # Check inputs consistency
        self.X_, y = check_X_y(X, y, order='F')
        check_classification_targets(y)

        # Store the appearing classes and the class index for each sample
        self.labels_, self.y_ = np.unique(y, return_inverse=True)
        self.classes_ = np.arange(len(self.labels_))

        # Check that the number of neighbors is achievable for all classes
        self.n_neighbors_ = self.check_n_neighbors(self.y_)
        # TODO: Notify superclass KNeighborsClassifier that n_neighbors
        # might have changed to n_neighbors_
        # super().set_params(n_neighbors=self.n_neighbors_)

        # Initialize transformer
        self.L_, self.n_features_out_ = self._init_transformer()

        # Prepare for saving if needed
        if self.save is not None:
            save_dir, save_file = os.path.split(self.save)
            if save_dir != '' and not os.path.exists(save_dir):
                os.mkdir(save_dir)
            save_file = self.save + '_' + str(self.n_funcalls_)
            np.save(save_file, self.L_)

        # Find target neighbors (fixed)
        self.targets_ = self._select_target_neighbors()

        # Compute gradient component of target neighbors (constant)
        self.grad_static_ = self._compute_grad_static()

        # Initialize number of optimizer iterations and objective funcalls
        self.n_iter_ = 0
        self.n_funcalls_ = 0

        # Set some default values for the optimizer
        max_corr = 100
        max_funcalls = 2*self.max_iter

        # Call optimizer
        try:
            L, loss, info = optimize.fmin_l_bfgs_b(func=self._loss_grad,
                                                   x0=self.L_,
                                                   m=max_corr,
                                                   pgtol=self.tol,
                                                   maxiter=self.max_iter,
                                                   maxfun=max_funcalls,
                                                   disp=self.disp,
                                                   callback=self._cb)
        except TypeError:
            # Type Error caused in old versions of SciPy because of no
            # maxiter argument (<= 0.9).
            self.L_ = np.asfortranarray(self.L_)
            try:
                L, loss, info = optimize.fmin_l_bfgs_b(func=self._loss_grad,
                                                       x0=self.L_,
                                                       m=max_corr,
                                                       maxfun=max_funcalls,
                                                       disp=self.disp,
                                                       pgtol=self.tol)
            except Exception:
                raise Exception('lbfgs does not work as expected in this '
                                'version of SciPy. Probably it is too old.')

        self.n_iter_ = info.get('nit', self.n_iter_)

        # Reshape result from optimizer
        self.L_ = L.reshape(self.n_features_out_, L.size //
                            self.n_features_out_)

        # Store output to return
        self.details_ = info
        self.details_['loss'] = loss

        # Fit a simple nearest neighbor classifier with the learned metric
        super(LargeMarginNearestNeighbor, self).fit(self.transform(), y)

        return self

    def transform(self, X=None):
        """Applies the learned transformation to the inputs.

        Parameters
        ----------
        X : array-like
            An array of data samples with shape (n_samples, n_features_in)
            (default: None, defined when fit is called).

        Returns
        -------
        array-like
            An array of transformed data samples with shape (n_samples,
            n_features_out).

        """
        if X is None:
            X = self.X_
        else:
            X = check_array(X)

        return X.dot(self.L_.T)

    def predict(self, X):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X : array-like, shape (n_query, n_features)
            Test samples.

        Returns
        -------
        y_pred : array of shape [n_query]
            Class labels for each data sample.
        """

        # Check if fit had been called
        check_is_fitted(self, ['X_', 'y_'])
        y_pred = super(LargeMarginNearestNeighbor, self).predict(
            self.transform(X))

        return y_pred

    def predict_proba(self, X):
        """Return probability estimates for the test data X.

        Parameters
        ----------
        X : array-like, shape (n_query, n_features)
            Test samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            of such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are ordered
            by lexicographic order.
        """

        # Check if fit had been called
        check_is_fitted(self, ['X_', 'y_'])
        probabilities = super(LargeMarginNearestNeighbor, self).predict_proba(
            self.transform(X))

        return probabilities

    def check_n_neighbors(self, y, n_neighbors=None):
        """Check if all classes have enough samples to query the specified
        number of neighbors."""

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        min_class_size = np.bincount(y).min()
        if min_class_size < 2:
            raise ValueError('At least one class has less than 2 ({}) '
                             'training samples.'.format(min_class_size))

        max_neighbors = min_class_size - 1
        if n_neighbors > max_neighbors:
            warnings.warn('n_neighbors(={}) too high. Setting to {}\n'.
                          format(n_neighbors, max_neighbors))

        return min(n_neighbors, max_neighbors)

    def _init_transformer(self):
        """Initialize the linear transformation by setting to user specified
        parameter, loading from a file, applying PCA or setting to identity."""

        if self.L is not None:
            L = self.L
        elif self.load is not None:
            L = np.load(self.load)
        elif self.use_pca and self.X_.shape[1] > 1:
            L = pca_fit(self.X_, return_transform=False)
        else:
            L = np.eye(self.X_.shape[1])

        n_features_out = L.shape[0] if self.n_features_out is None else \
            self.n_features_out
        n_features_in = self.X_.shape[1]

        if L.shape[1] != n_features_in:
            raise ValueError('Dimensionality of the given transformation and '
                             'the inputs don\'t match ({},{}).'.
                             format(L.shape[1], n_features_in))

        if n_features_out > n_features_in:
            warnings.warn('Outputs dimensionality ({}) cannot be larger than '
                          'inputs dimensionality, setting n_features_out to '
                          '{}!'.format(n_features_out, n_features_in),
                          DataDimensionalityWarning)
            n_features_out = n_features_in

        if L.shape[0] > n_features_out:
            L = L[:n_features_out]

        return L, n_features_out

    def _select_target_neighbors(self):
        """Find the target neighbors of each sample, that stay fixed during
        training.

        Returns
        -------
        array-like
            An array of neighbors indices for each sample with shape
            (n_samples, n_neighbors).

        """

        target_neighbors = np.empty((self.X_.shape[0], self.n_neighbors_),
                                    dtype=int)
        for class_ in self.classes_:
            class_ind, = np.where(np.equal(self.y_, class_))
            dist = euclidean_distances(self.X_[class_ind], squared=True)
            np.fill_diagonal(dist, np.inf)
            neigh_ind = argpartition(dist, self.n_neighbors_ - 1, axis=1)
            neigh_ind = neigh_ind[:, :self.n_neighbors_]
            # argpartition doesn't guarantee sorted order, so we sort again
            # but only the k neighbors
            row_ind = np.arange(len(class_ind))[:, None]
            neigh_ind = neigh_ind[row_ind,
                                  np.argsort(dist[row_ind, neigh_ind])]
            target_neighbors[class_ind] = class_ind[neigh_ind]

        return target_neighbors

    def _compute_grad_static(self):
        """Compute the gradient component due to the target neighbors that
        stays fixed throughout training

        Returns
        -------
        array-like
            An array with the sum of all weighted outer products with shape
            (n_features_in, n_features_in).

        """

        n_samples, n_neighbors = self.targets_.shape
        rows = np.repeat(np.arange(n_samples), n_neighbors)
        cols = self.targets_.flatten()
        targets_sparse = sparse.csr_matrix((np.ones(n_samples * n_neighbors),
                                            (rows, cols)),
                                           shape=(n_samples, n_samples))

        return sum_outer_products(self.X_, targets_sparse)

    def _cb(self, L):
        """Callback function called after every iteration of the optimizer.
        The intermediate transformations are saved to files if a valid
        `save` parameter was passed.

        Parameters
        ----------
        L : array-like
            The (flattened) linear transformation in the current iteration.

        """
        if self.save is not None:
            save_file = self.save + '_' + str(self.n_iter_)
            L = L.reshape(self.n_features_out_, L.size // self.n_features_out_)
            np.save(save_file, L)

        self.n_iter_ += 1

    def _loss_grad(self, L, X):
        """Compute the loss under a given linear transformation `L` and the
        loss gradient w.r.t. `L`.

        Parameters
        ----------
        L : array-like
            The current (flattened) linear transformation with shape
            (n_features_out x n_features_in,).

        Returns
        -------
        tuple
            float: The new loss.
            array-like: The new (flattened) gradient with shape
            (n_features_out x n_features_in,).

        """

        n_samples, n_features_in = X.shape
        self.L_ = L.reshape(self.n_features_out_, n_features_in)
        self.n_funcalls_ += 1

        Lx = self.transform()

        # Compute distances to target neighbors under L (plus margin)
        dist_tn = np.zeros((n_samples, self.n_neighbors_))
        for k in range(self.n_neighbors_):
            dist_tn[:, k] = np.sum(np.square(Lx - Lx[self.targets_[:, k]]),
                                   axis=1) + 1

        # Compute distances to impostors under L
        margin_radii = np.add(dist_tn[:, -1], 2)

        imp1, imp2, dist_imp = self._find_impostors(Lx, margin_radii,
                                                    use_sparse=self.use_sparse)

        loss = 0
        A0 = sparse.csr_matrix((n_samples, n_samples))
        for k in reversed(range(self.n_neighbors_)):
            loss1 = np.maximum(dist_tn[imp1, k] - dist_imp, 0)
            act, = np.where(loss1 != 0)
            A1 = sparse.csr_matrix((2*loss1[act], (imp1[act], imp2[act])),
                                   (n_samples, n_samples))

            loss2 = np.maximum(dist_tn[imp2, k] - dist_imp, 0)
            act, = np.where(loss2 != 0)
            A2 = sparse.csr_matrix((2*loss2[act], (imp1[act], imp2[act])),
                                   (n_samples, n_samples))

            vals = np.squeeze(np.asarray(A2.sum(0) + A1.sum(1).T))
            A0 = A0 - A1 - A2 + sparse.csr_matrix(
                     (vals, (range(n_samples), self.targets_[:, k])),
                     (n_samples, n_samples))
            loss = loss + np.sum(loss1 ** 2) + np.sum(loss2 ** 2)

        grad_new = sum_outer_products(self.X_, A0, remove_zero=True)
        df = self.L_.dot(self.grad_static_ + grad_new)
        df *= 2
        loss = loss + (self.grad_static_ * (self.L_.T.dot(self.L_))).sum()

        return loss, df.flatten()

    def _find_impostors(self, Lx, margin_radii, use_sparse=True):
        """Compute all impostor pairs exactly.

        Parameters
        ----------
        Lx : array-like
            An array of transformed samples with shape (n_samples,
            n_features_out).
        margin_radii : array-like
            An array of distances to the farthest target neighbors + margin,
            with shape (n_samples,).
        use_sparse : bool
            Whether to use a sparse matrix for storing the impostor pairs
            (default: True).

        Returns
        -------
        tuple: (array-like, array-like, array-like)

            imp1 : array-like
                An array of sample indices with shape (n_impostors,).
            imp2 : array-like
                An array of sample indices that violate a margin with shape
                (n_impostors,).
            dist : array-like
                An array of pairwise distances of (imp1, imp2) with shape
                (n_impostors,).

        """
        n_samples = Lx.shape[0]

        if use_sparse:
            # Initialize impostors matrix
            impostors_sp = sparse.csr_matrix((n_samples, n_samples),
                                             dtype=np.int8)

            for class_ in self.classes_[:-1]:
                imp1, imp2 = [], []
                ind_in, = np.where(np.equal(self.y_, class_))
                ind_out, = np.where(np.greater(self.y_, class_))

                # Subdivide ind_out x ind_in to chunks of a size that is
                # fitting in memory
                ii, jj = self._find_impostors_batch(
                    Lx[ind_out], Lx[ind_in], margin_radii[ind_out],
                    margin_radii[ind_in])
                if len(ii):
                    imp1.extend(ind_out[ii])
                    imp2.extend(ind_in[jj])
                    new_imps = sparse.csr_matrix(([1] * len(imp1),
                                                  (imp1, imp2)),
                                                 shape=(n_samples, n_samples),
                                                 dtype=np.int8)
                    impostors_sp = impostors_sp + new_imps

            imp1, imp2 = impostors_sp.nonzero()
            # subsample constraints if they are too many
            if impostors_sp.nnz > self.max_constr:
                random_state = check_random_state(self.random_state)
                ind_subsample = random_state.choice(impostors_sp.nnz,
                                                    self.max_constr,
                                                    replace=False)
                imp1, imp2 = imp1[ind_subsample], imp2[ind_subsample]

            dist = pairs_distances_batch(Lx, imp1, imp2)
        else:
            # Initialize impostors vectors
            imp1, imp2, dist = [], [], []
            for class_ in self.classes_[:-1]:
                ind_in, = np.where(np.equal(self.y_, class_))
                ind_out, = np.where(np.greater(self.y_, class_))

                # Subdivide idx_out x idx_in to chunks of a size that is
                # fitting in memory
                ii, jj, dd = self._find_impostors_batch(
                    Lx[ind_out], Lx[ind_in], margin_radii[ind_out],
                    margin_radii[ind_in], return_dist=True)
                if len(ii):
                    imp1.extend(ind_out[ii])
                    imp2.extend(ind_in[jj])
                    dist.extend(dd)

            ind_unique = unique_pairs(imp1, imp2, n_samples)

            # subsample constraints if they are too many
            if len(ind_unique) > self.max_constr:
                random_state = check_random_state(self.random_state)
                ind_unique = random_state.choice(ind_unique, self.max_constr,
                                                 replace=False)

            imp1 = np.asarray(imp1)[ind_unique]
            imp2 = np.asarray(imp2)[ind_unique]
            dist = np.asarray(dist)[ind_unique]

        return imp1, imp2, dist

    @staticmethod
    def _find_impostors_batch(x1, x2, t1, t2, return_dist=False,
                              batch_size=500):
        """Find impostor pairs in chunks to avoid large memory usage

        Parameters
        ----------
        x1 : array-like
            An array of transformed data samples with shape (n_samples,
            n_features).
        x2 : array-like
            An array of transformed data samples with shape (m_samples,
            n_features) where m_samples < n_samples.
        t1 : array-like
            An array of distances to the margins with shape (n_samples,).
        t2 : array-like
            An array of distances to the margins with shape (m_samples,).
        batch_size : int (Default value = 500)
            The size of each chunk of x1 to compute distances to.
        return_dist : bool (Default value = False)
            Whether to return the distances to the impostors.

        Returns
        -------
        tuple: (array-like, array-like, [array-like])

            imp1 : array-like
                An array of sample indices with shape (n_impostors,).
            imp2 : array-like
                An array of sample indices that violate a margin with shape
                (n_impostors,).
            dist : array-like, optional
                An array of pairwise distances of (imp1, imp2) with shape
                (n_impostors,).

        """

        n_samples = len(t1)
        imp1, imp2, dist = [], [], []
        for chunk in gen_batches(n_samples, batch_size):
            dist_out_in = euclidean_distances(x1[chunk], x2, squared=True)
            i1, j1 = np.where(dist_out_in < t1[chunk, None])
            i2, j2 = np.where(dist_out_in < t2[None, :])
            if len(i1):
                imp1.extend(i1 + chunk.start)
                imp2.extend(j1)
                if return_dist:
                    dist.extend(dist_out_in[i1, j1])
            if len(i2):
                imp1.extend(i2 + chunk.start)
                imp2.extend(j2)
                if return_dist:
                    dist.extend(dist_out_in[i2, j2])

        if return_dist:
            return imp1, imp2, dist
        else:
            return imp1, imp2


##########################
# Some helper functions #
#########################

def pca_fit(X, var_ratio=1, return_transform=True):
    """Do PCA and keep as many components as needed to explain the given
    variance ratio.

    Parameters
    ----------
    X : array-like
        An array of data samples with shape (n_samples, n_features).
    var_ratio : float
        The variance ratio to be captured (Default value = 1).
    return_transform : bool
        Whether to apply the transformation to the given data.

    Returns
    -------
    array-like
        If return_transform is True, an array with shape
        (n_samples, n_components) which is the input samples projected
        onto `n_components` principal components. Otherwise the first
        `n_components` eigenvectors of the covariance matrix corresponding to
        the `n_components` largest eigenvalues are returned as rows.

    """

    cov_ = np.cov(X, rowvar=False)  # Mean is removed
    evals, evecs = np.linalg.eigh(cov_)
    evecs = np.fliplr(evecs)

    if var_ratio == 1:
        L = evecs.T
    else:
        evals = np.flip(evals, axis=0)
        var_exp = np.cumsum(evals)
        var_exp = var_exp / var_exp[-1]
        n_components = np.argmax(np.greater_equal(var_exp, var_ratio))
        L = evecs.T[:n_components]

    if return_transform:
        return X.dot(L.T)
    else:
        return L


def sum_outer_products(X, weights, remove_zero=False):
    """Computes the sum of weighted outer products using a sparse weights
    matrix

    Parameters
    ----------
    X : array-like
        An array of data samples with shape (n_samples, n_features_in).
    weights : csr_matrix
        A sparse weights matrix (indicating target neighbors) with shape
        (n_samples, n_samples).
    remove_zero : bool
        Whether to remove rows and columns of the symmetrized weights matrix
        that are zero (default: False).

    Returns
    -------
    array-like
        An array with the sum of all weighted outer products with shape
        (n_features_in, n_features_in).

    """
    weights_sym = weights + weights.T
    if remove_zero:
        _, cols = weights_sym.nonzero()
        ind = np.unique(cols)
        weights_sym = weights_sym.tocsc()[:, ind].tocsr()[ind, :]
        X = X[ind]

    n = weights_sym.shape[0]
    diag = sparse.spdiags(weights_sym.sum(axis=0), 0, n, n)
    laplacian = diag.tocsr() - weights_sym
    sodw = X.T.dot(laplacian.dot(X))

    return sodw


def pairs_distances_batch(X, ind_a, ind_b, batch_size=500):
    """Equivalent to  np.sum(np.square(x[ind_a] - x[ind_b]), axis=1)

    Parameters
    ----------
    X : array-like
        An array of data samples with shape (n_samples, n_features_in).
    ind_a : array-like
        An array of samples indices with shape (n_indices,).
    ind_b : array-like
        Another array of samples indices with shape (n_indices,).
    batch_size :
        Size of each chunk of X to compute distances for (default: 500)

    Returns
    -------
    array-like
        An array of pairwise distances with shape (n_indices,).

    """
    n_indices = len(ind_a)
    res = np.zeros(n_indices)
    for chunk in gen_batches(n_indices, batch_size):
        res[chunk] = np.sum(np.square(X[ind_a[chunk]] - X[ind_b[chunk]]),
                            axis=1)

    return res


def unique_pairs(ind_a, ind_b, n_samples=None):
    """Find the unique pairs contained in zip(ind_a, ind_b)

    Parameters
    ----------
    ind_a : list
        A list with indices of reference samples of length m.
    ind_b : list
        A list with indices of impostor samples of length m.
    n_samples : int, optional
        The total number of samples (= maximum sample index + 1). If None it
        will be inferred from the indices.

    Returns
    -------
    array-like
         An array of indices of unique pairs with shape (k,) where k <= m.

    """
    # First generate a hash array
    if n_samples is None:
        n_samples = max(np.max(ind_a), np.max(ind_b))

    h = np.array([i * n_samples + j for i, j in zip(ind_a, ind_b)],
                 dtype=np.uint32)

    # Get the indices of the unique elements in the hash array
    _, ind_unique = np.unique(h, return_index=True)

    return ind_unique
