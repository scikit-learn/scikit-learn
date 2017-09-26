"""
Multi-dimensional Scaling (MDS)
"""

# author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD

import numpy as np
from numpy import linalg as la

import warnings

from ..base import BaseEstimator
from ..metrics import euclidean_distances
from ..utils import check_random_state, check_array, check_symmetric
from ..externals.joblib import Parallel
from ..externals.joblib import delayed
from ..isotonic import IsotonicRegression


def _smacof_single(dissimilarities, metric=True, n_components=2, init=None,
                   max_iter=300, verbose=0, eps=1e-3, random_state=None):
    """Computes multidimensional scaling using SMACOF algorithm

    Parameters
    ----------
    dissimilarities : ndarray, shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Must be symmetric.

    metric : boolean, optional, default: True
        Compute metric or nonmetric SMACOF algorithm.

    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities. If an
        ``init`` array is provided, this option is overridden and the shape of
        ``init`` is used to determine the dimensionality of the embedding
        space.

    init : ndarray, shape (n_samples, n_components), optional, default: None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.

    random_state : int, RandomState instance or None, optional, default: None
        The generator used to initialize the centers.  If int, random_state is
        the seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.

    Returns
    -------
    X : ndarray, shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).

    n_iter : int
        The number of iterations corresponding to the best stress.
    """
    dissimilarities = check_symmetric(dissimilarities, raise_exception=True)

    n_samples = dissimilarities.shape[0]
    random_state = check_random_state(random_state)

    sim_flat = ((1 - np.tri(n_samples)) * dissimilarities).ravel()
    sim_flat_w = sim_flat[sim_flat != 0]
    if init is None:
        # Randomly choose initial configuration
        X = random_state.rand(n_samples * n_components)
        X = X.reshape((n_samples, n_components))
    else:
        # overrides the parameter p
        n_components = init.shape[1]
        if n_samples != init.shape[0]:
            raise ValueError("init matrix should be of shape (%d, %d)" %
                             (n_samples, n_components))
        X = init

    old_stress = None
    ir = IsotonicRegression()
    for it in range(max_iter):
        # Compute distance and monotonic regression
        dis = euclidean_distances(X)

        if metric:
            disparities = dissimilarities
        else:
            dis_flat = dis.ravel()
            # dissimilarities with 0 are considered as missing values
            dis_flat_w = dis_flat[sim_flat != 0]

            # Compute the disparities using a monotonic regression
            disparities_flat = ir.fit_transform(sim_flat_w, dis_flat_w)
            disparities = dis_flat.copy()
            disparities[sim_flat != 0] = disparities_flat
            disparities = disparities.reshape((n_samples, n_samples))
            disparities *= np.sqrt((n_samples * (n_samples - 1) / 2) /
                                   (disparities ** 2).sum())

        # Compute stress
        stress = ((dis.ravel() - disparities.ravel()) ** 2).sum() / 2

        # Update X using the Guttman transform
        dis[dis == 0] = 1e-5
        ratio = disparities / dis
        B = - ratio
        B[np.arange(len(B)), np.arange(len(B))] += ratio.sum(axis=1)
        X = 1. / n_samples * np.dot(B, X)

        dis = np.sqrt((X ** 2).sum(axis=1)).sum()
        if verbose >= 2:
            print('it: %d, stress %s' % (it, stress))
        if old_stress is not None:
            if(old_stress - stress / dis) < eps:
                if verbose:
                    print('breaking at iteration %d with stress %s' % (it,
                                                                       stress))
                break
        old_stress = stress / dis

    return X, stress, it + 1


def smacof(dissimilarities, metric=True, n_components=2, init=None, n_init=8,
           n_jobs=1, max_iter=300, verbose=0, eps=1e-3, random_state=None,
           return_n_iter=False):
    """Computes multidimensional scaling using the SMACOF algorithm.

    The SMACOF (Scaling by MAjorizing a COmplicated Function) algorithm is a
    multidimensional scaling algorithm which minimizes an objective function
    (the *stress*) using a majorization technique. Stress majorization, also
    known as the Guttman Transform, guarantees a monotone convergence of
    stress, and is more powerful than traditional techniques such as gradient
    descent.

    The SMACOF algorithm for metric MDS can summarized by the following steps:

    1. Set an initial start configuration, randomly or not.
    2. Compute the stress
    3. Compute the Guttman Transform
    4. Iterate 2 and 3 until convergence.

    The nonmetric algorithm adds a monotonic regression step before computing
    the stress.

    Parameters
    ----------
    dissimilarities : ndarray, shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Must be symmetric.

    metric : boolean, optional, default: True
        Compute metric or nonmetric SMACOF algorithm.

    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities. If an
        ``init`` array is provided, this option is overridden and the shape of
        ``init`` is used to determine the dimensionality of the embedding
        space.

    init : ndarray, shape (n_samples, n_components), optional, default: None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.

    n_init : int, optional, default: 8
        Number of times the SMACOF algorithm will be run with different
        initializations. The final results will be the best output of the runs,
        determined by the run with the smallest final stress. If ``init`` is
        provided, this option is overridden and a single run is performed.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If multiple
        initializations are used (``n_init``), each run of the algorithm is
        computed in parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging. For ``n_jobs`` below -1,
        (``n_cpus + 1 + n_jobs``) are used. Thus for ``n_jobs = -2``, all CPUs
        but one are used.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.

    random_state : int, RandomState instance or None, optional, default: None
        The generator used to initialize the centers.  If int, random_state is
        the seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.

    return_n_iter : bool, optional, default: False
        Whether or not to return the number of iterations.

    Returns
    -------
    X : ndarray, shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).

    n_iter : int
        The number of iterations corresponding to the best stress. Returned
        only if ``return_n_iter`` is set to ``True``.

    Notes
    -----
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)
    """

    dissimilarities = check_array(dissimilarities)
    random_state = check_random_state(random_state)

    if hasattr(init, '__array__'):
        init = np.asarray(init).copy()
        if not n_init == 1:
            warnings.warn(
                'Explicit initial positions passed: '
                'performing only one init of the MDS instead of %d'
                % n_init)
            n_init = 1

    best_pos, best_stress = None, None

    if n_jobs == 1:
        for it in range(n_init):
            pos, stress, n_iter_ = _smacof_single(
                dissimilarities, metric=metric,
                n_components=n_components, init=init,
                max_iter=max_iter, verbose=verbose,
                eps=eps, random_state=random_state)
            if best_stress is None or stress < best_stress:
                best_stress = stress
                best_pos = pos.copy()
                best_iter = n_iter_
    else:
        seeds = random_state.randint(np.iinfo(np.int32).max, size=n_init)
        results = Parallel(n_jobs=n_jobs, verbose=max(verbose - 1, 0))(
            delayed(_smacof_single)(
                dissimilarities, metric=metric, n_components=n_components,
                init=init, max_iter=max_iter, verbose=verbose, eps=eps,
                random_state=seed)
            for seed in seeds)
        positions, stress, n_iters = zip(*results)
        best = np.argmin(stress)
        best_stress = stress[best]
        best_pos = positions[best]
        best_iter = n_iters[best]

    if return_n_iter:
        return best_pos, best_stress, best_iter
    else:
        return best_pos, best_stress


class MDS(BaseEstimator):
    """Multidimensional scaling

    Read more in the :ref:`User Guide <multidimensional_scaling>`.

    Parameters
    ----------
    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities.

    metric : boolean, optional, default: True
        If ``True``, perform metric MDS; otherwise, perform nonmetric MDS.

    n_init : int, optional, default: 4
        Number of times the SMACOF algorithm will be run with different
        initializations. The final results will be the best output of the runs,
        determined by the run with the smallest final stress.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If multiple
        initializations are used (``n_init``), each run of the algorithm is
        computed in parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging. For ``n_jobs`` below -1,
        (``n_cpus + 1 + n_jobs``) are used. Thus for ``n_jobs = -2``, all CPUs
        but one are used.

    random_state : int, RandomState instance or None, optional, default: None
        The generator used to initialize the centers.  If int, random_state is
        the seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.

    dissimilarity : 'euclidean' | 'precomputed', optional, default: 'euclidean'
        Dissimilarity measure to use:

        - 'euclidean':
            Pairwise Euclidean distances between points in the dataset.

        - 'precomputed':
            Pre-computed dissimilarities are passed directly to ``fit`` and
            ``fit_transform``.

    Attributes
    ----------
    embedding_ : array-like, shape (n_components, n_samples)
        Stores the position of the dataset in the embedding space.

    stress_ : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).


    References
    ----------
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    "Out-of-sample extensions for lle, isomap, mds, eigenmaps, and spectral
    clustering." Bengio, Y. et al. Advances in neural information processing
    systems 16 (2004): 177-184.

    """
    def __init__(self, n_components=2, metric=True, n_init=4,
                 max_iter=300, verbose=0, eps=1e-3, n_jobs=1,
                 random_state=None, dissimilarity="euclidean",
                 extendible=False):
        self.n_components = n_components
        self.dissimilarity = dissimilarity
        self.metric = metric
        self.n_init = n_init
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.extendible = extendible
        self.X_train = None
        self.D_XX = None

    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def fit(self, X, y=None, init=None):
        """Fit the model on X

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        y: Ignored

        init : ndarray, shape (n_samples,), optional, default: None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.
        """
        if not self.extendible:
            self.fit_transform(X, init=init)
        else:
            self.X_train = X
            if self.dissimilarity == 'precomputed':
                D = X
            elif self.dissimilarity == 'euclidean':
                D = euclidean_distances(X)
                self.D_XX = euclidean_distances(self.X_train, self.X_train)
            else:
                raise ValueError("Proximity must be 'precomputed' or"
                                 "'euclidean'."
                                 " Got %s instead" % str(self.dissimilarity))
                return self

            # Normalising similarities
            K = np.zeros(np.shape(D))
            n = len(D)
            D_sq = np.square(D)
            P = np.eye(n) - 1 / n * np.ones((n, n))
            K = -0.5 * np.dot(np.dot(P, D_sq), P)

            # Sorting e-vectors and e-values according to e-val
            e_vals, e_vecs = la.eigh(K)
            ind_sort = np.argsort(e_vals)[::-1]
            self.e_vecs = e_vecs[:, ind_sort]
            self.e_vals = e_vals[ind_sort]
        return self

    def fit_transform(self, X, y=None, init=None):
        """Fit the data from X, and returns the embedded coordinates

        Parameters
        ----------
        X : array, shape=[n_samples, n_features], or [n_samples, n_samples] \
                if dissimilarity='precomputed'

        init : ndarray, shape (n_samples,), optional, default: None
            Should only be used with the non-extendible MDS (SMACOF).
            If None, randomly chooses the initial configuration
            if ndarray, initialize the SMACOF algorithm with this array.
        """
        if self.extendible:
            if init is not None:
                raise ValueError("Init is only for the non-extendible MDS.")
            ret = self._fit_transform_ext(X)
        else:
            ret = self._fit_transform(X, init)
        return ret

    def transform(self, X):
        """Apply the transformation on X

        If dissimilarity is Euclidean, apply the transformation on X.
        If dissimilarity is precomputed, X is the similarity matrix to be used
        between new (out-of-sample) points with old ones.
        The new points (X if Euclidean, or with X similarity matrix if
        precomputed) are projected in the same space as the training set.

        Parameters
        ----------
        X : array, shape [n_samples, n_features], or \
                [n_samples, n_train_samples] if dissimilarity='precomputed'
            New data, where n_samples is the number of samples
            and n_features is the number of features for "euclidean"
            dissimilarity. Else, similarity matrix (e.g. Euclidean distances
            between new and training points).

            NB: similarity matrix has to be centered, use the
            make_euclidean_similarities function to create it.


        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """

        if not self.extendible:
            raise ValueError("Method only available if extendible is True.")
        if self.dissimilarity == 'precomputed':
            D_new = X
        elif self.dissimilarity == 'euclidean':
            if self.X_train is None or self.D_XX is None:
                raise ValueError("Extendible MDS with Euclidean Similarities "
                                 "must be fit first. use ``MDS.fit()``")
            else:
                D_aX = euclidean_distances(X, self.X_train)
                D_new = self.center_similarities(D_aX, self.D_XX)
        else:
            raise ValueError("Dissimilarity not set properly: 'precomputed' "
                             "and 'euclidean' allowed.")
        X_new = self._mds_project(D_new, k=self.n_components)
        return X_new

    def _fit_transform(self, X, init=None):
        X = check_array(X)
        if X.shape[0] == X.shape[1] and self.dissimilarity != "precomputed":
            warnings.warn("The MDS API has changed. ``fit`` now constructs an"
                          " dissimilarity matrix from data. To use a custom "
                          "dissimilarity matrix, set "
                          "``dissimilarity='precomputed'``.")

        if self.dissimilarity == "precomputed":
            self.dissimilarity_matrix_ = X
        elif self.dissimilarity == "euclidean":
            self.dissimilarity_matrix_ = euclidean_distances(X)
        else:
            raise ValueError("Proximity must be 'precomputed' or 'euclidean'."
                             " Got %s instead" % str(self.dissimilarity))

        self.embedding_, self.stress_, self.n_iter_ = smacof(
            self.dissimilarity_matrix_, metric=self.metric,
            n_components=self.n_components, init=init, n_init=self.n_init,
            n_jobs=self.n_jobs, max_iter=self.max_iter, verbose=self.verbose,
            eps=self.eps, random_state=self.random_state,
            return_n_iter=True)

        return self.embedding_

    def _fit_transform_ext(self, X):
        self.fit(X)
        X_new = np.dot(self.e_vecs[:, :self.n_components],
                       np.diag(np.sqrt(self.e_vals[:self.n_components])))
        return X_new

    def center_similarities(self, D_aX, D_XX):
        """Centers similarities D_aX around D_XX

        Parameters
        ----------
        D_aX : array, shape=[n_new_samples, n_train_samples]
            Dissimilarity matrix of new and training data.
        D_XX : array, shape=[n_train_samples, n_train_samples]
            Dissimilarity matrix of training data.

        Returns
        -------
        new_similarities : array-like, shape=[n_new_samples, n_train_samples]
        """
        D_aX = np.square(D_aX)
        D_XX = np.square(D_XX)
        N = len(D_XX)
        M = len(D_aX)
        I_NN = np.ones((N, N))
        I_MN = np.ones((M, N))
        Exp_XX = np.sum(D_XX) / N ** 2

        new_similarities = -0.5 * (D_aX - (np.dot(D_aX, I_NN) +
                                           np.dot(I_MN, D_XX)
                                           ) / N + Exp_XX)
        return new_similarities

    def _mds_project(self, new_similarities, k=None):
        if k is None:
            k = self.n_components

        e_projections = np.zeros((len(new_similarities), k))

        for i in range(len(new_similarities)):
            for j in range(k):
                e_projections[i, j] = ((np.dot(self.e_vecs[:, j],
                                               new_similarities[i]) /
                                        np.sqrt(self.e_vals[j])))
        e_projections = np.dot(new_similarities,
                               np.dot(self.e_vecs[:, :k],
                                      np.diag(1/np.sqrt(self.e_vals[:k]))))
        return e_projections
