"""
Multi-dimensional Scaling (MDS)
"""

# author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD

import numpy as np

import warnings

from ..base import BaseEstimator
from ..metrics import euclidean_distances
from ..utils import check_random_state, check_array, check_symmetric
from ..externals.joblib import Parallel
from ..externals.joblib import delayed
from ..isotonic import IsotonicRegression


def _e_vect(n, N):
    """
    Computes vector of length N, such that n-th element is one,
    and zero otherwise.

    Parameters
    ----------
    n: int
        n-th element of array is one

    N: int
        size of array

    Returns
    -------
    out: ndarray (N,)
        array with N zeros, n-th element is one

    """
    out = np.zeros((1, N)).T
    out[n, 0] = 1

    return np.array(out)


def _makeVinv(N, W=None):
    """
    Computes the pseudo-inverse of V matrix, used in the
    Guttman transform: X = Vinv * B * X

    Parameters
    ----------
    N: int
        Size of V matrix (in smacof, is n_samples)

    W: ndarray (N, N), optional, default None
        weighting matrix of similarities, default considers all weights to one

    """
    if W is None:
        Vinv = np.identity(N) - np.ones((N, N)) / N
        Vinv /= N
    else:
        V = np.zeros((N, N))
        for nn in range(N):
            for mm in range(nn, N):
                V += W[nn, mm] * np.dot(_e_vect(nn, N) - _e_vect(mm, N),
                                        (_e_vect(nn, N) - _e_vect(mm, N)).T)

        Vinv = np.linalg.pinv(V)

    return Vinv


def _smacof_single(similarities, metric=True, n_components=2,
                   init=None, weight=None, max_iter=300, verbose=0,
                   eps=1e-3, random_state=None):
    """
    Computes multidimensional scaling using SMACOF algorithm

    Parameters
    ----------
    similarities: symmetric ndarray, shape [n * n]
        similarities between the points

    metric: boolean, optional, default: True
        compute metric or nonmetric SMACOF algorithm

    n_components: int, optional, default: 2
        number of dimension in which to immerse the similarities
        overwritten if initial array is provided.

    init: {None or ndarray}, optional
        if None, randomly chooses the initial configuration
        if ndarray, initialize the SMACOF algorithm with this array

    weight: symmetric ndarray, shape [n * n], optional, default: None
        weighting matrix of similarities. Default considers all weights to 1.

    max_iter: int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run

    verbose: int, optional, default: 0
        level of verbosity

    eps: float, optional, default: 1e-6
        relative tolerance w.r.t stress to declare converge

    random_state: integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Returns
    -------
    X: ndarray (n_samples, n_components), float
               coordinates of the n_samples points in a n_components-space

    stress_: float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)

    n_iter : int
        Number of iterations run.

    """
    similarities = check_symmetric(similarities, raise_exception=True)

    n_samples = similarities.shape[0]
    random_state = check_random_state(random_state)

    sim_flat = ((1 - np.tri(n_samples)) * similarities).ravel()
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

    Vinv = _makeVinv(n_samples, W=weight)

    old_stress = None
    ir = IsotonicRegression()
    for it in range(max_iter):
        # Compute distance and monotonic regression
        dis = euclidean_distances(X)

        if metric:
            disparities = similarities
        else:
            dis_flat = dis.ravel()
            # similarities with 0 are considered as missing values
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
        if weight is None:
            ratio = disparities / dis
        else:
            ratio = weight * disparities / dis

        B = - ratio
        B[np.arange(len(B)), np.arange(len(B))] += ratio.sum(axis=1)
        X = np.dot(Vinv, np.dot(B, X))

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


def smacof(similarities, metric=True, n_components=2, init=None, weight=None,
           n_init=8, n_jobs=1, max_iter=300, verbose=0, eps=1e-3,
           random_state=None, return_n_iter=False):
    """
    Computes multidimensional scaling using SMACOF (Scaling by Majorizing a
    Complicated Function) algorithm

    The SMACOF algorithm is a multidimensional scaling algorithm: it minimizes
    a objective function, the *stress*, using a majorization technique. The
    Stress Majorization, also known as the Guttman Transform, guarantees a
    monotone convergence of Stress, and is more powerful than traditional
    techniques such as gradient descent.

    The SMACOF algorithm for metric MDS can summarized by the following steps:

    1. Set an initial start configuration, randomly or not.
    2. Compute the stress
    3. Compute the Guttman Transform
    4. Iterate 2 and 3 until convergence.

    The nonmetric algorithm adds a monotonic regression steps before computing
    the stress.

    Parameters
    ----------
    similarities : symmetric ndarray, shape (n_samples, n_samples)
        similarities between the points

    metric : boolean, optional, default: True
        compute metric or nonmetric SMACOF algorithm

    n_components : int, optional, default: 2
        number of dimension in which to immerse the similarities
        overridden if initial array is provided.

    init : {None or ndarray of shape (n_samples, n_components)}, optional
        if None, randomly chooses the initial configuration
        if ndarray, initialize the SMACOF algorithm with this array

    weight: symmetric ndarray, shape [n * n], optional, default: None
        weighting matrix of similarities. Default considers all weights to 1.

    n_init : int, optional, default: 8
        Number of time the smacof algorithm will be run with different
        initialisation. The final results will be the best output of the
        n_init consecutive runs in terms of stress.

    n_jobs : int, optional, default: 1

        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging. For n_jobs below -1,
        (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one
        are used.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run

    verbose : int, optional, default: 0
        level of verbosity

    eps : float, optional, default: 1e-6
        relative tolerance w.r.t stress to declare converge

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    return_n_iter : bool
        Whether or not to return the number of iterations.

    Returns
    -------
    X : ndarray (n_samples,n_components)
        Coordinates of the n_samples points in a n_components-space

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)

    n_iter : int
        The number of iterations corresponding to the best stress.
        Returned only if `return_n_iter` is set to True.

    Notes
    -----
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)
    """
    similarities = check_array(similarities)
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
                similarities, metric=metric,
                n_components=n_components, init=init, weight=weight,
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
                similarities, metric=metric, n_components=n_components,
                init=init, weight=weight, max_iter=max_iter, verbose=verbose,
                eps=eps, random_state=seed)
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
    metric : boolean, optional, default: True
        compute metric or nonmetric SMACOF (Scaling by Majorizing a
        Complicated Function) algorithm

    n_components : int, optional, default: 2
        number of dimension in which to immerse the similarities
        overridden if initial array is provided.

    n_init : int, optional, default: 4
        Number of time the smacof algorithm will be run with different
        initialisation. The final results will be the best output of the
        n_init consecutive runs in terms of stress.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run

    verbose : int, optional, default: 0
        level of verbosity

    eps : float, optional, default: 1e-6
        relative tolerance w.r.t stress to declare converge

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging. For n_jobs below -1,
        (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one
        are used.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    dissimilarity : string
        Which dissimilarity measure to use.
        Supported are 'euclidean' and 'precomputed'.


    Attributes
    ----------
    embedding_ : array-like, shape [n_components, n_samples]
        Stores the position of the dataset in the embedding space

    stress_ : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)


    References
    ----------
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    """
    def __init__(self, n_components=2, metric=True, n_init=4,
                 max_iter=300, verbose=0, eps=1e-3, n_jobs=1,
                 random_state=None, dissimilarity="euclidean"):
        self.n_components = n_components
        self.dissimilarity = dissimilarity
        self.metric = metric
        self.n_init = n_init
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state = random_state

    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def fit(self, X, y=None, init=None, weight=None):
        """
        Computes the position of the points in the embedding space

        Parameters
        ----------
        X : array, shape=[n_samples, n_features], or [n_samples, n_samples] \
                if dissimilarity='precomputed'
            Input data.

        init : {None or ndarray, shape (n_samples,)}, optional
            If None, randomly chooses the initial configuration
            if ndarray, initialize the SMACOF algorithm with this array.

        weight: symmetric ndarray, shape [n * n], optional, default: None
            weighting matrix of similarities. In default, all weights are 1.
        """
        self.fit_transform(X, init=init, weight=weight)
        return self

    def fit_transform(self, X, y=None, init=None, weight=None):
        """
        Fit the data from X, and returns the embedded coordinates

        Parameters
        ----------
        X : array, shape=[n_samples, n_features], or [n_samples, n_samples] \
                if dissimilarity='precomputed'
            Input data.

        init : {None or ndarray, shape (n_samples,)}, optional
            If None, randomly chooses the initial configuration
            if ndarray, initialize the SMACOF algorithm with this array.

        weight: symmetric ndarray, shape [n * n], optional, default: None
            weighting matrix of similarities. In default, all weights are 1.
        """
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
            n_components=self.n_components, init=init, weight=weight,
            n_init=self.n_init, n_jobs=self.n_jobs, max_iter=self.max_iter,
            verbose=self.verbose, eps=self.eps, random_state=self.random_state,
            return_n_iter=True)

        return self.embedding_
