"""
Multi-dimensional Scaling (MDS)
"""

# author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# Licence: BSD

import numpy as np

import warnings

from ..base import BaseEstimator
from ..metrics import euclidean_distances
from ..utils import check_random_state, check_arrays
from ..externals.joblib import Parallel
from ..externals.joblib import delayed


def pool_adjacent_violators(distances, similarities, max_iter=300,
                            verbose=0):
    """
    Pool adjancent violators

    Computes an isotonic regression of distances on similarities.

    Parameters
    ----------
    distances: ndarray, shape (n, 1)
        array to fit

    similarities: ndarray, shape (n, 1)
        array on which to fit

    max_iter: int, optional, default:300
        Set the maximum number of iteration

    verbose: int, optional, default: 0
        set the level of verbosity

    Returns
    -------
    distances: ndarray, shape (n, 1)

    Notes
    -----
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)
    """
    # First approach for ties: ignore them. The multidimensional scaling won't
    # enforce that points with equal similarity be at equal distance.
    indxs = np.lexsort((distances, similarities))

    new_blocks = range(len(indxs))

    block = []
    sort = True
    for it in range(max_iter):
        sort = False
        blocks = new_blocks[:]
        new_blocks = []
        block = []
        dis = distances[indxs[blocks[:-1]]] <= distances[indxs[blocks[1:]]] + \
             np.finfo(np.float).resolution
        for i, element in enumerate(dis):
            if not element:
                sort = True
                block.append(blocks[i])
            elif element and block:
                tmp = np.arange(block[0], blocks[i + 1])
                distances[indxs[tmp]] = distances[indxs[tmp]].mean()
                new_blocks.append(block[0])
                block = []
            else:
                new_blocks.append(blocks[i])
        # The last element
        if block:
            tmp = np.arange(block[0], len(similarities))
            distances[indxs[tmp]] = distances[indxs[tmp]].mean()
            new_blocks.append(block[0])
        else:
            new_blocks.append(len(similarities) - 1)
        if not sort:
            if verbose:
                print "Breaking at iteration %d" % it
            break

    return distances


def _smacof_single(similarities, metric=True, n_components=2, init=None,
           max_iter=300, verbose=0, eps=1e-3, random_state=None):
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

    init: {None or ndarray}
        if None, randomly chooses the initial configuration
        if ndarray, initialize the SMACOF algorithm with this array

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

    """
    n_samples = similarities.shape[0]
    random_state = check_random_state(random_state)

    if similarities.shape[0] != similarities.shape[1]:
        raise ValueError("similarities must be a square array (shape=%d)" % \
                            n_samples)
    res = 100 * np.finfo(np.float).resolution
    if np.any((similarities - similarities.T) > res):
        raise ValueError("similarities must be symmetric")

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
            raise ValueError("init matrix should be of shape (%d, %d)" % \
                                 (n_samples, n_components))
        X = init

    old_stress = None
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
            disparities_flat = pool_adjacent_violators(dis_flat_w,
                                                       sim_flat_w)
            disparities = dis_flat.copy()
            disparities[sim_flat != 0] = disparities_flat
            disparities = disparities.reshape((n_samples, n_samples))
            disparities *= np.sqrt((n_samples * (n_samples - 1) / 2) / \
                           (disparities ** 2).sum())

        # Compute stress
        stress = ((dis.ravel() - \
                    disparities.ravel()) ** 2).sum() / 2

        # Update X using the Guttman transform
        dis[dis == 0] = 1e-5
        ratio = disparities / dis
        B = - ratio
        B[np.arange(len(B)), np.arange(len(B))] += ratio.sum(axis=1)
        X = 1. / n_samples * np.dot(B, X)

        dis = np.sqrt((X ** 2).sum(axis=1)).sum()
        if verbose == 2:
            print 'it: %d, stress %s' % (it, stress)
        if old_stress is not None:
            if(old_stress - stress / dis) < eps:
                if verbose:
                    print 'breaking at iteration %d with stress %s' % (it,
                                                                       stress)
                break
        old_stress = stress / dis

    return X, stress


def smacof(similarities, metric=True, n_components=2, init=None, n_init=8,
           n_jobs=1, max_iter=300, verbose=0, eps=1e-3, random_state=None):
    """
    Computes multidimensional scaling using SMACOF (Scaling by Majorizing a
    Complicated Function) algorithm

    The SMACOF algorithm is a multidimensional scaling algorithm: it minimizes
    a objective function, the *stress*, using a majorization technique. The
    Stress Majorization, also known as the Guttman Transform, guarantees a
    monotone convergence of Stress, and is more powerful than traditional
    technics such as gradient descent.

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

    init : {None or ndarray of shape (n_samples, n_components)}
        if None, randomly chooses the initial configuration
        if ndarray, initialize the SMACOF algorithm with this array

    n_init : int, optional, default: 8
        Number of time the smacof algorithm will be run with different
        initialisation. The final results will be the best output of the
        n_init consecutive runs in terms of stress.

    n_jobs : int, optional, default: 1

        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debuging. For n_jobs below -1,
        (n_cpus + 1 - n_jobs) are used. Thus for n_jobs = -2, all CPUs but one
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

    Returns
    -------
    X : ndarray (n_samples,n_components)
        Coordinates of the n_samples points in a n_components-space

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)

    Notes
    -----
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)
    """

    similarities, = check_arrays(similarities, sparse_format='dense')
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
            pos, stress = _smacof_single(similarities, metric=metric,
                                        n_components=n_components,
                                        init=init, max_iter=max_iter,
                                        verbose=verbose, eps=eps,
                                        random_state=random_state)
            if best_stress is None or stress < best_stress:
                best_stress = stress
                best_pos = pos.copy()
    else:
        seeds = random_state.randint(np.iinfo(np.int32).max, size=n_init)
        results = Parallel(n_jobs=n_jobs, verbose=max(verbose - 1, 0))(
            delayed(_smacof_single)(
                        similarities, metric=metric, n_components=n_components,
                        init=init, max_iter=max_iter,
                        verbose=verbose, eps=eps,
                        random_state=seed)
                for seed in seeds)
        positions, stress = zip(*results)
        best = np.argmin(stress)
        best_stress = stress[best]
        best_pos = positions[best]
    return best_pos, best_stress


class MDS(BaseEstimator):
    """
    Multidimensional scaling

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
        used at all, which is useful for debuging. For n_jobs below -1,
        (n_cpus + 1 - n_jobs) are used. Thus for n_jobs = -2, all CPUs but one
        are used.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.


    Attributes
    ----------
    ``embedding_`` : array-like, shape [n_components, n_samples]
        Stores the position of the dataset in the embedding space

    ``stress_`` : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)


    Notes
    -----
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    """
    def __init__(self, n_components=2, metric=True, n_init=4,
                 max_iter=300, verbose=0, eps=1e-3, n_jobs=1,
                 random_state=None):
        self.n_components = n_components
        self.metric = metric
        self.n_init = n_init
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state = None

    def fit(self, X, init=None, y=None):
        """
        Computes the position of the points in the embedding space

        Parameters
        ----------
        X: array, shape=[n_samples, n_samples], symetric
            Proximity matrice

        init: {None or ndarray, shape (n_samples,)}
            if None, randomly chooses the initial configuration
            if ndarray, initialize the SMACOF algorithm with this array
        """
        self.fit_transform(X, init=init)
        return self

    def fit_transform(self, X, init=None, y=None):
        """
        Fit the data from X, and returns the embedded coordinates

        Parameters
        ----------
        X: array, shape=[n_samples, n_samples], symetric
            Proximity matrice

        init: {None or ndarray, shape (n_samples,)}
            if None, randomly chooses the initial configuration
            if ndarray, initialize the SMACOF algorithm with this array

        """
        self.embedding_, self.stress_ = smacof(X, metric=self.metric,
                                     n_components=self.n_components,
                                     init=init,
                                     n_init=self.n_init,
                                     n_jobs=self.n_jobs,
                                     max_iter=self.max_iter,
                                     verbose=self.verbose,
                                     eps=self.eps,
                                     random_state=self.random_state)

        return self.embedding_
