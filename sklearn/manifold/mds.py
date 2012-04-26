"""
Multi-dimensional Scaling (MDS)
"""

# author: Nelle Varoquaux <nelle.varoquaux@gmail.com>

import numpy as np

import warnings

from ..base import BaseEstimator
from ..metrics import euclidean_distances
from ..utils import check_random_state
from ..externals.joblib import Parallel
from ..externals.joblib import delayed


def pool_adjacent_violators(distances, similarities):
    """
    Pool adjancent violators

    Computes an isotonic regression of distances on similarities.

    Parameters
    ----------
    distances: ndarray, shape (n, 1)
        array to fit

    similarities: ndarray, shape (n, 1)
        array on which to fit

    Returns
    -------
    distances: ndarray, shape (n, 1)
    """
    # First approach for ties: ignore them. The multidimensional scaling won't
    # enforce that points with equal similarity be at equal distance.
    indxs = np.lexsort((distances, similarities))

    new_blocks = range(len(indxs))

    block = []
    sort = True
    while sort:
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

    return distances


def _smacof_single(similarities, metric=True, out_dim=2, init=None,
           max_iter=300, verbose=0, eps=1e-3, random_state=None):
    """
    Computes multidimensional scaling using SMACOF algorithm

    Parameters
    ----------
    similarities: symmetric ndarray, shape [n * n]
        similarities between the points

    metric: boolean, optional, default: True
        compute metric or nonmetric SMACOF algorithm

    out_dim: int, optional, default: 2
        number of dimension in which to immerse the similarities
        overridden if initial array is provided.

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
    X: ndarray (n_samples, out_dim), float
               coordinates of the n_samples points in a out_dim-space

    stress_: float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)

    """
    n_samples = similarities.shape[0]
    random_state = check_random_state(random_state)

    if similarities.shape[0] != similarities.shape[1]:
        raise ValueError("similarities must be a square array (shape=%d)" % \
                            n_samples)

    if np.any(similarities != similarities.T):
        raise ValueError("similarities must be symmetric")

    sim_flat = ((1 - np.tri(n_samples)) * similarities).flatten()
    sim_flat_w = sim_flat[sim_flat != 0]
    if init is None:
        # Randomly choose initial configuration
        X = random_state.rand(n_samples * out_dim)
        X = X.reshape((n_samples, out_dim))
    else:
        # overrides the parameter p
        out_dim = init.shape[1]
        if n_samples != init.shape[0]:
            raise ValueError("init matrix should be of shape (%d, %d)" % \
                                 (n_samples, out_dim))
        X = init

    old_stress = None
    for it in range(max_iter):
        # Compute distance and monotonic regression
        dis = euclidean_distances(X)

        if metric:
            disparities = similarities
        else:
            dis_flat = dis.flatten()
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
        stress = ((dis.flatten() - \
                    disparities.flatten()) ** 2).sum() / 2

        # Update X using the Guttman transform
        ratio = disparities / dis
        ratio[np.isinf(ratio) | np.isnan(ratio)] = 0
        B = - ratio + np.diag(ratio.sum(axis=1))
        X = 1. / n_samples * np.dot(B, X)
        if verbose == 2:
            print 'it: %d, stress %s' % (it, stress)
        if old_stress is not None:
            if(old_stress - stress) < eps:
                if verbose:
                    print 'breaking at iteration %d with stress %s' % (it,
                                                                       stress)
                break
        old_stress = stress

    return X, stress


def smacof(similarities, metric=True, out_dim=2, init=None, n_init=8, n_jobs=1,
           max_iter=300, verbose=0, eps=1e-3, random_state=None):
    """
    Computes multidimensional scaling using SMACOF algorithm

    Parameters
    ----------
    similarities: symmetric ndarray, shape (n_samples, n_samples)
        similarities between the points

    metric: boolean, optional, default: True
        compute metric or nonmetric SMACOF algorithm

    out_dim: int, optional, default: 2
        number of dimension in which to immerse the similarities
        overridden if initial array is provided.

    init: {None or ndarray of shape (n_samples, out_dim)}
        if None, randomly chooses the initial configuration
        if ndarray, initialize the SMACOF algorithm with this array

    n_init: int, optional, default: 8
        Number of time the smacof algorithm will be run with different
        initialisation. The final results will be the best output of the
        n_init consecutive runs in terms of stress.

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
    X: ndarray (n_samples,out_dim)
        Coordinates of the n_samples points in a out_dim-space

    stress: float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points)
    """

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
                                        out_dim=out_dim,
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
                        similarities, metric=metric, out_dim=out_dim,
                        init=init, max_iter=max_iter,
                        verbose=verbose, eps=eps,
                        random_state=seed)
                for seed in seeds)
        positions, stress = zip(results)
        best = np.argmin(stress)
        best_stress = stress[best]
        best_pos = positions[best]
    return best_pos, best_stress


class MDS(BaseEstimator):
    """
    Multidimensional scaling

    Parameters
    ----------
    metric: boolean, optional, default: True
        compute metric or nonmetric SMACOF algorithm

    out_dim: int, optional, default: 2
        number of dimension in which to immerse the similarities
        overridden if initial array is provided.

    n_init: int, optional, default: 8
        Number of time the smacof algorithm will be run with different
        initialisation. The final results will be the best output of the
        n_init consecutive runs in terms of stress.

    max_iter: int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run

    verbose: int, optional, default: 0
        level of verbosity

    eps: float, optional, default: 1e-6
        relative tolerance w.r.t stress to declare converge

    Attributes
    ----------
    positions_: array-like, shape [out_dim, n_samples]
        Stores the position of the dataset in the embedding space

    stress_: float
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
    def __init__(self, out_dim=2, metric=True, n_init=8,
                 max_iter=300, verbose=0, eps=1e-3, n_jobs=1):
        self.out_dim = out_dim
        self.metric = metric
        self.n_init = n_init
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose
        self.n_jobs = n_jobs

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
        self.positions_, self.stress_ = smacof(X, metric=self.metric,
                                     out_dim=self.out_dim,
                                     init=init,
                                     n_init=self.n_init,
                                     max_iter=self.max_iter,
                                     verbose=self.verbose,
                                     eps=self.eps)
        return self
