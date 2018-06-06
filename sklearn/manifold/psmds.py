import numpy as np
from ..base import BaseEstimator
from ..utils import check_array, check_random_state

from _mds_fast import (
    distance_matrix,
    update_distance_matrix,
    c_pertub_error as best_pertubation,
    mse as mse1d,
    mse2 as mse2d,
)


def _radius_update(radius, error, prev_error, tolerance=1e-4):
    """Updates the search radius based on the error decrease

    If the error increases or (error_decrease / error < tolerance) the search
    radius is halved. Else it stays the same
    Parameters
    ----------
    radius: float
        Search radius in current iteration
    error: float
        Error in current iteration
    prev_error: float
        Error in previous iteration
    tolerance: float, optional, default: 1e-4
        Tolerance for relative decrease in error before halving the radius
    Returns
    -------
    radius: float
        The Updated search radius
    """
    if error >= prev_error or prev_error - error <= error * tolerance:
        return radius * 0.5
    return radius


def _point_sampling(points, keep_percent=1.0, turn=-1, recalculate_each=-1):
    """Sample a percentage of points to explore

    Parameters
    ----------
    points: ndarray, shape (n_samples, n_components)
        The data points in ``n_components``-space for the current iteration.
    keep_percent: float, optional, default: 1.0
        Percentage of points to explore
    turn: int, optional, default=-1
        The current iteration
    recalculate_each: int, optional, default=-1
        If this value is greater than 0, the algorithm explores all the data
        points every ``recalculate_each`` iterations

    Returns
    -------
        points: ndarray, (floor(n_samples * keep_percent), n_components)
            The sampled data points
    """
    if keep_percent > 1.0 or 1.0 - keep_percent < 1e-5:
        return points
    if turn > 0 and recalculate_each > 0 and turn % recalculate_each == 0:
        return points
    keep = int(points.shape[0] * keep_percent)
    return np.random.choice(points, size=keep, replace=False)


def pattern_search_mds(
        d_goal, init=None, n_components=2, starting_radius=1.0, max_iter=1000,
        radius_barrier=1e-3, explore_dim_percent=1.0, sample_points=1.0,
        radius_update_tolerance=1e-4, verbose=0, n_jobs=1, random_state=None):
    """Computes multidimensional scaling using the pattern search MDS
    algorithm.

    The pattern search MDS algorithm is a multidimensional scaling algorithm
    that is inspired by derivative-free optimization techniques (General
    Pattern Search). This algorithm has higher per-epoch computational
    complexity than SMACOF (O(N^2 * L) vs O(N^2)), but converges in a
    significantly smaller number of epochs. The goal is to minimize the
    objective function (MSE between produced and target distance matrices) by
    greedily moving each data point towards the direction that produces the
    maximum decrease in error. The points can be moved on the surface of a
    hypershpere of radius r, in the positive and negative directions of the
    axis in the orthonormal base of the space (2 * n_components search
    directions). The search radius is halved when no global error decrease is
    produced and the search is continued in the refined search space.

    Pattern search MDS can summarized by the following steps:
    1. Set an initial start configuration, randomly or not and define starting
       radius r
    2. Compute the error
    3. For all points
        i. Try to move each point in all search dimensions by r
        ii. Select the move (pertubation) that produces the maximum error
           decrease
    4. If the global error increases or the error decrease is to small cut the
       search radius r by half
    5. If the radius becomes too small terminate
    6. Go to 1.

    Parameters
    ----------
    d_goal: ndarray, shape (n_samples, n_samples)
        Pairwise dissimilarities between the points.
    init : ndarray, shape (n_samples, n_components), optional, default: None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.
    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities.
    starting_radius : float, optional, default: 1.0
        The radius around the data points where the algorithm will start to
        search. If set too small, convergence will be slow. If set too high
        the error will overshoot and will affect the speed of convergence.
        In a sense it is analogous to the learning rate for gradient descent
        based algorithms.
    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.
    radius_barrier: float, optional, default: 1e-3
        When search radius becomes smaller than radius_barrier the algorithm
        terminates.
    explore_dim_percent: float, optional, default: 1.0
        Percentage of the search dimensions the algorithm will explore
    sample_points: float, optional, default: 1.0
        Percentage of points the algorithm will try to move
    radius_update_tolerance: float, optional, default: 1e-4
        When the relative error decrease becomes smaller than
        radius_update_tolerance, the search radius is halved and the
        exploration is continued  in the reduced search space
    verbose : int, optional, default: 0
        Level of verbosity.
    n_jobs : int, optional, default: 1
        The number of threads to use for the computation. The search across
         dimensions is parallelized using OpenMP. This parameter specifies the
         number of threads to use. An upper bound for this parameter is
         2 * n_components
    random_state : int, RandomState instance or None, optional, default: None
        The generator used to initialize the centers.  If int, random_state is
        the seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.
    Returns
    -------
    xs : ndarray, shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.
    error : float
        The final value of the error (MSE between produced and target distance
        matrices)
    n_iter : int
        The number of iterations until convergence
    Notes
    -----
    "Pattern Search MDS" Paraskevopoulos, G.; Tzinis, E.
    under review for JMLR, (2018)
    """
    n_samples = d_goal.shape[0]
    random_state = check_random_state(random_state)
    xs = (init if init is not None
          else random_state.rand(n_samples, n_components))
    d_current = distance_matrix(xs)
    points = np.arange(xs.shape[0])

    radius = starting_radius
    turn = 0
    error = mse2d(d_goal, d_current)
    prev_error = np.Inf
    if verbose:
        print("Starting Error: {}".format(error))

    while turn <= max_iter and radius > radius_barrier:
        turn += 1
        radius = _radius_update(
            radius, error, prev_error, tolerance=radius_update_tolerance)
        prev_error = error
        filtered_points = _point_sampling(points, keep_percent=sample_points)
        for point in filtered_points:
            point_error = mse1d(d_goal[point], d_current[point])
            optimum_error, optimum_k, optimum_step = best_pertubation(
                xs, radius, d_current, d_goal, point,
                percent=explore_dim_percent, n_jobs=n_jobs)
            error -= (point_error - optimum_error)
            d_current = update_distance_matrix(
                xs, d_current, point, optimum_step, optimum_k)
            xs[point, optimum_k] += optimum_step
        if verbose >= 2:
            print("Turn {0}: Radius {1}: (prev, error decrease, error): "
                  "({2}, {3}, {4})"
                  .format(turn, radius, prev_error, prev_error - error, error))
    if verbose:
        print("Ending Error: {}".format(error))
    return xs, error, turn


class MDS(BaseEstimator):
    """Multidimensional scaling using General Pattern Search
    Read more in the :ref:`User Guide <multidimensional_scaling>`.
    Parameters
    ----------
    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities.
    starting_radius : float, optional, default: 1.0
        The radius around the data points where the algorithm will start to
        search. If set too small, convergence will be slow. If set too high
        the error will overshoot and will affect the speed of convergence.
        In a sense it is analogous to the learning rate for gradient descent
        based algorithms.
    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.
    radius_barrier: float, optional, default: 1e-3
        When search radius becomes smaller than radius_barrier the algorithm
        terminates.
    explore_dim_percent: float, optional, default: 1.0
        Percentage of the search dimensions the algorithm will explore
    sample_points: float, optional, default: 1.0
        Percentage of points the algorithm will try to move
    radius_update_tolerance: float, optional, default: 1e-4
        When the relative error decrease becomes smaller than
        radius_update_tolerance, the search radius is halved and the
        exploration is continued  in the reduced search space
    verbose : int, optional, default: 0
        Level of verbosity.
    n_jobs : int, optional, default: 1
        The number of threads to use for the computation. The search across
         dimensions is parallelized using OpenMP. This parameter specifies the
         number of threads to use. An upper bound for this parameter is
         2 * n_components
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
    error_ : float
        The final value of the error (MSE between the distance matrix of the
        real data points and the one of the final self.embedding_)
    n_iter_: int
        The number of iterations until convergence
    References
    ----------
    "Pattern Search MDS" Paraskevopoulos, G.; Tzinis, E.
    under review for JMLR, (2018)
    """
    def __init__(self,
                 n_components=2,
                 starting_radius=1.0,
                 max_iter=1000,
                 radius_barrier=1e-3,
                 explore_dim_percent=1.0,
                 sample_points=1.0,
                 radius_update_tolerance=1e-4,
                 verbose=0,
                 n_jobs=1,
                 random_state=None,
                 dissimilarity='euclidean'):
        self.radius_update_tolerance = radius_update_tolerance
        self.sample_points = sample_points
        self.n_components = n_components
        self.starting_radius = starting_radius
        self.max_iter = max_iter
        self.radius_barrier = radius_barrier
        self.explore_dim_percent = explore_dim_percent
        self.num_epochs = 0
        self.verbose = verbose
        self.random_state = random_state
        self.n_jobs = 1
        self.dissimilarity = dissimilarity
        self.n_jobs = n_jobs

    def fit_transform(self, X, init=None):
        X = X.astype(np.float64)
        X = check_array(X)
        if self.dissimilarity == 'euclidean':
            d_goal = distance_matrix(X)
        elif self.dissimilarity == 'precomputed':
            d_goal = X
        else:
            raise ValueError("Proximity must be 'precomputed' or 'euclidean'."
                             " Got %s instead" % str(self.dissimilarity))
        self.embedding_, self.error_, self.n_iter_ = pattern_search_mds(
            d_goal, init=init, n_components=self.n_components,
            starting_radius=self.starting_radius, max_iter=self.max_iter,
            sample_points=self.sample_points,
            explore_dim_percent=self.explore_dim_percent,
            radius_update_tolerance=self.radius_update_tolerance,
            radius_barrier=self.radius_barrier,
            n_jobs=self.n_jobs, verbose=self.verbose,
            random_state=self.random_state
        )
        return self.embedding_

    def fit(self, X, init=None):
        self.fit_transform(X, init=init)
        return self
