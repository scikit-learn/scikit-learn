"""Algorithms to support fitting routines in seaborn plotting functions."""
from __future__ import division
import numpy as np
from scipy import stats
from .external.six.moves import range


def bootstrap(*args, **kwargs):
    """Resample one or more arrays with replacement and store aggregate values.

    Positional arguments are a sequence of arrays to bootstrap along the first
    axis and pass to a summary function.

    Keyword arguments:
        n_boot : int, default 10000
            Number of iterations
        axis : int, default None
            Will pass axis to ``func`` as a keyword argument.
        units : array, default None
            Array of sampling unit IDs. When used the bootstrap resamples units
            and then observations within units instead of individual
            datapoints.
        smooth : bool, default False
            If True, performs a smoothed bootstrap (draws samples from a kernel
            destiny estimate); only works for one-dimensional inputs and cannot
            be used `units` is present.
        func : callable, default np.mean
            Function to call on the args that are passed in.
        random_seed : int | None, default None
            Seed for the random number generator; useful if you want
            reproducible resamples.

    Returns
    -------
    boot_dist: array
        array of bootstrapped statistic values

    """
    # Ensure list of arrays are same length
    if len(np.unique(list(map(len, args)))) > 1:
        raise ValueError("All input arrays must have the same length")
    n = len(args[0])

    # Default keyword arguments
    n_boot = kwargs.get("n_boot", 10000)
    func = kwargs.get("func", np.mean)
    axis = kwargs.get("axis", None)
    units = kwargs.get("units", None)
    smooth = kwargs.get("smooth", False)
    random_seed = kwargs.get("random_seed", None)
    if axis is None:
        func_kwargs = dict()
    else:
        func_kwargs = dict(axis=axis)

    # Initialize the resampler
    rs = np.random.RandomState(random_seed)

    # Coerce to arrays
    args = list(map(np.asarray, args))
    if units is not None:
        units = np.asarray(units)

    # Do the bootstrap
    if smooth:
        return _smooth_bootstrap(args, n_boot, func, func_kwargs)

    if units is not None:
        return _structured_bootstrap(args, n_boot, units, func,
                                     func_kwargs, rs)

    boot_dist = []
    for i in range(int(n_boot)):
        resampler = rs.randint(0, n, n)
        sample = [a.take(resampler, axis=0) for a in args]
        boot_dist.append(func(*sample, **func_kwargs))
    return np.array(boot_dist)


def _structured_bootstrap(args, n_boot, units, func, func_kwargs, rs):
    """Resample units instead of datapoints."""
    unique_units = np.unique(units)
    n_units = len(unique_units)

    args = [[a[units == unit] for unit in unique_units] for a in args]

    boot_dist = []
    for i in range(int(n_boot)):
        resampler = rs.randint(0, n_units, n_units)
        sample = [np.take(a, resampler, axis=0) for a in args]
        lengths = map(len, sample[0])
        resampler = [rs.randint(0, n, n) for n in lengths]
        sample = [[c.take(r, axis=0) for c, r in zip(a, resampler)]
                  for a in sample]
        sample = list(map(np.concatenate, sample))
        boot_dist.append(func(*sample, **func_kwargs))
    return np.array(boot_dist)


def _smooth_bootstrap(args, n_boot, func, func_kwargs):
    """Bootstrap by resampling from a kernel density estimate."""
    n = len(args[0])
    boot_dist = []
    kde = [stats.gaussian_kde(np.transpose(a)) for a in args]
    for i in range(int(n_boot)):
        sample = [a.resample(n).T for a in kde]
        boot_dist.append(func(*sample, **func_kwargs))
    return np.array(boot_dist)


def randomize_corrmat(a, tail="both", corrected=True, n_iter=1000,
                      random_seed=None, return_dist=False):
    """Test the significance of set of correlations with permutations.

    By default this corrects for multiple comparisons across one side
    of the matrix.

    Parameters
    ----------
    a : n_vars x n_obs array
        array with variables as rows
    tail : both | upper | lower
        whether test should be two-tailed, or which tail to integrate over
    corrected : boolean
        if True reports p values with respect to the max stat distribution
    n_iter : int
        number of permutation iterations
    random_seed : int or None
        seed for RNG
    return_dist : bool
        if True, return n_vars x n_vars x n_iter

    Returns
    -------
    p_mat : float
        array of probabilites for actual correlation from null CDF

    """
    if tail not in ["upper", "lower", "both"]:
        raise ValueError("'tail' must be 'upper', 'lower', or 'both'")

    rs = np.random.RandomState(random_seed)

    a = np.asarray(a, np.float)
    flat_a = a.ravel()
    n_vars, n_obs = a.shape

    # Do the permutations to establish a null distribution
    null_dist = np.empty((n_vars, n_vars, n_iter))
    for i_i in range(n_iter):
        perm_i = np.concatenate([rs.permutation(n_obs) + (v * n_obs)
                                 for v in range(n_vars)])
        a_i = flat_a[perm_i].reshape(n_vars, n_obs)
        null_dist[..., i_i] = np.corrcoef(a_i)

    # Get the observed correlation values
    real_corr = np.corrcoef(a)

    # Figure out p values based on the permutation distribution
    p_mat = np.zeros((n_vars, n_vars))
    upper_tri = np.triu_indices(n_vars, 1)

    if corrected:
        if tail == "both":
            max_dist = np.abs(null_dist[upper_tri]).max(axis=0)
        elif tail == "lower":
            max_dist = null_dist[upper_tri].min(axis=0)
        elif tail == "upper":
            max_dist = null_dist[upper_tri].max(axis=0)

        cdf = lambda x: stats.percentileofscore(max_dist, x) / 100.

        for i, j in zip(*upper_tri):
            observed = real_corr[i, j]
            if tail == "both":
                p_ij = 1 - cdf(abs(observed))
            elif tail == "lower":
                p_ij = cdf(observed)
            elif tail == "upper":
                p_ij = 1 - cdf(observed)
            p_mat[i, j] = p_ij

    else:
        for i, j in zip(*upper_tri):

            null_corrs = null_dist[i, j]
            cdf = lambda x: stats.percentileofscore(null_corrs, x) / 100.

            observed = real_corr[i, j]
            if tail == "both":
                p_ij = 2 * (1 - cdf(abs(observed)))
            elif tail == "lower":
                p_ij = cdf(observed)
            elif tail == "upper":
                p_ij = 1 - cdf(observed)
            p_mat[i, j] = p_ij

    # Make p matrix symettrical with nans on the diagonal
    p_mat += p_mat.T
    p_mat[np.diag_indices(n_vars)] = np.nan

    if return_dist:
        return p_mat, null_dist
    return p_mat
