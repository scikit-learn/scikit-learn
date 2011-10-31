"""
Robust location and covariance estimators.

Here are implemented estimators that are resistant to outliers.

"""
# Author: Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD Style.
import warnings
import numpy as np
from scipy import linalg
from scipy.stats import chi2
from . import empirical_covariance, EmpiricalCovariance
from ..utils.extmath import fast_logdet


###############################################################################
### Minimum Covariance Determinant
#   Implementing of an algorithm by Rousseeuw & Van Driessen described in
#   (A Fast Algorithm for the Minimum Covariance Determinant Estimator,
#   1999, American Statistical Association and the American Society
#   for Quality, TECHNOMETRICS)
###############################################################################
def c_step(X, h, remaining_iterations=30, initial_estimates=None,
           verbose=False, cov_computation_method=empirical_covariance):
    """C_step procedure described in [1] aiming at computing the MCD

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      Data set in which we look for the h observations whose scatter matrix
      has minimum determinant
    h: int, > n_samples / 2
      Number of observations to compute the ribust estimates of location
      and covariance from.
    remaining_iterations: int
      Number of iterations to perform.
      According to Rousseeuw [1], two iterations are sufficient to get close
      to the minimum, and we never need more than 30 to reach convergence.
    initial_estimates: 2-tuple
      Initial estimates of location and shape from which to run the c_step
      procedure:
      - initial_estimates[0]: an initial location estimate
      - initial_estimates[1]: an initial covariance estimate
    verbose: boolean
      Verbose mode

    Returns
    -------
    location: array-like, shape (n_features,)
      Robust location estimates
    covariance: array-like, shape (n_features, n_features)
      Robust covariance estimates
    support: array-like, shape (n_samples,)
      A mask for the `h` observations whose scatter matrix has minimum
      determinant

    Notes
    -----
    References:
    [1] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
        1999, American Statistical Association and the American Society
        for Quality, TECHNOMETRICS

    """
    n_samples, n_features = X.shape

    # Initialisation
    if initial_estimates is None:
        # compute initial robust estimates from a random subset
        support = np.zeros(n_samples).astype(bool)
        support[np.random.permutation(n_samples)[:h]] = True
        location = X[support].mean(0)
        covariance = cov_computation_method(X[support])
    else:
        # get initial robust estimates from the function parameters
        location = initial_estimates[0]
        covariance = initial_estimates[1]
        # run a special iteration for that case (to get an initial support)
        precision = linalg.pinv(covariance)
        X_centered = X - location
        dist = (np.dot(X_centered, precision) * X_centered).sum(1)
        # compute new estimates
        support = np.zeros(n_samples).astype(bool)
        support[np.argsort(dist)[:h]] = True
        location = X[support].mean(0)
        covariance = cov_computation_method(X[support])
        det = fast_logdet(covariance)
    previous_det = np.inf

    # Iterative procedure for Minimum Covariance Determinant computation
    det = fast_logdet(covariance)
    while (det < previous_det) and (remaining_iterations > 0):
        # compute a new support from the full data set mahalanobis distances
        precision = linalg.pinv(covariance)
        X_centered = X - location
        dist = (np.dot(X_centered, precision) * X_centered).sum(1)
        # save old estimates values
        previous_location = location
        previous_covariance = covariance
        previous_det = det
        previous_support = support
        # compute new estimates
        support = np.zeros(n_samples).astype(bool)
        support[np.argsort(dist)[:h]] = True
        location = X[support].mean(0)
        covariance = cov_computation_method(X[support])
        det = np.log(linalg.det(covariance))
        # update remaining iterations for early stopping
        remaining_iterations = remaining_iterations - 1

    # Check convergence
    if np.allclose(det, previous_det):
        # c_step procedure converged
        if verbose:
            print "Optimal couple (location, covariance) found before" \
                "ending iterations (%d left)" % (remaining_iterations)
        results = location, covariance, det, support
    elif det > previous_det:
        # determinant has increased (should not happen)
        warnings.warn("Warning! det > previous_det (%.15f > %.15f)" \
                          % (det, previous_det), RuntimeWarning)
        results = previous_location, previous_covariance, \
            previous_det, previous_support

    # Check early stopping
    if remaining_iterations == 0:
        if verbose:
            print 'Maximum number of iterations reached'
        det = fast_logdet(covariance)
        results = location, covariance, det, support

    return results


def select_candidates(X, h, n_trials, select=1, n_iter=30, verbose=False,
                      cov_computation_method=empirical_covariance):
    """Finds the best pure subset of observations to compute MCD from it.

    The purpose of this function is to find the best sets of h
    observations with respect to a minimization of their covariance
    matrix determinant. Equivalently, it removes n_samples-h
    observations to construct what we call a pure data set (i.e. not
    containing outliers). The list of the observations of the pure
    data set is referred to as the `support`.

    Starting from a random support, the pure data set is found by the
    c_step procedure introduced by Rousseeuw and Van Driessen in [1].

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      Data (sub)set in which we look for the h purest observations
    h: int, [(n + p + 1)/2] < h < n
      The number of samples the pure data set must contain.
    select: int, int > 0
      Number of best candidates results to return.
    n_trials: int, nb_trials > 0 or 2-tuple
      Number of different initial sets of observations from which to
      run the algorithm.
      Instead of giving a number of trials to perform, one can provide a
      list of initial estimates that will be used to iteratively run
      c_step procedures. In this case:
      - n_trials[0]: array-like, shape (n_trials, n_features)
        is the list of `n_trials` initial location estimates
      - n_trials[1]: array-like, shape (n_trials, n_features, n_features)
        is the list of `n_trials` initial covariances estimates
    n_iter: int, nb_iter > 0
      Maximum number of iterations for the c_step procedure.
      (2 is enough to be close to the final solution. "Never" exceeds 20)

    See
    ---
    `c_step` function

    Returns
    -------
    best_locations: array-like, shape (select, n_features)
      The `select` location estimates computed from the `select` best
      supports found in the data set (`X`)
    best_covariances: array-like, shape (select, n_features, n_features)
      The `select` covariance estimates computed from the `select`
      best supports found in the data set (`X`)
    best_supports: array-like, shape (select, n_samples)
      The `select` best supports found in the data set (`X`)

    Notes
    -----
    References:
    [1] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
        1999, American Statistical Association and the American Society
        for Quality, TECHNOMETRICS

    """
    n_samples, n_features = X.shape

    if isinstance(n_trials, int):
        run_from_estimates = False
    elif isinstance(n_trials, tuple):
        run_from_estimates = True
        estimates_list = n_trials
        n_trials = estimates_list[0].shape[0]
    else:
        raise Exception('Bad \'n_trials\' parameter (wrong type)')

    # compute `n_trials` location and shape estimates candidates in the subset
    all_estimates = []
    if not run_from_estimates:
        # perform `n_trials` computations from random initial supports
        for j in range(n_trials):
            all_estimates.append(
                c_step(
                    X, h, remaining_iterations=n_iter, verbose=verbose,
                    cov_computation_method=cov_computation_method))
    else:
        # perform computations from every given initial estimates
        for j in range(n_trials):
            initial_estimates = (estimates_list[0][j], estimates_list[1][j])
            all_estimates.append(c_step(
                    X, h, remaining_iterations=n_iter,
                    initial_estimates=initial_estimates, verbose=verbose,
                    cov_computation_method=cov_computation_method))
    all_locations_sub, all_covariances_sub, all_dets_sub, all_supports_sub = \
                                     zip(*all_estimates)
    # find the `n_best` best results among the `n_trials` ones
    index_best = np.argsort(all_dets_sub)[:select]
    best_locations = np.asarray(all_locations_sub)[index_best]
    best_covariances = np.asarray(all_covariances_sub)[index_best]
    best_supports = np.asarray(all_supports_sub)[index_best]

    return best_locations, best_covariances, best_supports


def fast_mcd(X, h=None, correction="empirical", reweighting="rousseeuw",
             cov_computation_method=empirical_covariance):
    """Estimates the Minimum Covariance Determinant matrix.

    The FastMCD algorithm has been introduced by Rousseuw and Van Driessen
    in "A Fast Algorithm for the Minimum Covariance Determinant Estimator,
    1999, American Statistical Association and the American Society
    for Quality, TECHNOMETRICS".
    The principle is to compute robust estimates and random subsets before
    pooling them into a larger subsets, and finally into the full data set.
    Depending on the size of the initial sample, we have one, two or three
    such computation levels.

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      The data matrix, with p features and n samples.
    h: float, 0 < h < 1
          The proportion of points to be included in the support of the raw
          MCD estimate. Default is None, which implies that the minimum
          value of h will be used within the algorithm:
          [n_sample + n_features + 1] / 2
    correction: str
      Improve the covariance estimator consistency at gaussian models
        - "empirical" (default): correction using the empirical correction
          factor suggested by Rousseeuw and Van Driessen in [1]
        - "theoretical": correction using the theoretical correction factor
          derived in [2]
        - else: no correction
    reweighting: str
      Computation of a reweighted estimator:
        - "rousseeuw" (default): Reweight observations using Rousseeuw's
          method (equivalent to deleting outlying observations from the
          data set before computing location and covariance estimates)
        - else: no reweighting

    Notes
    -----
    References:
    [1] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
        1999, American Statistical Association and the American Society
        for Quality, TECHNOMETRICS
    [2] R. W. Butler, P. L. Davies and M. Jhun,
        Asymptotics For The Minimum Covariance Determinant Estimator,
        The Annals of Statistics, 1993, Vol. 21, No. 3, 1385-1400

    Returns
    -------
    location: array-like, shape (n_features,)
      Robust location of the data
    covariance: array-like, shape (n_features, n_features)
      Robust covariance of the features
    support: array-like, type boolean, shape (n_samples,)
      a mask of the observations that have been used to compute
      the robust location and covariance estimates of the data set

    """
    X = np.asanyarray(X)
    if X.ndim == 1:
        X = np.reshape(X, (1, -1))
        warnings.warn("Only one sample available. " \
                          "You may want to reshape your data array")
        n_samples = 1
        n_features = X.size
    else:
        n_samples, n_features = X.shape

    # minimum breakdown value
    if h is None:
        h = int(np.ceil(0.5 * (n_samples + n_features + 1)))
    else:
        h = int(h * n_samples)

    # 1-dimensional case quick computation
    # (Rousseeuw, P. J. and Leroy, A. M. (2005) References, in Robust
    #  Regression and Outlier Detection, John Wiley & Sons, chapter 4)
    if n_features == 1:
        # find the sample shortest halves
        X_sorted = np.sort(np.ravel(X))
        diff = X_sorted[h:] - X_sorted[:(n_samples - h)]
        halves_start = np.where(diff == np.min(diff))[0]
        # take the middle points' mean to get the robust location estimate
        location = 0.5 * \
            (X_sorted[h + halves_start] + X_sorted[halves_start]).mean()
        support = np.zeros(n_samples).astype(bool)
        support[np.argsort(np.abs(X - location), axis=0)[:h]] = True
        covariance = np.asarray([[np.var(X[support])]])
        location = np.array([location])

    ### Starting FastMCD algorithm for p-dimensional case
    if (n_samples > 500) and (n_features > 1):
        ## 1. Find candidate supports on subsets
        # a. split the set in subsets of size ~ 300
        n_subsets = n_samples / 300
        n_samples_subsets = n_samples / n_subsets
        samples_shuffle = np.random.permutation(n_samples)
        h_subset = np.ceil(n_samples_subsets * (h / float(n_samples)))
        # b. perform a total of 500 trials
        n_trials_tot = 500
        n_trials = n_trials_tot / n_subsets
        # c. select 10 best (location, covariance) for each subset
        n_best_sub = 10
        n_best_tot = n_subsets * n_best_sub
        all_best_locations = np.zeros((n_best_tot, n_features))
        all_best_covariances = np.zeros((n_best_tot, n_features, n_features))
        for i in range(n_subsets):
            low_bound = i * n_samples_subsets
            high_bound = low_bound + n_samples_subsets
            current_subset = X[samples_shuffle[low_bound:high_bound]]
            best_locations_sub, best_covariances_sub, _ = select_candidates(
                current_subset, h_subset, n_trials,
                select=n_best_sub, n_iter=2,
                cov_computation_method=cov_computation_method)
            subset_slice = np.arange(i * n_best_sub, (i + 1) * n_best_sub)
            all_best_locations[subset_slice] = best_locations_sub
            all_best_covariances[subset_slice] = best_covariances_sub
        ## 2. Pool the candidate supports into a merged set
        ##    (possibly the full dataset)
        n_samples_merged = min(1500, n_samples)
        h_merged = np.ceil(n_samples_merged * (h / float(n_samples)))
        if n_samples > 1500:
            n_best_merged = 10
        else:
            n_best_merged = 1
        # find the best couples (location, covariance) on the merged set
        locations_merged, covariances_merged, supports_merged = \
            select_candidates(
            X[np.random.permutation(n_samples)[:n_samples_merged]],
            h_merged, n_trials=(all_best_locations, all_best_covariances),
            select=n_best_merged,
            cov_computation_method=cov_computation_method)
        ## 3. Finally get the overall best (locations, covariance) couple
        if n_samples < 1500:
            # directly get the best couple (location, covariance)
            location = locations_merged[0]
            covariance = covariances_merged[0]
            support = supports_merged[0]
        else:
            # select the best couple on the full dataset
            locations_full, covariances_full, supports_full = \
                select_candidates(
                X, h, n_trials=(locations_merged, covariances_merged),
                select=1,
                cov_computation_method=cov_computation_method)
            location = locations_full[0]
            covariance = covariances_full[0]
            support = supports_full[0]
    elif n_features > 1:
        ## 1. Find the 10 best couples (location, covariance)
        ## considering two iterations
        n_trials = 30
        n_best = 10
        locations_best, covariances_best, _ = select_candidates(
            X, h, n_trials=n_trials, select=n_best, n_iter=2,
            cov_computation_method=cov_computation_method)
        ## 2. Select the best couple on the full dataset amongst the 10
        locations_full, covariances_full, supports_full = select_candidates(
            X, h, n_trials=(locations_best, covariances_best), select=1,
            cov_computation_method=cov_computation_method)
        location = locations_full[0]
        covariance = covariances_full[0]
        support = supports_full[0]

    ## 4. Obtain consistency at Gaussian models
    covariance_corrected = _correct(
        location, covariance, support, correction, data=X)

    ## 5. Reweight estimate
    location_reweighted, covariance_reweighted, support = _reweight(
        location, covariance_corrected, support, reweighting, data=X,
        cov_computation_method=cov_computation_method)

    return location_reweighted, covariance_reweighted, support


def _correct(location, covariance, support, correction="empirical", data=None):
    """Apply a correction to a raw Minimum Covariance Determinant estimate

    Parameters
    ----------
    location: array-like, shape (n_features,)
      Raw robust location estimate.
    covariance: array-like, shape (n_features, n_features)
      Raw robust covariance estimate.
    support: array-like, type boolean, shape (n_samples,)
      A mask of the observations that have been used to compute
      the raw robust location and covariance estimates.
    correction: str
      Improve the covariance estimator consistency at gaussian models
        - "empirical" (default): correction using the empirical correction
          factor suggested by Rousseeuw and Van Driessen in [1]
        - "theoretical": correction using the theoretical correction factor
          derived in [2]
        - else: no correction
    data: array-like, shape (n_samples, n_features)
      The data matrix, with p features and n samples.
      The data set must be the one which was used to compute the raw estimates.

    Returns
    -------
    covariance_corrected: array-like, shape (n_features, n_features)
      Corrected robust covariance estimate.

    """
    h = support.sum()
    n_samples = support.size
    n_features = covariance.shape[0]
    if correction == "theoretical":
        # theoretical correction
        inliers_ratio = h / float(n_samples)
        covariance_corrected = covariance * (inliers_ratio) \
            / chi2(n_features + 2).cdf(chi2(n_features).ppf(inliers_ratio))
    elif correction == "empirical":
        # empirical correction
        if data is None:
            raise ValueError("Need `data` for empirical correction")
        else:
            X_centered = data - location
            dist = np.sum(
                    np.dot(X_centered, linalg.pinv(covariance)) * X_centered,
                    1)
            covariance_corrected = covariance * \
                              (np.median(dist) / chi2(n_features).isf(0.5))
    else:
        # no correction
        covariance_corrected = covariance

    return covariance_corrected


def _reweight(location, covariance, support,
              reweighting="rousseeuw", data=None,
              cov_computation_method=empirical_covariance):
    """Reweight a raw Minimum Covariance Determinant estimate

    Parameters
    ----------
    location: array-like, shape (n_features,)
      Raw robust location estimate.
    covariance: array-like, shape (n_features, n_features)
      Raw robust covariance estimate.
    support: array-like, type boolean, shape (n_samples,)
      A mask of the observations that have been used to compute
      the raw robust location and covariance estimates.
    reweighting: str
      Computation of a reweighted estimator:
        - "rousseeuw" (default): Reweight observations using Rousseeuw's
          method (equivalent to deleting outlying observations from the
          data set before computing location and covariance estimates)
        - else: no reweighting
    data: array-like, shape (n_samples, n_features)
      The data matrix, with p features and n samples.
      The data set must be the one which was used to compute the raw estimates.

    Returns
    -------
    location_reweighted: array-like, shape (n_features, )
      Reweighted robust location estimate.
    covariance_reweighted: array-like, shape (n_features, n_features)
      Reweighted robust covariance estimate.
    support: array-like, type boolean, shape (n_samples,)
      A mask of the observations that have been used to compute
      the reweighted robust location and covariance estimates.

    """
    n_features = covariance.shape[0]
    if reweighting == "rousseeuw":
        if data is None:
            raise ValueError("Need data for Rousseeuw reweighting")
        else:
            n_samples = data.shape[0]
            X_centered = data - location
            dist = np.sum(
                np.dot(X_centered, linalg.pinv(covariance)) * X_centered,
                1)
            mask = dist < chi2(n_features).isf(0.025)
            location_reweighted = data[mask].mean(0)
            covariance_reweighted = cov_computation_method(data[mask])
            support = np.zeros(n_samples).astype(bool)
            support[mask] = True
    else:
        location_reweighted = location
        covariance_reweighted = covariance
        support = support

    return location_reweighted, covariance_reweighted, support


class MinCovDet(EmpiricalCovariance):
    """Minimum Covariance Determinant (MCD) robust estimator of covariance

    The Minimum Covariance Determinant estimator is a robust estimator
    of a data set's covariance introduced by P.J.Rousseuw in [1].
    The idea is to find a given proportion (h) of "good" observations which
    are not outliers and compute their empirical covariance matrix.
    This empirical covariance matrix is then rescaled to compensate the
    performed selection of observations ("consistency step").
    Having computed the Minimum Covariance Determinant estimator, one
    can give weights to observations according to their Mahalanobis
    distance, leading the a reweighted estimate of the covariance
    matrix of the data set.

    Rousseuw and Van Driessen [2] developed the FastMCD algorithm in order
    to compute the Minimum Covariance Determinant. This algorithm is used
    when fitting an MCD object to data.
    The FastMCD algorithm also computes a robust estimate of the data set
    location at the same time.

    Attributes
    ----------
    `location_`: array-like, shape (n_features,)
        Estimated robust location

    `covariance_`: array-like, shape (n_features, n_features)
        Estimated robust covariance matrix

    `precision_`: array-like, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)

    `support_`: array-like, shape (n_samples,)
        A mask of the observations that have been used to compute
        the robust estimates of location and shape.

    Notes
    -----
    References:
    [1] P. J. Rousseeuw. Least median of squares regression. J. Am
        Stat Ass, 79:871, 1984.
    [2] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
        1999, American Statistical Association and the American Society
        for Quality, TECHNOMETRICS
    [3] R. W. Butler, P. L. Davies and M. Jhun,
        Asymptotics For The Minimum Covariance Determinant Estimator,
        The Annals of Statistics, 1993, Vol. 21, No. 3, 1385-1400

    """
    def __init__(self, store_precision=True, assume_centered=False,
                 h=None, correction="empirical", reweighting="rousseeuw"):
        """

        Parameters
        ----------
        store_precision: bool
          Specify if the estimated precision is stored
        assume_centered: Boolean
          If True, the support of robust location and covariance estimates
          is computed, and a covariance estimate is recomputed from it,
          without centering the data.
          Useful to work with data whose mean is significantly equal to
          zero but is not exactly zero.
          If False, the robust location and covariance are directly computed
          with the FastMCD algorithm without additional treatment.
        h: float, 0 < h < 1
          The proportion of points to be included in the support of the raw
          MCD estimate. Default is None, which implies that the minimum
          value of h will be used within the algorithm:
          [n_sample + n_features + 1] / 2
        correction: str
          Improve the covariance estimator consistency at Gaussian models
            - "empirical" (default): correction using the empirical correction
              factor suggested by Rousseeuw and Van Driessen in [1]
            - "theoretical": correction using the theoretical correction factor
              derived in [2]
            - else: no correction
        reweighting: str
          Computation of a reweighted estimator:
            - "rousseeuw" (default): Reweight observations using Rousseeuw's
              method (equivalent to deleting outlying observations from the
              data set before computing location and covariance estimates)
            - else: no re-weighting
        """
        self.store_precision = store_precision
        self.assume_centered = assume_centered
        self.h = h
        self.correction = correction
        self.reweighting = reweighting

    def fit(self, X):
        """Fits a Minimum Covariance Determinant with the FastMCD algorithm.

        Parameters
        ----------
        X: array-like, shape = [n_samples, n_features]
          Training data, where n_samples is the number of samples
          and n_features is the number of features.

        Returns
        -------
        self: object
          Returns self.

        """
        n_samples, n_features = X.shape
        # compute and store raw estimates
        raw_location, raw_covariance, raw_support = fast_mcd(
                X, h=self.h, correction=None, reweighting=None,
                cov_computation_method=self._nonrobust_covariance)
        if self.assume_centered:
            raw_location = np.zeros(n_features)
            raw_covariance = self._nonrobust_covariance(
                    X[raw_support], assume_centered=True)
        if self.h is None:
            self.h = int(np.ceil(0.5 * (n_samples + n_features + 1)))
        self.raw_location_ = raw_location
        self.raw_covariance_ = raw_covariance
        self.raw_support_ = raw_support
        self.location_ = raw_location
        self.support_ = raw_support
        # obtain consistency at normal models
        covariance = self.correct(self.correction, X, permanent=True)
        # reweight estimator
        location, covariance, support = self.reweight(
                self.reweighting, X, permanent=True)

        return self

    def _nonrobust_covariance(self, data, assume_centered=False):
        """Non-robust estimation of the covariance to be used within MCD.

        Parameters
        ----------
        data: array_like, shape (n_samples, n_features)
          Data for which to compute the non-robust covariance matrix.
        assume_centered: Boolean
          Whether or not the observations should be considered as centered.

        Returns
        -------
        nonrobust_covariance: array_like, shape (n_features, n_features)
          The non-robust covariance of the data.

        """
        return empirical_covariance(data, assume_centered=assume_centered)

    def correct(self, correction="empirical", data=None, permanent=False):
        """Apply a correction to raw Minimum Covariance Determinant estimates.

        Parameters
        ----------
        correction: str
          Improve the covariance estimator consistency at gaussian models
            - "empirical" (default): correction using the empirical correction
              factor suggested by Rousseeuw and Van Driessen in [1]
            - "theoretical": correction using the theoretical correction factor
              derived in [2]
            - else: no correction
        data: array-like, shape (n_samples, n_features)
          The data matrix, with p features and n samples.
          The data set must be the one which was used to compute
          the raw estimates.
        permanent: bool
          If True, the `reweighting` parameter is taken a the new object's
          `reweighting` parameter and the reweighting is persistent.

        Returns
        -------
        covariance: array-like, shape (n_features, n_features)
          Corrected robust covariance estimate.

        Notes
        -----
        References:
        [1] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
            1999, American Statistical Association and the American Society
            for Quality, TECHNOMETRICS
        [2] R. W. Butler, P. L. Davies and M. Jhun,
            Asymptotics For The Minimum Covariance Determinant Estimator,
            The Annals of Statistics, 1993, Vol. 21, No. 3, 1385-1400

        """
        covariance = _correct(self.raw_location_, self.raw_covariance_,
                     self.raw_support_, correction, data=data)
        if permanent:
            self._set_estimates(covariance)
            self.correction = correction
        return covariance

    def reweight(self, reweighting="rousseeuw", data=None, permanent=False):
        """Reweight raw Minimum Covariance Determinant estimates.

        Parameters
        ----------
        reweighting: str
          Computation of a reweighted estimator:
            - "rousseeuw" (default): Reweight observations using Rousseeuw's
              method (equivalent to deleting outlying observations from the
              data set before computing location and covariance estimates)
            - else: no re-weighting
        data: array-like, shape (n_samples, n_features)
          The data matrix, with p features and n samples.
          The data set must be the one which was used to compute
          the raw estimates.
        permanent: bool
          If True, the `reweighting` parameter is taken a the new object's
          `reweighting` parameter and the reweighting is persistent.

        Returns
        -------
        location: array-like, shape (n_features, )
          Reweighted robust location estimate.
        covariance: array-like, shape (n_features, n_features)
          Reweighted robust covariance estimate.
        support: array-like, type boolean, shape (n_samples,)
          A mask of the observations that have been used to compute
          the reweighted robust location and covariance estimates.

        """
        location, covariance, support = _reweight(
                self.location_, self.covariance_, self.support_,
                reweighting, data=data,
                cov_computation_method=self._nonrobust_covariance)
        if permanent:
            self._set_estimates(covariance)
            self.location_ = location
            self.support_ = support
            self.reweighting = reweighting
        return location, covariance, support
