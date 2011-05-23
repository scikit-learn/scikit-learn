"""


"""
# Author: Virgile Fritsch <virgile.fritsch@inria.fr>
#         implementing a method by Rousseeuw & Van Driessen described in
#         (A Fast Algorithm for the Minimum Covariance Determinant Estimator,
#         1999, American Statistical Association and the American Society
#         for Quality, TECHNOMETRICS)
#
# License: BSD Style.
import numpy as np
from scipy import linalg
from scipy.stats import chi2
from scikits.learn.covariance import empirical_covariance
from ..utils.extmath import fast_logdet as exact_logdet


def c_step(X, h, remaining_iterations=30, initial_estimates=None,
           initial_support=None, verbose=False):
    """

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      ?
    support: ?
      ?
    remaining_iterations: int
      ?
    initial_estimates: ?
      ?
    verbose: boolean
      ?

    Returns
    -------
    ?
    
    """
    # TODO: function parameters checking
    n_samples, n_features = X.shape
    
    # Initialisation
    if initial_estimates is None:
        # compute initial robust location and covariance estimates
        if initial_support is None:
            support_aux = np.zeros(n_samples).astype(bool)
            support_aux[np.random.permutation(n_samples)[:h]] = True
            support = support_aux
        else:
            if support.dtype == 'bool':
                support = initial_support
            else:
                support_aux = np.zeros(n_samples).astype(bool)
                support_aux[initial_support] = True
                support = support_aux
            h = support.size
        T = X[support].mean(0)
        S = empirical_covariance(X[support])
    else:
        # get estimates from function parameters
        T = initial_estimates[0]
        S = initial_estimates[1]
        # run a special iteration for that case
        # (in order to get a meaningful support and to avoid the algorithm
        #  to be traped in a local minimum)
        inv_S = linalg.pinv(S)
        X_centered = X - T
        dist = (np.dot(X_centered,inv_S) * X_centered).sum(1)
        # compute new estimates
        support_aux = np.zeros(n_samples).astype(bool)
        support_aux[np.argsort(dist)[:h]] = True
        support = support_aux
        T = X[support].mean(0)
        S = empirical_covariance(X[support])
        detS = exact_logdet(S)
    previous_detS = np.inf
    
    # Iterative procedure for Minimum Covariance Determinant computation
    detS = exact_logdet(S)
    while (detS < previous_detS) and (remaining_iterations > 0):
        # compute a new support from the full data set mahalanobis distances
        inv_S = linalg.pinv(S)
        X_centered = X - T
        dist = (np.dot(X_centered,inv_S) * X_centered).sum(1)
        # save old estimates values
        previous_T = T
        previous_S = S
        previous_detS = detS
        previous_support = support
        # compute new estimates
        support_aux = np.zeros(n_samples).astype(bool)
        support_aux[np.argsort(dist)[:h]] = True
        support = support_aux
        T = X[support].mean(0)
        S = empirical_covariance(X[support])
        detS = np.log(linalg.det(S))
        remaining_iterations = remaining_iterations - 1
    
    # Check convergence
    if np.allclose(detS, previous_detS):
        # c_step procedure converged
        if verbose:
            print 'Optimal couple (T,S) found before ending iterations' \
                '(%d left)' %(remaining_iterations)
        results = T, S, detS, support
    elif detS > previous_detS:
        # determinant has increased (should not happen)
        # FIXME: Transform in real warning
        print "Warning! detS > previous_detS (%f > %f)" %(detS, previous_detS)
        1/0
        results = previous_T, previous_S, previous_detS, previous_support
    
    # Check early stopping
    if remaining_iterations == 0:
        if verbose:
            print 'Maximum number of iterations reached'
        detS = exact_logdet(S)
        results = T, S, detS, support

    return results


def select_candidates(X, h, n_trials, select=1, n_iter=30, verbose=False):
    """Finds the best pure subset of observations to compute MCD from it.
    
    The purpose of this function is to find the best sets of h
    observations with respect to a minimization of their covariance
    matrix determinant. Equivalently, it removes n_samples-h
    observations to construct what we call a pure data set (i.e. not
    containing outliers). The list of the observations of the pure
    data set is referred to as the `support`.
    
    Starting from a random support, the pure data set is found by the
    c_step procedure introduced by Rousseeuw and Van Driessen in [1].
    
    [1] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
        1999, American Statistical Association and the American Society
        for Quality, TECHNOMETRICS

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      Data set in which we look for the h purest observations
    h: int, [(n + p + 1)/2] < h < n
      The number of samples the pure data set must contain.
    select: int, int > 0
      Number of best candidates results to return.
    n_trials: int, nb_trials > 0
      Number of different initial sets of observations from which to
      run the algorithm.
    n_iter: int, nb_iter > 0
      Maximum number of iterations for the c_step procedure.
      (2 is enough to be close to the final solution. "Never" exceeds 20)

    Returns
    -------
    best_T: array-like, shape (select, n_features)
      The `select` location estimates computed from the `select` best
      supports found in the data set (`X`)
    best_S: array-like, shape (select, n_features, n_features)
      The `select` covariance estimates computed from the `select`
      best supports found in the data set (`X`)
    best_H: array-like, shape (select, n_samples)
      The `select` best supports found in the data set (`X`)

    """
    X = np.asanyarray(X)
    if X.ndim <= 1:
        X = X.reshape((-1,1))
    n_samples, n_features = X.shape
    
    if isinstance(n_trials, int):
        run_from_estimates = False
    elif isinstance(n_trials, tuple):
        run_from_estimates = True
        estimates_list = n_trials
        n_trials = estimates_list[0][0].shape[0]
    else:
        raise Exception('Bad \'n_trials\' parameter (wrong type)')
    
    # find `select` best location and shape estimates candidates in the subset
    all_T_sub = np.zeros((n_trials, n_features))
    all_S_sub = np.zeros((n_trials, n_features, n_features))
    all_detS_sub = np.zeros(n_trials)
    all_H_sub = np.zeros((n_trials, n_samples)).astype(bool)
    # perform `n_trials` randomly seeded trials on the current subset
    if not run_from_estimates:
        for j in range(n_trials):
            all_T_sub[j], all_S_sub[j], all_detS_sub[j], all_H_sub[j] = c_step(
                X, h, remaining_iterations=n_iter, verbose=verbose)
    else:
        for j in range(n_trials):
            all_T_sub[j], all_S_sub[j], all_detS_sub[j], all_H_sub[j] = c_step(
                X, h, remaining_iterations=n_iter,
                initial_estimates=(estimates_list[0][j], estimates_list[1][j]))
    
    # find the `n_best` best results among the `n_trials` ones
    mask_best = np.argsort(all_detS_sub)[:select]
    best_T = all_T_sub[mask_best]
    best_S = all_S_sub[mask_best]
    best_H = all_H_sub[mask_best]
    
    return best_T, best_S, best_H


def fast_mcd(X, correction="empirical", reweight="rousseeuw"):
    """Estimates the Minimum Covariance Determinant matrix.
    
    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      The data matrix, with p features and n samples.
    reweight: int (0 < reweight < 2)
      Computation of a reweighted estimator:
        - "rousseeuw" (default): Reweight observations using Rousseeuw's
          method (equivalent to deleting outlying observations from the
          data set before computing location and covariance estimates)
        - else: no reweighting
    correction : str
      Improve the covariance estimator consistency at gaussian models
        - "empirical" (default): correction using the empirical correction
          factor derived in [1]
        - "theoretical": correction using the theoretical correction factor
          derived in [2]
        - else: no correction

    Returns
    -------
    T: array-like, shape (n_features,)
      Robust location of the data
    S: array-like, shape (n_features, n_features)
      Robust covariance of the features
    support: array-like, type boolean, shape (n_samples,)
      a mask of the observations that have been used to compute
      the robust location and covariance estimates of the data set
    
    """
    X = np.asanyarray(X)
    if X.ndim <= 1:
        X = X.reshape((-1,1))
    n_samples, n_features = X.shape

    # minimum breakdown value
    h = np.ceil(0.5 * (n_samples + n_features + 1))

    ### Starting FastMCD algorithm
    if n_samples > 500:
        ## 1. Find candidate supports on subsets
        # a. split the set in subsets of size ~ 300
        n_subsets = n_samples / 300
        n_samples_subsets = n_samples / n_subsets
        samples_shuffle = np.random.permutation(n_samples)
        h_subset = np.ceil(n_samples_subsets * (h/float(n_samples)))
        # b. perform a total of 500 trials
        n_trials_tot = 500
        n_trials = n_trials_tot / n_subsets
        # c. select 10 best (T,S) for each subset
        n_best_sub = 10
        n_best_tot = n_subsets*n_best_sub
        all_best_T = np.zeros((n_best_tot, n_features))
        all_best_S = np.zeros((n_best_tot, n_features, n_features))
        for i in range(n_subsets):
            low_bound = i*n_samples_subsets
            high_bound = low_bound + n_samples_subsets
            current_subset = X[samples_shuffle[low_bound:high_bound]]
            best_T_sub, best_S_sub, _ = select_candidates(
                current_subset, h_subset, n_trials, select=n_best_sub, n_iter=2)
            subset_slice = np.arange(i*n_best_sub, (i+1)*n_best_sub)
            all_best_T[subset_slice] = best_T_sub
            all_best_S[subset_slice] = best_S_sub
        ## 2. Pool the candidate supports into a merged set
        ##    (possibly the full dataset)
        n_samples_merged = min(1500, n_samples)
        h_merged = np.ceil(n_samples_merged * (h/float(n_samples)))
        if n_samples > 1500:
            n_best_merged = 10
        else:
            n_best_merged = 1
        # find the best couples (T,S) on the merged set
        T_merged, S_merged, H_merged = select_candidates(
            X[np.random.permutation(n_samples)[:n_samples_merged]],
            h_merged, n_trials=(all_best_T, all_best_S), select=n_best_merged)
        ## 3. Finally get the overall best (T,S) couple
        if n_samples < 1500:
            # directly get the best couple (T,S)
            T = T_merged[0]
            S = S_merged[0]
            H = H_merged[0]
        else:
            # select the best couple on the full dataset
            T_full, S_full, H_full = select_candidates(
                X, h, n_trials=(T_merged, S_merged), select=1)
            T = T_full[0]
            S = S_full[0]
            H = H_full[0]
    else:
        ## 1. Find the 10 best couples (T,S) considering two iterations
        n_trials = 30
        n_best = 10
        T_best, S_best, _ = select_candidates(
            X, h, n_trials=n_trials, select=n_best, n_iter=2)
        ## 2. Select the best couple on the full dataset amongst the 10
        T_full, S_full, H_full = select_candidates(
            X, h, n_trials=(T_best, S_best), select=1)
        T = T_full[0]
        S = S_full[0]
        H = H_full[0]

    ## 4. Obtain consistency at Gaussian models
    if correction == "empirical":
        # c1 correction (theoretical)
        S_corrected = S * (h/float(n_samples)) / chi2(n_features+2).cdf(chi2(n_features).ppf(h/float(n_samples)))
    elif correction == "theoretical":
        # c2 correction (empirical)
        X_centered = X - T
        dist = (np.dot(X_centered,linalg.inv(S)) * X_centered).sum(1)
        S_corrected = S * (np.median(dist)/chi2(n_features).isf(0.5))
    else:
        # no correction
        S_corrected = S
    
    ## 5. Reweight estimate
    if reweight == "rousseeuw":
        X_centered = X - T
        dist = (np.dot(X_centered, linalg.inv(S_corrected)) * X_centered).sum(1)
        mask = np.where(dist < chi2(n_features).isf(0.025))[0]
        T_reweighted = X[mask].mean(0)
        S_reweighted = empirical_covariance(X[mask])
        support_aux = np.zeros(n_samples).astype(bool)
        support_aux[np.where(mask)[0]] = True
        support = support_aux
    else:
        T_reweighted = T
        S_reweighted = S_corrected
        support = H
    
    return T_reweighted, S_reweighted, support
