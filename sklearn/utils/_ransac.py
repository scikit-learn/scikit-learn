import numpy as np


def ransac(X, y, estimator, min_n_samples, residual_threshold,
           is_data_valid=None, is_model_valid=None, max_trials=100,
           stop_n_inliers=np.inf, stop_score=np.inf):
    """Fit a model to data with the RANSAC (random sample consensus) algorithm.

    Parameters
    ----------
    X : numpy array or sparse matrix of shape [n_samples, n_features]
        Training data.

    y : numpy array of shape [n_samples, n_targets]
        Target values.

    estimator : object
        Estimator object which implements the following methods:
        * `fit(X, y)`: Fit model to given  training data and target values.
        * `predict(X)`: Predict using the estimated model.
        * `score(X)`: Returns the mean accuracy on the given test data.

    residual_threshold : float
        Maximum residual for a data sample to be classified as an inlier.

    is_data_valid : function, optional
        This function is called with the randomly selected data before the
        model is fitted to it: `is_data_valid(X, y)`.

    is_model_valid : function, optional
        This function is called with the estimated model and the randomly
        selected data: `is_model_valid(model, X, y)`, .

    max_trials : int, optional
        Maximum number of iterations for random sample selection.

    stop_n_inliers : int, optional
        Stop iteration if at least this number of inliers are found.

    stop_score : float, optional
        Stop iteration if score is greater equal than this threshold.

    Returns
    -------
    n_trials : int
        Number of random selection trials.

    inlier_mask : bool array of shape [n_samples]
        Boolean mask of inliers classified as ``True``.

    Raises
    ------
    ValueError: If no valid consensus set could be found.

    Notes
    -----
    RANSAC is an iterative algorithm for the robust estimation of parameters
    from a subset of inliers from the complete data set. Each iteration
    performs the following steps:

    1. Select `min_n_samples` random samples from the original data and check
       whether the set of data is valid (see `is_data_valid`).
    2. Fit a model to the random subset (`estimator.fit`) and check whether
       the estimated model is valid (see `is_model_valid`).
    3. Classify all data as inliers or outliers by calculating the residuals
       to the estimated model (`estimator.predict(X) - y`) - all data samples
       with absolute residuals smaller than the `residual_threshold` are
       considered as inliers.
    4. Save fitted model as best model if number of inlier samples is
       maximal. In case the current estimated model has the same number of
       inliers, it is only considered as the best model if it has better score.

    These steps are performed either a maximum number of times (`max_trials`)
    or until one of the special stop criteria are met (see `stop_n_inliers` and
    `stop_score`). The final model is estimated using all inlier samples of the
    previously determined best model.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/RANSAC
    .. [2] http://www.cs.columbia.edu/~belhumeur/courses/compPhoto/ransac.pdf
    .. [3] http://www.bmva.org/bmvc/2009/Papers/Paper355/Paper355.pdf

    """

    best_n_inliers = 0
    best_score = np.inf
    best_inlier_mask = None
    best_inlier_X = None
    best_inlier_y = None

    # number of data samples
    n_samples = X.shape[0]

    for n_trials in range(max_trials):

        # choose random sample set
        random_idxs = np.random.randint(0, n_samples, min_n_samples)
        rsample_X = X[random_idxs]
        rsample_y = y[random_idxs]

        # check if random sample set is valid
        if is_data_valid is not None and not is_data_valid(X, y):
            continue

        # fit model for current random sample set
        estimator.fit(rsample_X, rsample_y)

        # check if estimated model is valid
        if is_model_valid is not None and not is_model_valid(estimator,
                                                             rsample_X,
                                                             rsample_y):
            continue

        # residuals of all data for current random sample model
        rsample_residuals = np.abs(estimator.predict(X) - y)

        # classify data into inliers and outliers
        rsample_inlier_mask = rsample_residuals < residual_threshold
        rsample_n_inliers = np.sum(rsample_inlier_mask)

        # less inliers -> skip current random sample
        if rsample_n_inliers < best_n_inliers \
           or rsample_n_inliers == best_n_inliers == 0:
            continue

        # extract inlier data set
        rsample_inlier_X = X[rsample_inlier_mask]
        rsample_inlier_y = y[rsample_inlier_mask]

        # score of inlier data set
        rsample_score = estimator.score(rsample_inlier_X, rsample_inlier_y)

        # same number of inliers but worse score -> skip current random sample
        if rsample_n_inliers == best_n_inliers and rsample_score < best_score:
            continue

        # save current random sample as best sample
        best_n_inliers = rsample_n_inliers
        best_score = rsample_score
        best_inlier_mask = rsample_inlier_mask
        best_inlier_X = rsample_inlier_X
        best_inlier_y = rsample_inlier_y

        # break if sufficient number of inliers or score is reached
        if best_n_inliers >= stop_n_inliers or best_score >= stop_score:
            break

    # if none of the iterations met the required criteria
    if best_inlier_mask is None:
        raise ValueError("RANSAC could not find valid consensus set.")

    # estimate final model using all inliers
    estimator.fit(best_inlier_X, best_inlier_y)

    return n_trials + 1, best_inlier_mask
