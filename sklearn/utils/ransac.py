import numpy as np


def ransac(X, y, estimator_cls, min_samples, residual_threshold,
           is_data_valid=None, is_model_valid=None, max_trials=100,
           stop_n_inliers=np.inf, stop_score=1, estimator_kwargs={}):
    """Fit a model to data with the RANSAC (random sample consensus) algorithm.

    """

    best_estimator = None
    best_n_inliers = 0
    best_score = np.inf
    best_inlier_mask = None
    best_inlier_X = None
    best_inlier_y = None

    # estimator used for all iterations and for output
    estimator = estimator_cls(**estimator_kwargs)

    # number of data samples
    n_samples = X.shape[0]

    for _ in range(max_trials):

        # choose random sample set
        random_idxs = np.random.randint(0, n_samples, min_samples)
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
        if rsample_n_inliers < best_n_inliers:
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
        return None, None

    # estimate final model using all inliers
    estimator.fit(best_inlier_X, best_inlier_y)

    return estimator, best_inlier_mask
