"""
TODO: Documentation

AISTAT 2023 `Sampling uncertainties on the Precision-Recall curve`
Baak, Max; Collot, Stéphane; Fridman Rojas, Ilan; Urlus, Ralph E.Q.


"""

import numpy as np
from scipy.stats import chi2, norm
from scipy.special import xlogy


def phat_PR(rec, prec, x_tp, x_fp, x_tn, x_fn):
    """Fit probability parameters of confusion matrix under the constraint of
    fixed recall and precision
    """
    n4 = x_tp + x_fp + x_tn + x_fn
    n3 = x_tp + x_fp + x_fn
    p_tp = n3 / (n4*(1/rec + 1/prec - 1))  # p_tp hat
    p_fn = ((1-rec)/rec) * p_tp  # remarq: rec >= epilson
    p_fp = ((1-prec)/prec) * p_tp  # remarq: prec >= epilson
    p_tn = 1. - p_fn - p_fp - p_tp
    # prevent negative values to due machine level noise
    if isinstance(p_tn, np.ndarray):
        p_tn[p_tn < 0] = 0
    elif isinstance(p_tn, float) and p_tn < 0:
        p_tn = 0.
    return p_tp, p_fp, p_tn, p_fn


def phat_ROC(fpr, tpr, x_tp, x_fp, x_tn, x_fn):
    """Fit probability parameters of confusion matrix under the constraint of
    fixed FPR and TPR
    """
    n4 = x_tp + x_fp + x_tn + x_fn
    p_tp = (tpr*(x_fn+x_tp)) / n4
    p_fn = (1-tpr)/tpr * p_tp
    p_fp = (fpr*(tpr-p_tp)) / tpr
    p_tn = 1. - p_fn - p_fp - p_tp
    # prevent negative values to due machine level noise
    if isinstance(p_tn, np.ndarray):
        p_tn[p_tn < 0] = 0
    elif isinstance(p_tn, float) and p_tn < 0:
        p_tn = 0.
    return p_tp, p_fp, p_tn, p_fn


def nll(rec, prec, x_tp, x_fp, x_tn, x_fn, phat_fnc):
    """Return -2logp of multinomial distribution fixed at certain point on the curve
    either for precision-recall, or for ROC by choosing the corresponding phat_fnc()

    Two steps:
    1. Fit with fixed recall and precision
    2. Fit with all probability parameters free

    Return the difference in -2 log L
    """
    # optimal fit of x
    n4 = x_tp + x_fp + x_tn + x_fn
    p_fn0 = x_fn / n4
    p_tp0 = x_tp / n4
    p_fp0 = x_fp / n4
    p_tn0 = x_tn / n4
    nll_minimum = -2 * xlogy(x_tp, p_tp0) - 2 * xlogy(x_fp, p_fp0) - 2 * xlogy(x_fn, p_fn0) - 2 * xlogy(x_tn, p_tn0)

    # fit of x constrained to recall and precision
    p_tp, p_fp, p_tn, p_fn = phat_fnc(rec, prec, x_tp, x_fp, x_tn, x_fn)
    nll_value = -2 * xlogy(x_tp, p_tp) - 2 * xlogy(x_fp, p_fp) - 2 * xlogy(x_fn, p_fn) - 2 * xlogy(x_tn, p_tn)

    # return the difference
    return nll_value - nll_minimum


def get_range_axis(term1, term2, epsilon, n_sigma):
    """
    Works for all of those: Recall, Precision, FPR and TPR = term1/(term1+term2)

    FPR       = x_fp / (x_fp + x_tn) # x-axis
    TPR       = x_tp / (x_tp + x_fn) # y-axis == recall
    recall    = x_tp / (x_tp + x_fn) # x-axis
    precision = x_tp / (x_tp + x_fp) # y-avis

    # Sigma estimation based on the covariance matrix first-order approximation:
    sigma FPR       = (x_fp*x_tn) / (x_fp + x_tn)**3
    sigma TPR       = (x_tp*x_fn) / (x_tp + x_fn)**3
    sigma recall    = (x_tp*x_fn) / (x_tp + x_fn)**3 == TPR
    sigma precision = (x_tp*x_fp) / (x_tp + x_fp)**3
    In all these case we can notice that:
    sigma XXXX      = (term1 * term2) / (term1 + term2)**3

    Remark: if you observe at least one fp, precision cannot be 100%
    """
    V = term1/(term1+term2)

    # Get sigma estimation based on the covariance matrix first-order approximation
    # If we get one the term of the product to be zero (like x_fp * x_tn)
    # we set it in order to have at least one x_fp, or x_tn, so that we have a non-zero sigma
    if term1 == 0:
        term1 = 1
    if term2 == 0:
        term2 = 1
    sigma_V = np.sqrt((term1 * term2) / (term1+term2)**3)

    # We introduce an epsilon to prevent division by zero at the edge because in phat() we have some divisions by precision, recall and TPR
    # but most importantly, we need an epsilon to have nice contour plots: (for X and Y)
    # if X == 1, its means that we can have a grid including value 1, and the contour will have finite values
    # if X <  1, its means that the probability to have value 1 is 0, and therefore the nll is infinity,
    #            and therefore we have a plotting issue for the contour, so we set it to 1-epsilon so the contour knows how to extrapolate
    max_V_clip = 1 if V == 1 else 1-epsilon

    # Ranges of values for the axis to scan, with clipping to not draw outside of the square (0,1)
    V_max = min(V + n_sigma * sigma_V, max_V_clip)  # max_V_clip to have nice contours
    V_min = max(V - n_sigma * sigma_V, epsilon)  # epsilon to prevent division by 0

    return V_max, V_min


def get_scaling_factor(norm_n_std):
    # Get the scale for 2 degrees of freedom confidence interval
    # We use chi2 because the equation of an ellipse is a sum of squared variable,
    # more details here https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    # norm_n_std = 1  # number of standard deviation
    norm_pct = 2. * (norm.cdf(norm_n_std) - 0.5)
    chi2_quantile = chi2.ppf(norm_pct, 2)
    scale = np.sqrt(chi2_quantile)
    return scale

def get_confusion_matrix(y_true, y_pred, thresholds):
    N = len(y_true)

    # remark: computing them with metrics.confusion_matrix() takes too much time
    P = np.array([sum(y_true)] * len(thresholds))
    # we use ">= thr" like in precision_recall_curve():
    TP = np.array([((y_pred >= thr) & y_true).sum() for thr in thresholds])
    PP = np.array([(y_pred >= thr).sum() for thr in thresholds])
    FN = P - TP
    FP = PP - TP
    TN = N - TP - FP - FN

    return TP, FP, TN, FN


def compute_sampling_uncertainty(curve, y_true, y_pred, thresholds, norm_n_std, n_bins, epsilon=1e-4):
    """Compute sampling uncertainty, with profile likelihoods based on Wilks’ theorem.
    It consists of the following steps:

    1. Get the curve
    2. Get the confusion matrix of each point of the curve
    3. For each observed point of the curve, estimate a surrounding 6 (i.e. more than the desired number) sigmas uncertainty grid rectangle (based on first-order approximation of the covariance matrix, with the bivariate normal distribution assumption)
    4. For each of these hypothesis point in the grid, compute the test static with the observed point, called the profile log likelihood ratio (using the fact that the confusion matrix follows a multinomial distribution).

    More details in:
    AISTAT 2023 `Sampling uncertainties on the Precision-Recall curve`
    Baak, Max; Collot, Stéphane; Fridman Rojas, Ilan; Urlus, Ralph E.Q.


    Parameters
    ----------
    curve : str
        Name of the curve, supported values: `precision_recall` and `ROC`.

    y_true : array-like of shape (n_samples,)
        True binary labels.

    y_pred : array-like of shape (n_samples,)
        Estimated probabilities or output of decision function.

    thresholds : ndarray of shape (n_thresholds,)
        Thresholds corresponding to the points to plot.

    norm_n_std : int
        Number of standard deviation to plot for sampling uncertainty level.

    n_bins : int
        Number of bins to use for the 2D grid to compute uncertainty for each point.

    epsilon : float, default=1e-4
        Small number to avoid division by zero and plotting edge cases
    """

    TP, FP, TN, FN = get_confusion_matrix(y_true, y_pred, thresholds)

    # Create the grid, used for all the curve points
    ls = np.linspace(0+epsilon, 1-epsilon, n_bins)  #[1:-1]
    RX, PY = np.meshgrid(ls, ls)
    CHI2 = np.full_like(RX, np.inf)

    # For each point in the curve...
    for x_tp, x_fp, x_tn, x_fn in zip(TP, FP, TN, FN):

        if curve == "precision_recall":
            phat = phat_PR
            # x-axis: rec = x_tp / (x_tp + x_fn)
            X_term1 = x_tp
            X_term2 = x_fn

            # y-axis: prec = x_tp / (x_tp + x_fp)  # same as TPR
            Y_term1 = x_tp
            Y_term2 = x_fp
        elif curve == "ROC":
            phat = phat_ROC
            # x-axis: FPR = x_fp / (x_fp + x_tn)
            X_term1 = x_fp
            X_term2 = x_tn

            # y-axis: TPR = x_tp / (x_tp + x_fn)  # same as rec
            Y_term1 = x_tp
            Y_term2 = x_fn
        else:
            raise ValueError(f"Unknown curve: {curve}")

        # Get a rough range of values to evaluate via approximation to cover a larger norm_n_std + 3 sigmas
        x_max, x_min = get_range_axis(X_term1, X_term2, epsilon, norm_n_std+3)
        y_max, y_min = get_range_axis(Y_term1, Y_term2, epsilon, norm_n_std+3)

        # Get the matrix position index
        x_min_idx = np.searchsorted(ls, x_min, side='left')
        x_max_idx = np.searchsorted(ls, x_max, side='right')
        y_min_idx = np.searchsorted(ls, y_min, side='left')
        y_max_idx = np.searchsorted(ls, y_max, side='right')

        # Get a submatrix view
        RX_v = RX[y_min_idx:y_max_idx, x_min_idx:x_max_idx]
        PY_v = PY[y_min_idx:y_max_idx, x_min_idx:x_max_idx]

        # Compute chi2 for this view
        CHI2_v = nll(RX_v, PY_v, x_tp, x_fp, x_tn, x_fn, phat)

        # Aggregate with the entire grid, and keep always the minimum chi2
        CHI2[y_min_idx:y_max_idx, x_min_idx:x_max_idx] = np.minimum(CHI2[y_min_idx:y_max_idx, x_min_idx:x_max_idx], CHI2_v)

    return RX, PY, CHI2, norm_n_std


def plot_sampling_uncertainty(ax, sampling_uncertainty):
    """Plot the contour (i.e. isoline) for the observed points.
    Using Wilks’ theorem stating that the profile log likelihood ratio is described asymptotically by a chi2 distribution.

    Parameters
    ----------
    ax : Matplotlib Axes
        Axes object to plot on.

    sampling_uncertainty : tuple (RX, RY, CHI2)
        The sampling uncertainty for each point on the curve.

    norm_n_std : int
        Number of standard deviation to plot for sampling uncertainty level.
        Relevant only if plot_uncertainty = True.
    """
    if sampling_uncertainty is None:
        raise ValueError("Sampling uncertainty is None, so we cannot plot it.")
    
    RX, PY, CHI2, norm_n_std = sampling_uncertainty

    scale = get_scaling_factor(norm_n_std)
    levels = [scale**2]
    levels = [0.0] + levels#.tolist()

    # For each point in the curve plot a contour based on uncertainty level
    ax.contourf(RX, PY, CHI2, levels=levels, alpha=0.50, colors='lightblue')
