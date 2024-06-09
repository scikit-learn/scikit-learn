from numbers import Integral

import numpy as np

from ._scorer import _BaseScorer


def _threshold_scores_to_class_labels(y_score, threshold, classes, pos_label):
    """Threshold `y_score` and return the associated class labels."""
    if pos_label is None:
        map_thresholded_score_to_label = np.array([0, 1])
    else:
        pos_label_idx = np.flatnonzero(classes == pos_label)[0]
        neg_label_idx = np.flatnonzero(classes != pos_label)[0]
        map_thresholded_score_to_label = np.array([neg_label_idx, pos_label_idx])

    return classes[map_thresholded_score_to_label[(y_score >= threshold).astype(int)]]


class _CurveScorer(_BaseScorer):
    """Scorer taking a continuous response and output a score for each threshold.

    Parameters
    ----------
    score_func : callable
        The score function to use. It will be called as
        `score_func(y_true, y_pred, **kwargs)`.

    sign : int
        Either 1 or -1 to returns the score with `sign * score_func(estimator, X, y)`.
        Thus, `sign` defined if higher scores are better or worse.

    kwargs : dict
        Additional parameters to pass to the score function.

    thresholds : int or array-like
        Related to the number of decision thresholds for which we want to compute the
        score. If an integer, it will be used to generate `thresholds` thresholds
        uniformly distributed between the minimum and maximum predicted scores. If an
        array-like, it will be used as the thresholds.

    response_method : str
        The method to call on the estimator to get the response values.
    """

    def __init__(self, score_func, sign, kwargs, thresholds, response_method):
        super().__init__(
            score_func=score_func,
            sign=sign,
            kwargs=kwargs,
            response_method=response_method,
        )
        self._thresholds = thresholds

    @classmethod
    def from_scorer(cls, scorer, response_method, thresholds):
        """Create a continuous scorer from a normal scorer."""
        instance = cls(
            score_func=scorer._score_func,
            sign=scorer._sign,
            response_method=response_method,
            thresholds=thresholds,
            kwargs=scorer._kwargs,
        )
        # transfer the metadata request
        instance._metadata_request = scorer._get_metadata_request()
        return instance

    def _score(self, method_caller, estimator, X, y_true, **kwargs):
        """Evaluate predicted target values for X relative to y_true.

        Parameters
        ----------
        method_caller : callable
            Returns predictions given an estimator, method name, and other
            arguments, potentially caching results.

        estimator : object
            Trained estimator to use for scoring.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Test data that will be fed to estimator.predict.

        y_true : array-like of shape (n_samples,)
            Gold standard target values for X.

        **kwargs : dict
            Other parameters passed to the scorer. Refer to
            :func:`set_score_request` for more details.

        Returns
        -------
        scores : ndarray of shape (thresholds,)
            The scores associated to each threshold.

        potential_thresholds : ndarray of shape (thresholds,)
            The potential thresholds used to compute the scores.
        """
        pos_label = self._get_pos_label()
        y_score = method_caller(
            estimator, self._response_method, X, pos_label=pos_label
        )

        scoring_kwargs = {**self._kwargs, **kwargs}
        if isinstance(self._thresholds, Integral):
            potential_thresholds = np.linspace(
                np.min(y_score), np.max(y_score), self._thresholds
            )
        else:
            potential_thresholds = np.asarray(self._thresholds)
        score_thresholds = [
            self._sign
            * self._score_func(
                y_true,
                _threshold_scores_to_class_labels(
                    y_score, th, estimator.classes_, pos_label
                ),
                **scoring_kwargs,
            )
            for th in potential_thresholds
        ]
        return np.array(score_thresholds), potential_thresholds
