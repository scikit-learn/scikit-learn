import numpy as np

from . import (r2_score, mean_squared_error, accuracy_score, f1_score,
               auc_score, average_precision_score, precision_score,
               recall_score)

from .cluster import adjusted_rand_score


class AsScorer(object):
    """Flexible scores for any estimator.

    This class wraps estimator scoring functions for the use in GridSearchCV
    and cross_val_score. It takes a score function, such as ``accuracy_score``,
    ``mean_squared_error``, ``adjusted_rand_index`` or ``average_precision``
    and provides a call method.

    Parameters
    ----------
    score_func : callable,
        Score function (or loss function) with signature
        ``score_func(y, y_pred, **kwargs)``.

    greater_is_better : boolean, default=True
        Whether score_func is a score function (default), meaning high is good,
        or a loss function, meaning low is good.

    needs_threshold : bool, default=False
        Whether score_func takes a continuous decision certainty.
        For example ``average_precision`` or the area under the roc curve
        can not be computed using predictions alone, but need the output of
        ``decision_function`` or ``predict_proba``.

    Examples
    --------
    >>> from sklearn.metrics import fbeta_score, AsScorer
    >>> ftwo_scorer = AsScorer(fbeta_score, beta=2)
    >>> from sklearn.grid_search import GridSearchCV
    >>> from sklearn.svm import LinearSVC
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]},
    ...                     scoring=ftwo_scorer)
    """
    def __init__(self, score_func, greater_is_better=True,
                 needs_threshold=False, **kwargs):
        self.score_func = score_func
        self.greater_is_better = greater_is_better
        self.needs_threshold = needs_threshold
        self.kwargs = kwargs

    def __call__(self, estimator, X, y):
        """Score X and y using the provided estimator.

        Parameters
        ----------
        estimator : object
            Trained estimator to use for scoring.
            If ``needs_threshold`` is True, estimator needs
            to provide ``decision_function`` or ``predict_proba``.
            Otherwise, estimator needs to provide ``predict``.

        X : array-like or sparse matrix
            Test data that will be scored by the estimator.

        y : array-like
            True prediction for X.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        if self.needs_threshold:
            if len(np.unique(y)) > 2:
                raise ValueError("This classification score only "
                                 "supports binary classification.")
            try:
                y_pred = estimator.decision_function(X).ravel()
            except (NotImplementedError, AttributeError):
                y_pred = estimator.predict_proba(X)[:, 1]
            return self.score_func(y, y_pred, **self.kwargs)
        else:
            y_pred = estimator.predict(X)
            return self.score_func(y, y_pred, **self.kwargs)


# Standard regression scores
r2_scorer = AsScorer(r2_score)
mse_scorer = AsScorer(mean_squared_error, greater_is_better=False)

# Standard Classification Scores
accuracy_scorer = AsScorer(accuracy_score)
f1_scorer = AsScorer(f1_score)

# Score functions that need decision values
auc_scorer = AsScorer(auc_score, True, True)
average_precision_scorer = AsScorer(average_precision_score,
                                    needs_threshold=True)
precision_scorer = AsScorer(precision_score)
recall_scorer = AsScorer(recall_score)

# Clustering scores
ari_scorer = AsScorer(adjusted_rand_score)

scorers = dict(r2=r2_scorer, mse=mse_scorer, accuracy=accuracy_scorer,
               f1=f1_scorer, roc_auc=auc_scorer,
               average_precision=average_precision_scorer,
               precision=precision_scorer, recall=recall_scorer,
               ari=ari_scorer)


