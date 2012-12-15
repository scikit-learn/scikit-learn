import numpy as np

from .metrics import (r2_score, mean_squared_error, zero_one_score, f1_score,
                      auc_score, average_precision_score, precision_score)


class AsScorer(object):
    def __init__(self, score_func, greater_is_better=True,
                 needs_threshold=False, **kwargs):
        self.score_func = score_func
        self.greater_is_better = greater_is_better
        self.needs_threshold = needs_threshold
        self.kwargs = kwargs

    def __call__(self, estimator, X, y):
        if self.needs_threshold:
            if len(np.unique(y)) > 2:
                raise ValueError("This classification score only "
                                 "supports binary classification.")
            try:
                y_pred = estimator.decision_function(X).ravel()
            except (NotImplementedError, AttributeError):
                y_pred = estimator.predict_proba(X)[:, 1]
            return self.score_func(y, y_pred)
        else:
            y_pred = estimator.predict(X)
            return self.score_func(y, y_pred)


# Standard regression scores
R2Scorer = AsScorer(r2_score)
MSEScorer = AsScorer(mean_squared_error, greater_is_better=False)

# Standard Classification Scores
ZeroOneScorer = AsScorer(zero_one_score)
F1Scorer = AsScorer(f1_score)

# Score functions that need decision values
AUCScorer = AsScorer(auc_score, True, True)
AveragePrecisionScorer = AsScorer(average_precision_score,
                                  needs_threshold=True)
PrecisionScorer = AsScorer(precision_score)

scorers = dict(r2=R2Scorer, mse=MSEScorer, zero_one=ZeroOneScorer, f1=F1Scorer,
               auc=AUCScorer, ap=AveragePrecisionScorer,
               precision=PrecisionScorer)
