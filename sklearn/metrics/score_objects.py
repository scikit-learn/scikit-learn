import numpy as np

from .metrics import (r2_score, mean_squared_error, zero_one_score, f1_score,
                      auc_score, average_precision_score)


def _score_obj_from_func(score_func, greater_is_better_=True,
                         needs_threshold=False):
    # construct a scoring object from a simple scoring function
    if needs_threshold:
        class ScoreObj(object):
            greater_is_better = greater_is_better_

            def __call__(self, estimator, X, y):
                if len(np.unique(y)) > 2:
                    raise ValueError("This classification score only "
                                     "supports binary classification.")
                try:
                    y_pred = estimator.decision_function(X).ravel()
                except (NotImplementedError, AttributeError):
                    y_pred = estimator.predict_proba(X)[:, 1]
                return score_func(y, y_pred)
    else:
        class ScoreObj(object):
            greater_is_better = greater_is_better_

            def __call__(self, estimator, X, y):
                y_pred = estimator.predict(X)
                return score_func(y, y_pred)
    return ScoreObj


def _score_obj_from_string(score_func):
    # construct scoring object from string
    # TODO
    pass

# Standard regression scores
R2Score = _score_obj_from_func(r2_score, True)
MSEScore = _score_obj_from_func(mean_squared_error, False)

# Standard Classification Scores
ZeroOneScore = _score_obj_from_func(zero_one_score, False)
F1Score = _score_obj_from_func(f1_score, False)

# Score functions that need decision values
AUCScore = _score_obj_from_func(auc_score, True, True)
AveragePrecisionScore = _score_obj_from_func(average_precision_score, True,
                                             True)
