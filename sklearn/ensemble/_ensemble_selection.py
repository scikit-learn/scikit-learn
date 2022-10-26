"""Ensemble selection of classifiers and regressors"""

# Authors: Pierrick Pochelu <pierrick.pochelu@gmail.com>

import math
import numpy as np
from joblib import Parallel, delayed
from ._base import _BaseHeterogeneousEnsemble
from ..base import ClassifierMixin
from ..base import RegressorMixin


class EnsembleSelection(ClassifierMixin, RegressorMixin, _BaseHeterogeneousEnsemble):
    """
    An ensemble classifier built by greedy stepwise selection.
    # Bagged Ensemble Selection (section 2.3) of the original paper is not implemented.

    Parameters
    ----------
    estimators : list of (string, estimator) tuples
        The estimators from which the ensemble selection classifier is built.
        These estimators must be fitted.

    score_metric : callable.
        Classification or regression

    min_estimators: integer, optional (default=1)
        The minimum number of base estimators.
        The calibration dataset may be overfitted. We may hypothesis big ensembles
        generalize better. So min_estimators value allows to regularize.

    max_estimators : integer, optional (default=50)
        The maximum number of base estimators. It allows to control the final
        ensemble size.

    pruning_factor : float, optional (default=0.6)
        The worst pruning_factor estimators are remove to reduce the number of
        combinations

    is_base_estimator_proba: bool, optional (default=True)
        If True, estimators call "predict_proba(X)", otherwise, "predict(X)"

    n_jobs : int, default=None
        The number of jobs to run in parallel for ``fit``.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose : bool, default=False
        Computing time may take a long time
        If True, the model's selection will be printed.

    References
    ----------
    .. [1] R. Caruana, A. Niculescu-Mizil, G. Crew, A. Ksikes, "Ensemble
           Selection from Libraries of Models", 2004.
    .. [2] R. Caruana, A. Munson, A. Niculescu-Mizil, "Getting the Most Out
           of Ensemble Selection", 2006.

    """

    _parameter_constraints: dict = {"estimators": [list], "score_metric": [callable]}

    def __init__(
        self,
        estimators=None,
        score_metric=None,
        score_direction="min",
        max_estimators=5,
        min_estimators=1,
        pruning_factor=0.4,
        n_jobs=None,
        with_replacement=True,
        is_base_estimator_proba=True,
        verbose=False,
    ):
        estimators_as_dict = {}
        for i, d in enumerate(estimators):
            estimators_as_dict[f"e{i}"] = d
        super().__init__(estimators_as_dict)
        self.max_estimators = max_estimators
        self.min_estimators = min_estimators
        self.pruning_factor = pruning_factor
        self.prune_fraction = pruning_factor
        self.verbose = verbose
        self.with_replacement = with_replacement
        self.score_metric_ = score_metric
        self.is_base_estimator_proba = is_base_estimator_proba
        self.score_direction = score_direction
        self.n_jobs = n_jobs

    def fit(self, X, y, sample_weight=None):
        """Conduct ensemble selection on the validation set (X, y).
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            calibration data
        y : array-like, shape = [n_samples]
            Target values.
        Returns
        -------
        self : object
        """
        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError(
                "Invalid `estimators` attribute, `estimators`"
                " should be a list of"
                " (string, fitted_estimator) tuples"
            )

        self.output_shape = list(y[0].shape)

        # Cache score and predict results on X
        self.base_model_score_ = {}
        self.base_model_preds_ = {}
        for name, est in self.named_estimators.items():
            score, pred = self._get_score_of_model(X, y, est)
            self.base_model_score_[name] = score
            self.base_model_preds_[name] = pred
        self.log(f"Base estimators score:{self.base_model_score_}")

        # Prune them
        is_reverse = self.score_direction == "max"
        sorted_estimators = list(
            sorted(
                self.base_model_score_.items(), key=lambda t: t[1], reverse=is_reverse
            )
        )
        nb_to_kept = math.ceil(1.0 - self.pruning_factor * len(sorted_estimators))
        pruned_estimator = sorted_estimators[
            :nb_to_kept
        ]  # contains [(name_A,score_A),(name_B,score_B),...]
        pruned_est_name = [est[0] for est in pruned_estimator]
        self.log(f"pruned_estimator:{pruned_estimator}")

        # Build the ensemble
        (
            self.ensemble,
            self.ensemble_score,
            self.ensemble_pred,
        ) = self._greedy_combi_ensemble(y, pruned_est_name)

        self.log(f"Ensemble :{self.ensemble}, score: {self.ensemble_score}")

        return self.ensemble

    def log(self, txt):
        if self.verbose:
            print(txt)

    def predict(self, X):
        ensemble_size = sum(self.ensemble.values())
        ens_pred = np.zeros([X.shape[0]] + self.output_shape)  # For instance

        for estimator_name, estimator_weight in self.ensemble.items():
            est = self.estimators[estimator_name]

            pred = self._estimator_predict(X, est)

            ens_pred += pred * estimator_weight

        ens_pred /= ensemble_size
        return ens_pred

    def fit_predict(self, X, y, sample_weight=None):
        self.fit(X, y, sample_weight)
        return self.ensemble_pred

    def _get_score_of_ensemble(self, X, y, ensemble):
        ensemble_size = sum(ensemble.values())
        ens_preds = np.zeros([X.shape[0]] + self.output_shape)

        for estimator_name, weight in ensemble.items():
            proba = self.base_model_preds_[estimator_name]
            ens_preds += proba * weight

        ens_preds /= ensemble_size
        ens_score = self.score_metric_(y, ens_preds)
        return ens_score, ens_preds

    def _estimator_predict(self, X, est):
        if self.is_base_estimator_proba:
            return est.predict_proba(X)
        else:
            return est.predict(X)
        # return est.predict(X)

    def _get_score_of_model(self, X, y, estimator):
        pred = self._estimator_predict(X, estimator)
        score = self.score_metric_(y, pred)
        return score, pred

    def _greedy_combi_ensemble(self, y, estimators_name):
        """
        Greedy Combinatorial Optimization to build ensembles
        Search the best subset estimators_name
        """

        # cur_ens_score, cur_ens_score = self._get_score_of_ensemble(X, y, ensemble)
        # n_clfs = sum(ensemble.values())

        # Init data structure
        current_ensemble = dict()
        for n in estimators_name:
            current_ensemble[n] = 0
        current_ensemble_score = -np.inf if self.score_direction == "max" else np.inf
        current_ensemble_pred = np.zeros(self.output_shape)

        with Parallel(n_jobs=self.n_jobs) as parallel:

            e = 0  # number of select model
            while e < self.max_estimators:

                # Scan all potential model
                if self.n_jobs == 1:
                    list_info = [
                        self.evaluate_candidate_ensemble(
                            current_ensemble_pred, e, estimator, y
                        )
                        for estimator in estimators_name
                    ]
                else:
                    list_info = parallel(
                        delayed(self.evaluate_candidate_ensemble)(
                            current_ensemble_pred, e, estimator, y
                        )
                        for estimator in estimators_name
                    )
                # synchronization barrier when parallel
                local_choice_score = {}
                for info in list_info:
                    local_choice_score.update(info)
                for name, info in local_choice_score.items():
                    self.log(f"----> Candidate: {name}, score:{info[0]}")

                # Get the best local ensemble candidate
                sorted_local_choice_sorted = list(
                    sorted(
                        local_choice_score.items(),
                        key=lambda t: t[1][0],
                        reverse=self.score_direction == "max",
                    )
                )  # [ [name_A,(pred_A,score_A)], [name_B,(pred_B,score_B)], ... ]
                best_candidate_info = sorted_local_choice_sorted[0]
                best_candidate_ens_name = best_candidate_info[0]
                best_candidate_ens_pred = best_candidate_info[1][1]
                best_candidate_ens_score = best_candidate_info[1][0]

                # Adding it ?
                self.log(
                    "--> RESULT:"
                    f" iteration:{e+1}/{self.max_estimators}"
                    f" add:{best_candidate_ens_name}"
                    f" score:{best_candidate_ens_score}"
                )

                if (
                    self.score_direction == "max"
                    and best_candidate_ens_score > current_ensemble_score
                    or self.score_direction == "min"
                    and best_candidate_ens_score < current_ensemble_score
                ):
                    self.log("ADDING IT. The model improves the score of the ensemble")
                    adding_it = True
                elif e < self.min_estimators:
                    self.log(
                        "ADDING IT. The model decreases the score of the ensemble but"
                        " the minimum of models is not yet reached"
                    )
                    adding_it = True
                else:
                    self.log(
                        "STOPPING. The result is not improved and minimum of models"
                        " reached"
                    )
                    e = np.inf  # We left the loop
                    adding_it = False

                # Adding it.
                if adding_it:
                    e += 1
                    current_ensemble[best_candidate_ens_name] += 1
                    current_ensemble_score = best_candidate_ens_score
                    current_ensemble_pred = best_candidate_ens_pred
                    self.log(f"--> We add {best_candidate_ens_name}")

                    if not self.with_replacement:
                        estimators_name.remove(best_candidate_ens_name)

        return current_ensemble, current_ensemble_score, current_ensemble_pred

    def evaluate_candidate_ensemble(
        self, current_ensemble_pred, prev_ensemble_size, added_estimator_name, y
    ):
        candidate_pred = self.base_model_preds_[added_estimator_name]
        candidate_ens_pred = (
            current_ensemble_pred * prev_ensemble_size + candidate_pred
        ) / (prev_ensemble_size + 1.0)
        new_ens_score = self.score_metric_(y, candidate_ens_pred)
        return {added_estimator_name: (new_ens_score, candidate_ens_pred)}
