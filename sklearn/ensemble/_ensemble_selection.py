"""Ensemble selection of classifiers and regressors"""

# Authors: Pierrick Pochelu <pierrick.pochelu@gmail.com>

from numbers import Integral
import math
import numpy as np
from joblib import Parallel, delayed

from ..utils._array_api import get_namespace
from ..base import is_classifier, is_regressor
from ._base import _fit_single_estimator
from ._base import _BaseHeterogeneousEnsemble
from ..utils.validation import column_or_1d
from ..exceptions import NotFittedError

from ..preprocessing import LabelEncoder
from ..utils.multiclass import type_of_target

from ..utils._param_validation import StrOptions
from ..utils.validation import check_is_fitted
from sklearn.metrics import log_loss, mean_squared_error
from sklearn import preprocessing

DEFAULT_REGRESSION_METRIC = mean_squared_error
DEFAULT_REGRESSION_DIRECTION = "min"
DEFAULT_CLASSFIER_METRIC = log_loss
DEFAULT_CLASSIFIER_DIRECTION = "min"


def check_and_convert_estimators(estimators):
    # convert estimators to a list of (name,estimator).
    # This bloc of code allows CI compliance.
    format = "named_est"
    if isinstance(estimators, list) and len(estimators) > 0:
        first_obj = estimators[0]
        if isinstance(first_obj, tuple) and len(first_obj) == 2:
            format = "named_est"
        else:
            format = "list_of_est"
    else:
        raise ValueError("Not expected container type")

    # conversion if required
    default_prefix = "est"
    if format == "list_of_est":
        named_estimators = []
        for i, e in enumerate(estimators):
            named_estimators.append((default_prefix + str(i), e))
        estimators = named_estimators
    return estimators


def get_task_type(estimators):
    if len(estimators) == 0:
        raise ValueError("No estimator")
    est = estimators[0][1]
    if is_classifier(est):
        return "classification"
    elif is_regressor(est):
        return "regression"
    else:
        raise ValueError("Unknown task type")


class EnsembleSelection(_BaseHeterogeneousEnsemble):
    """
    An ensemble classifier built by greedy stepwise selection.
    # Bagged Ensemble Selection (section 2.3) of the original paper is not implemented.

    Parameters
    ----------
    estimators: list of (string, estimator) tuples
        The estimators from which the ensemble selection classifier is built.
        These estimators must be fitted.

    score_metric: callable.
        Classification or regression

    min_estimators: integer, optional (default=1)
        The minimum number of base estimators.
        The ensemble selection estimator may overfit the calibration dataset.
        We may assume big ensemble generalize better. So min_estimators
        value allows to regularize.

    max_estimators: integer, optional (default=50)
        The maximum number of base estimators. It allows to control the final
        ensemble size, and the ensemble selection computing time.

    pruning_factor: float, optional (default=0.6)
        The worst pruning_factor estimators are remove to reduce the number of
        combinations

    with_replacement: bool, optional (defaut=True)
        If the same model can be replaced multiple times in the ensemble.
        When enabled, the prediction is the weighted averaging.
        The number of potential ensembles is given by the Combination formula:
        C(n,max_estimators) with n base estimators and min_estimators=0.
        When disabled, the  ensemble's prediction is the simple averaging.
        The number of potential ensembles is 2**n with n base estimators,
        min_estimators=0 and large max_estimators.


    is_base_estimator_proba: bool, optional (default=True)
        If True, estimators call "predict_proba(X)", otherwise, "predict(X)"

    n_jobs: int, default=None
        The number of jobs to run in parallel for ``fit``.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose: bool, default=False
        Computing time may take a long time
        If True, the model's selection will be printed.

    References
    ----------
    .. [1] R. Caruana, A. Niculescu-Mizil, G. Crew, A. Ksikes, "Ensemble
           Selection from Libraries of Models", 2004.
    .. [2] R. Caruana, A. Munson, A. Niculescu-Mizil, "Getting the Most Out
           of Ensemble Selection", 2006.

    """

    _parameter_constraints: dict = {
        "estimators": [list],
        "score_metric": [callable],
        "score_direction": [None, StrOptions({"min", "max"})],
        "max_estimators": [Integral],
        "min_estimators": [Integral],
        "pruning_factor": [Integral],
        "n_jobs": [None, Integral],
        "with_replacement": "boolean",
        "is_base_estimator_proba": "boolean",
        "verbose": "verbose",
    }

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
        super().__init__(estimators)
        self.estimators = check_and_convert_estimators(estimators)
        self.max_estimators = max_estimators
        self.min_estimators = min_estimators
        self.pruning_factor = pruning_factor
        self.verbose = verbose
        self.with_replacement = with_replacement
        self.score_metric = score_metric
        self.is_base_estimator_proba = is_base_estimator_proba
        self.score_direction = score_direction
        self.n_jobs = n_jobs

        # Detect the task
        task_type = get_task_type(self.estimators)

        # Default metric
        if self.score_metric is None or self.score_direction is None:
            if task_type == "regression":
                self.score_metric = DEFAULT_REGRESSION_METRIC
                self.score_direction = DEFAULT_REGRESSION_DIRECTION
            elif task_type == "classification":
                self.score_metric = DEFAULT_CLASSFIER_METRIC
                self.score_direction = DEFAULT_CLASSIFIER_DIRECTION
            else:
                raise ValueError("Unknown metric in this case")

    def _validate_estimators(self):
        """Overload the method of `_BaseHeterogeneousEnsemble` to be more
        lenient towards the type of `estimators`.

        Regressors can be accepted for some cases such as ordinal regression.
        """
        if len(self.estimators) == 0:
            raise ValueError(
                "Invalid 'estimators' attribute, 'estimators' should be a "
                "non-empty list of (string, estimator) tuples."
            )
        names, estimators = zip(*self.estimators)
        self._validate_names(names)

        has_estimator = any(est != "drop" for est in estimators)
        if not has_estimator:
            raise ValueError(
                "All estimators are dropped. At least one is required "
                "to be an estimator."
            )

        return names, estimators

    def fit(self, X, y, sample_weight=None):
        """Conduct ensemble selection on the validation set (X, y).
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            calibration data
        y : array-like, shape = [n_samples]
            Target values.

        sample_weight : array-like

        Returns
        -------
        self : object
        """
        # self._validate_params()
        # self._validate_estimators()

        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError(
                "Invalid `estimators` attribute, `estimators`"
                " should be a list of"
                " (string, fitted_estimator) tuples"
            )

        # TODO uses check_classification_targets(y) instead
        if len(np.array(y).shape) == 0:
            self.log("The labels are 0-dimensional")
            return {}

        # Fit estimator if not already done
        named_estimator_to_fit = []
        for name, estimator in self.estimators:
            try:
                check_is_fitted(estimator)
            except NotFittedError:
                named_estimator_to_fit.append((name, estimator))
        """
        fitted_estimators = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_single_estimator)(clone(est), X, y, sample_weight)
            for name, est in named_estimator_to_fit
        )
        """
        _ = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_single_estimator)(est, X, y, sample_weight)
            for name, est in named_estimator_to_fit
        )  # estimators are updated

        # converse sparse representation into one-hot vector
        if get_task_type(self.estimators) == "classification":
            y = self.y_preprocess_for_classification(y)
        else:
            y = self.y_preproces_for_regression(y)

        # Cache score and predict results on X
        self.base_model_score_ = {}
        self.base_model_preds_ = {}
        for name, est in self.estimators:
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
        nb_to_kept = math.ceil((1.0 - self.pruning_factor) * len(sorted_estimators))
        pruned_estimator = sorted_estimators[
            :nb_to_kept
        ]  # contains [(name_A,score_A),(name_B,score_B),...]
        pruned_est_name = [est[0] for est in pruned_estimator]
        self.log(f"pruned_estimator:{pruned_estimator}")

        # Build the ensemble
        (
            self.ensemble_,
            self.ensemble_score,
            self.ensemble_pred,
        ) = self._greedy_combi_ensemble(y, pruned_est_name)

        self.log(f"Ensemble :{self.ensemble_}, score: {self.ensemble_score}")

        return self

    def transform(self, X):
        preds_weights = []
        for estimator_name, estimator in self.estimators:
            estimator_weight = self.ensemble_.get(estimator_name, 0)
            if estimator_weight > 0:

                pred = self._estimator_predict(X, estimator)

                preds_weights.append((pred, estimator_weight))
        return preds_weights

    def decision_function(self, X):
        check_is_fitted(self)
        return self.transform(X)

    def log(self, txt):
        if self.verbose:
            print(txt)

    def _predict_proba_and_regression(self, X):
        preds_weights = self.transform(X)
        ensemble_preds = self._combination_rule(preds_weights)
        return ensemble_preds

    def predict_proba(self, X):
        return self._predict_proba_and_regression(X)

    def predict(self, X):
        task = get_task_type(self.estimators)
        if task == "classification":
            proba = self.predict_proba(X)
            pred = np.argmax(proba, axis=1)
        else:
            pred = self._predict_proba_and_regression(X)
        return pred

    def fit_predict(self, X, y, sample_weight=None):
        self.fit(X, y, sample_weight)
        return self.ensemble_pred

    def _get_score_of_ensemble(self, X, y, ensemble):
        preds_weights = []

        for estimator_name, weight in ensemble.items():
            proba = self.base_model_preds_[estimator_name]
            preds_weights.append((proba, weight))
        ens_preds = self._combination_rule(preds_weights)
        ens_score = self.score_metric(y, ens_preds)
        return ens_score, ens_preds

    def _estimator_predict(self, X, est):
        if self.is_base_estimator_proba and hasattr(est, "predict_proba"):
            return est.predict_proba(X)
        else:
            return est.predict(X)

    def _combination_rule(self, preds_weights):
        final_predictions = None
        sum_weights = 0
        for p, w in preds_weights:
            if final_predictions is None:
                final_predictions = p * w
            else:
                final_predictions += p * w
            sum_weights += w
        final_predictions /= sum_weights
        return final_predictions  # if final_predictions is not None else None

    def _sliding_combination_rule(
        self, current_ensemble_pred, prev_ensemble_size, added_estimator_name, y
    ):
        candidate_pred = self.base_model_preds_[added_estimator_name]

        if current_ensemble_pred is None:
            sliding_ensemble_pred = candidate_pred
        else:
            sliding_ensemble_pred = (
                current_ensemble_pred * prev_ensemble_size + candidate_pred
            ) / (prev_ensemble_size + 1.0)

        new_ens_score = self.score_metric(y, sliding_ensemble_pred)
        return {added_estimator_name: (new_ens_score, sliding_ensemble_pred)}

    def _get_score_of_model(self, X, y, estimator):
        pred = self._estimator_predict(X, estimator)
        score = self.score_metric(y, pred)
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
        current_ensemble_score = (
            -np.inf if self.score_direction == "max" else np.inf
        )  # worst
        current_ensemble_pred = None

        if len(estimators_name) == 0:
            return current_ensemble, current_ensemble_score, current_ensemble_pred

        with Parallel(n_jobs=self.n_jobs) as parallel:

            e = 0  # number of select model
            while (
                e < self.max_estimators
            ):  # while max iterations is not reached or no more improvement

                # Scan all potential model
                if self.n_jobs == 1:
                    list_info = [
                        self._sliding_combination_rule(
                            current_ensemble_pred, e, name, y
                        )
                        for name in estimators_name
                    ]
                else:
                    list_info = parallel(
                        delayed(self._sliding_combination_rule)(
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

    def get_classif_y_format(self, y):
        xp, _ = get_namespace(y)
        y = xp.asarray(y)
        shape = y.shape

        if len(shape) == 1:
            return "sparse"
        elif len(shape) == 2:
            if shape[1] == 1:
                return "sparse"
            elif shape[1] > 1:
                return "onehot"
            else:
                raise ValueError("y shape unexpected")
        else:
            raise ValueError("Y format not handled")

    def y_preproces_for_regression(self, y):
        return column_or_1d(y)

    """
    def y_preproces_for_classification(self, y):
        yfmt = self.get_classif_y_format(y)
        if yfmt == "sparse":
            classes=[]
            for _, est in self.estimators:
                classes.append(est.classes_)
            seen_classes=set.union( *list(set(c) for c in classes) )
            categ=list(range(max(seen_classes)+1))

            y = y.reshape((-1, 1))
            y = preprocessing.OneHotEncoder() \
                .fit_transform(y) \
                .toarray() \
                .squeeze()
        return y
    """

    def y_preprocess_for_classification(self, y):
        if type_of_target(y) == "multilabel-indicator":
            self._label_encoder = [LabelEncoder().fit(yk) for yk in y.T]
            self.classes_ = [le.classes_ for le in self._label_encoder]
            y = np.array(
                [
                    self._label_encoder[target_idx].transform(target)
                    for target_idx, target in enumerate(y.T)
                ]
            ).T
        else:
            # convert
            self._label_encoder = LabelEncoder().fit(y)
            self.classes_ = self._label_encoder.classes_
            y = self._label_encoder.transform(y)

            # TODO check if it is usefull
            y = y.reshape((-1, 1))
            y = preprocessing.OneHotEncoder().fit_transform(y).toarray().squeeze()
        return y
