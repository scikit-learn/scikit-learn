"""Ensemble selection of classifiers and regressors."""

# Authors: Pierrick Pochelu <pierrick.pochelu@gmail.com>

from abc import ABCMeta
from numbers import Integral, Real
import math
import numpy as np
import joblib

from ..base import clone
from ..base import is_classifier, is_regressor
from ..base import TransformerMixin

from ._base import _fit_single_estimator
from ._base import _BaseHeterogeneousEnsemble

from ..exceptions import NotFittedError
from ..preprocessing import LabelEncoder
from .. import preprocessing

from ..utils._array_api import get_namespace
from ..utils.validation import column_or_1d
from ..utils.multiclass import type_of_target
from ..utils._param_validation import StrOptions
from ..utils.validation import check_is_fitted

from ..metrics import log_loss, mean_squared_error

DEFAULT_REGRESSION_METRIC = mean_squared_error
DEFAULT_REGRESSION_DIRECTION = "min"
DEFAULT_CLASSFIER_METRIC = log_loss
DEFAULT_CLASSIFIER_DIRECTION = "min"


class EnsembleSelection(
    TransformerMixin, _BaseHeterogeneousEnsemble, metaclass=ABCMeta
):
    """An ensemble classifier built by greedy stepwise selection.

    Bagged Ensemble Selection (section 2.3) of the original paper is not implemented.

    Parameters
    ----------
    estimators : list of (str, Estimator) tuples
        The estimators from which the ensemble selection classifier is built.
        These estimators must be fitted.

    score_metric : callable, optional (default=None)
        Classification or regression 2 args function.
        If no value is given, log_loss or MSE is used automatically.

    score_direction : str in {"min","max"}, optional (default="min")
        If the score metric should be minimized or maximized.

    min_estimators : int, optional (default=1)
        The minimum number of base estimators.
        The ensemble selection estimator may overfit the calibration dataset.
        We may assume big ensemble generalize better. So min_estimators
        value allows to regularize.

    max_estimators : int, optional (default=50)
        The maximum number of base estimators. It allows to control the final
        ensemble size, and the ensemble selection computing time.

    pruning_factor : float, optional (default=0.6)
        The worst pruning_factor estimators are remove to reduce the number of
        combinations.

    with_replacement : bool, optional (defaut=True)
        If the same model can be replaced multiple times in the ensemble.
        When enabled, the prediction is the weighted averaging.
        The number of potential ensembles is given by the Combination formula:
        C(n,max_estimators) with n base estimators and min_estimators=0.
        When disabled, the  ensemble's prediction is the simple averaging.
        The number of potential ensembles is 2**n with n base estimators,
        min_estimators=0 and large max_estimators.

    is_base_estimator_proba : bool, optional (default=True)
        If True, estimators call `predict_proba(X)`, otherwise, `predict(X)`.
        In the doubt, let it to True.

    n_jobs : int, default=None
        The number of jobs to run in parallel for ``fit``.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose : bool, default=False
        Computing time may take a long time
        If True, the model's selection will be printed.

    Attributes
    ----------
    _estimator_type : string
        Type of the ensemble to build. Either "classifier" or "regressor".

    fitted_estimators_ : list of estimators
        The elements of the `estimators` parameter, having been fitted on the
        training data.
        `estimators` fitted outside EnsembleSelection are not fitted in the `fit`.
        `estimators` unfitted given to EnsembleSelection are fitted in the `fit`.

    score_metric_ : 2 args callable returning the score
        The first argument is the target and the second one the predictions

    score_direction_ : string in set("min", "max")
        Direction to optimize `score_metric_`.

    n_features_in_ : ndarray of shape (`n_features_in_`,)
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimators expose such an attribute when fit.

    base_model_score_ : dict[str, float]
        Associates the score to estimators name.

    base_model_preds_ : dict[str, ndarray(n_samples, n_outputs)]
        Associates the predictions to estimators name.

    ensemble_ : dict[str, int]
        associates the number of times an estimator is selected.
        If an estimator is not selected it is  either associated to
        0 value or is not in the dict.

    ensemble_score_ : float
        the score of the current ensemble `ensemble_` computed during :term:`fit`.

    ensemble_pred_ : ndarray(n_samples, n_outputs)
        cached predictions of the ensemble `ensemble_`  computed during :term:`fit`.

    See also
    ----------
    Stacking : Stack of estimators with a supervised learner as combination rule.

    References
    ----------
    .. [1] R. Caruana, A. Niculescu-Mizil, G. Crew, A. Ksikes, "Ensemble
           Selection from Libraries of Models", 2004.
    .. [2] R. Caruana, A. Munson, A. Niculescu-Mizil, "Getting the Most Out
           of Ensemble Selection", 2006.

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.neural_network import MLPClassifier
    >>> from sklearn.preprocessing import StandardScaler
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.ensemble import EnsembleSelection
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = load_iris(return_X_y=True)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 42)
    >>> rf=RandomForestClassifier(n_estimators=10, random_state=42)
    >>> mlp=make_pipeline(StandardScaler(), MLPClassifier(random_state=42))
    >>> estimators = [('rf', rf),('mlp', mlp)]
    >>> clf=EnsembleSelection(estimators)
    >>> score=clf.fit(X_train, y_train).score(X_test, y_test)
    >>> print(score)
    0.02...
    """

    _parameter_constraints: dict = {
        "estimators": [list],
        "score_metric": [None, callable],
        "score_direction": [None, StrOptions({"min", "max"})],
        "max_estimators": [Integral],
        "min_estimators": [Integral],
        "pruning_factor": [Real],
        "n_jobs": [None, Integral],
        "with_replacement": ["boolean"],
        "is_base_estimator_proba": ["boolean"],
        "verbose": ["verbose"],
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
        self.max_estimators = max_estimators
        self.min_estimators = min_estimators
        self.pruning_factor = pruning_factor
        self.verbose = verbose
        self.with_replacement = with_replacement
        self.score_metric = score_metric
        self.is_base_estimator_proba = is_base_estimator_proba
        self.score_direction = score_direction
        self.n_jobs = n_jobs

    @staticmethod
    def _check_and_convert_estimators(estimators):
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

    @staticmethod
    def _get_task_type(estimators):
        if len(estimators) == 0:
            raise ValueError("No estimator")
        est = estimators[0][1]
        if is_classifier(est):
            return "classifier"
        elif is_regressor(est):
            return "regressor"
        else:
            raise ValueError("Unknown task type")

    @staticmethod
    def _combination_rule(preds_weights):
        final_predictions = None
        sum_weights = 0
        for p, w in preds_weights:
            if final_predictions is None:
                final_predictions = p * w
            else:
                final_predictions += p * w
            sum_weights += w
        final_predictions /= sum_weights
        return final_predictions

    def fit(self, X, y, sample_weight=None):
        """Fit the estimator.

        If base estimators are not already fitted they are also fitted on X, y.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples, n_features)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.
            Note that this is supported only if all underlying estimators
            support sample weights.

        Returns
        -------
        self : object
            Returns a fitted instance of estimator.
        """
        super()._validate_params()

        if len(np.array(y).shape) == 0:
            self._log("The labels are 0-dimensional")
            return self

        self._init_estimators_and_metrics()

        # Convert y to array
        y = np.asarray(y)

        # sample_weight to array
        if sample_weight is not None:
            sample_weight = column_or_1d(sample_weight)
            sample_weight = np.asarray(sample_weight)
            if sample_weight.shape[0] != y.shape[0]:
                raise ValueError("sample_weight and y should be the same size")

        # Get X info
        if X.__class__.__name__ == "DataFrame":
            header = list(X.head())
            self.feature_names_in_ = np.array(header, dtype=object)
        elif X.__class__.__name__ == "ndarray":
            if len(X.shape) == 1:
                self.n_features_in_ = 1
            elif len(X.shape) == 2:
                self.n_features_in_ = X.shape[1]
            else:
                raise ValueError("Unexpected X shape")

        # Fit estimator if not already done
        named_estimator_to_fit = []
        for name, estimator in self.estimators:
            try:
                check_is_fitted(estimator)
            except NotFittedError:
                named_estimator_to_fit.append((name, estimator))

        self.fitted_estimators_ = self._fit_those_estimators_by_copy(
            named_estimator_to_fit, X, y, sample_weight
        )

        # converse sparse representation into one-hot vector
        estimator_type = EnsembleSelection._get_task_type(self.fitted_estimators_)
        if estimator_type == "classifier":
            y = self._y_preprocess_for_classification(y)
        elif estimator_type == "regressor":
            y = self._y_preproces_for_regression(y)
        else:
            raise ValueError(
                f"EnsembleSelection type not understood: '{estimator_type}'."
            )

        # Cache score and predict results on X
        self.base_model_score_ = {}
        self.base_model_preds_ = {}
        for name, est in self.fitted_estimators_:
            score, pred = self._get_score_of_model(X, y, est, sample_weight)
            self.base_model_score_[name] = score
            self.base_model_preds_[name] = pred
        self._log(f"Base estimators score:{self.base_model_score_}")

        # Prune them
        is_reverse = self.score_direction_ == "max"
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
        self._log(f"pruned_estimator:{pruned_estimator}")

        # Build the ensemble
        (
            self.ensemble_,
            self.ensemble_score_,
            self.ensemble_pred_,
        ) = self._greedy_combi_ensemble(y, pruned_est_name, sample_weight)

        self._log(f"Ensemble :{self.ensemble_}, score: {self.ensemble_score_}")

        return self

    def transform(self, X):
        """Return class labels or probabilities for X for each estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        y_preds : ndarray of shape (n_samples, n_estimators, n_outputs)
            Predictions for each estimator.
        """
        preds_weights = self._get_pred_weights(X)
        list_of_preds = np.array(
            [(pw[0]) for pw in preds_weights]
        )  # shape: n_estimators, n_sample, n_output
        preds = np.transpose(
            list_of_preds, axes=(1, 0, 2)
        )  # shape: n_sample, n_estimators, n_output
        return preds

    def decision_function(self, X):
        """Decision function for samples in `X` using the selected ensemble.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        decisions : ndarray of shape (n_samples, n_samples, n_outputs)
            Predictions for each estimator.
        """
        check_is_fitted(self)
        return self.transform(X)

    def predict_proba(self, X):
        """Predict class probabilities for `X` using the selected ensemble.

        If estimators are regressors, calling it is equivalent to `predict`.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        probabilities : ndarray of shape (n_samples, n_classes) or (n_output,)
            The class probabilities of the input samples.
        """
        return self._predict_proba_and_regression(X)

    def predict(self, X):
        """Predict target for X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            Predicted targets.
        """
        estimator_type = EnsembleSelection._get_task_type(self.estimators)
        if estimator_type == "classifier":
            proba = self.predict_proba(X)
            pred = np.argmax(proba, axis=1)
        elif estimator_type == "regressor":
            pred = self._predict_proba_and_regression(X)
        else:
            raise ValueError(f"Unexpected estimator type '{estimator_type}'")
        return pred

    def fit_predict(self, X, y, sample_weight=None):
        """Fit and predict target for X.

        It is faster than `fit(X,y).predict(X)` because predictions computed in `fit`
        are directly returned.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples, n_features)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.
            Note that this is supported only if all underlying estimators
            support sample weights.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            Predicted targets for X.
        """
        self.fit(X, y, sample_weight)
        return self.ensemble_pred_

    def score(self, X, y, sample_weight=None):
        """Score given input data `X` and target `y`.

        The metric used is given by `score_metric_`.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples, n_features)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.
            Note that this is supported only if all underlying estimators
            support sample weights.

        Returns
        -------
        score : float
            Score given by ``, samples importance are weighted with `sample_weight`.
        """
        check_is_fitted(self)
        pred = self._predict_proba_and_regression(X)
        score = self._weighted_score_metric(y, pred, sample_weight)
        return score

    def _log(self, txt):
        if self.verbose:
            print(txt)

    def _init_estimators_and_metrics(self):
        self.estimators = EnsembleSelection._check_and_convert_estimators(
            self.estimators
        )

        # Detect the task
        self._estimator_type = EnsembleSelection._get_task_type(self.estimators)

        # Get default metric
        if self.score_metric is None:
            if self._estimator_type == "regressor":
                self.score_metric_ = DEFAULT_REGRESSION_METRIC
                self.score_direction_ = DEFAULT_REGRESSION_DIRECTION
            elif self._estimator_type == "classifier":
                self.score_metric_ = DEFAULT_CLASSFIER_METRIC
                self.score_direction_ = DEFAULT_CLASSIFIER_DIRECTION
            else:
                raise ValueError("Unknown metric in this case")
        else:
            self.score_metric_ = self.score_metric
            self.score_direction_ = self.score_direction

    def _fit_one_estimator(self, est, X, y, sample_weight):
        try:
            _fit_single_estimator(est, X, y, sample_weight)
        except TypeError:
            # TODO clean check if compatible with sample_weight
            # For example Gaussian estimator cannot handle weights
            _fit_single_estimator(est, X, y)
        return est

    def _fit_those_estimators_by_copy(
        self, named_estimator_to_fit, X, y, sample_weight
    ):
        """`named_estimator_to_fit` is the list of (name, Estimator) with name a str.
        Fit the unfitted estimators.
        """

        # Parallel training of estimators
        new_fitted_est = joblib.Parallel(n_jobs=self.n_jobs)(
            joblib.delayed(self._fit_one_estimator)(clone(est), X, y, sample_weight)
            for name, est in named_estimator_to_fit
        )
        new_fitted_est_name = [named_est[0] for named_est in named_estimator_to_fit]

        new_named_est = dict()
        for name, est in zip(new_fitted_est_name, new_fitted_est):
            new_named_est[name] = est

        fitted_estimators = []
        for i, (n, old_estimator) in enumerate(self.estimators):
            if n in new_named_est:
                fitted_estimators.append((n, new_named_est[n]))
            else:
                fitted_estimators.append((n, old_estimator))
        return fitted_estimators

    def _predict_proba_and_regression(self, X):
        preds_weights = self._get_pred_weights(X)
        ensemble_preds = EnsembleSelection._combination_rule(preds_weights)
        return ensemble_preds

    def _get_score_of_ensemble(self, y, ensemble, sample_weight):
        """Not used in this file because an optimized version is implemented in"""
        preds_weights = []
        for estimator_name, weight in ensemble.items():
            proba = self.base_model_preds_[estimator_name]
            preds_weights.append((proba, weight))
        ens_preds = EnsembleSelection._combination_rule(preds_weights)
        ens_score = self._weighted_score_metric(y, ens_preds, sample_weight)
        return ens_score, ens_preds

    def _get_pred_weights(self, X):
        check_is_fitted(self)
        preds_weights = []
        for estimator_name, estimator in self.fitted_estimators_:
            estimator_weight = self.ensemble_.get(estimator_name, 0)
            if estimator_weight > 0:

                pred = self._estimator_predict(X, estimator)

                preds_weights.append((pred, estimator_weight))
        return preds_weights

    def _estimator_predict(self, X, est):
        if self.is_base_estimator_proba and hasattr(est, "predict_proba"):
            return est.predict_proba(X)
        else:
            return est.predict(X)

    def _progressive_metric(
        self,
        current_ensemble_pred,
        prev_ensemble_size,
        added_estimator_name,
        y,
        sample_weight,
    ):
        """
        Compute the score if we add one estimation to the current ensemble.

        returns : {added_estimator_name: (new_ens_score, sliding_ensemble_pred)}
        with `new_ens_score` the new score of the ensemble
        and `sliding_ensemble_pred` the new predictions
        """
        candidate_pred = self.base_model_preds_[added_estimator_name]

        if current_ensemble_pred is None:
            sliding_ensemble_pred = candidate_pred
        else:
            sliding_ensemble_pred = (
                current_ensemble_pred * prev_ensemble_size + candidate_pred
            ) / (prev_ensemble_size + 1.0)

        # uniform
        new_ens_score = self._weighted_score_metric(
            y, sliding_ensemble_pred, sample_weight
        )

        return {added_estimator_name: (new_ens_score, sliding_ensemble_pred)}

    def _get_score_of_model(self, X, y, estimator, sample_weight):
        pred = self._estimator_predict(X, estimator)
        score = self._weighted_score_metric(y, pred, sample_weight)
        return score, pred

    def _weighted_score_metric(self, y, pred, sample_weight):
        if sample_weight is None:
            score = self.score_metric_(y, pred)
        else:
            n = len(y)
            weighted_data_score = np.zeros((n,))
            for i in range(n):
                yi = np.array([y[i]])
                predi = np.array([pred[i]])
                wi = sample_weight[i]
                weighted_data_score[i] = self.score_metric_(yi, predi) * wi
            score = np.mean(weighted_data_score)
        return score

    def _greedy_combi_ensemble(self, y, estimators_name, sample_weight):
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
            -np.inf if self.score_direction_ == "max" else np.inf
        )  # worst
        current_ensemble_pred = None

        if len(estimators_name) == 0:
            return current_ensemble, current_ensemble_score, current_ensemble_pred

        with joblib.Parallel(n_jobs=self.n_jobs) as parallel:

            e = 0  # number of select model
            while (
                e < self.max_estimators
            ):  # while max iterations is not reached or no more improvement

                # Scan all potential model
                if self.n_jobs == 1:
                    list_info = [
                        self._progressive_metric(
                            current_ensemble_pred, e, name, y, sample_weight
                        )
                        for name in estimators_name
                    ]
                else:
                    list_info = parallel(
                        joblib.delayed(self._progressive_metric)(
                            current_ensemble_pred, e, estimator, y, sample_weight
                        )
                        for estimator in estimators_name
                    )
                # synchronization barrier when parallel
                local_choice_score = {}
                for info in list_info:
                    local_choice_score.update(info)
                for name, info in local_choice_score.items():
                    self._log(f"----> Candidate: {name}, score:{info[0]}")

                # Get the best local ensemble candidate
                sorted_local_choice_sorted = list(
                    sorted(
                        local_choice_score.items(),
                        key=lambda t: t[1][0],
                        reverse=self.score_direction_ == "max",
                    )
                )  # [ [name_A,(pred_A,score_A)], [name_B,(pred_B,score_B)], ... ]
                best_candidate_info = sorted_local_choice_sorted[0]
                best_candidate_ens_name = best_candidate_info[0]
                best_candidate_ens_pred = best_candidate_info[1][1]
                best_candidate_ens_score = best_candidate_info[1][0]

                # Adding it ?
                self._log(
                    "--> RESULT:"
                    f" iteration:{e+1}/{self.max_estimators}"
                    f" add:{best_candidate_ens_name}"
                    f" score:{best_candidate_ens_score}"
                )

                if (
                    self.score_direction_ == "max"
                    and best_candidate_ens_score > current_ensemble_score
                    or self.score_direction_ == "min"
                    and best_candidate_ens_score < current_ensemble_score
                ):
                    self._log("ADDING IT. The model improves the score of the ensemble")
                    adding_it = True
                elif e < self.min_estimators:
                    self._log(
                        "ADDING IT. The model decreases the score of the ensemble but"
                        " the minimum of models is not yet reached"
                    )
                    adding_it = True
                else:
                    self._log(
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
                    self._log(f"--> We add {best_candidate_ens_name}")

                    if not self.with_replacement:
                        estimators_name.remove(best_candidate_ens_name)

        return current_ensemble, current_ensemble_score, current_ensemble_pred

    def _get_classif_y_format(self, y):
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

    def _y_preproces_for_regression(self, y):
        return column_or_1d(y)

    def _y_preprocess_for_classification(self, y):
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

            y = y.reshape((-1, 1))
            y = preprocessing.OneHotEncoder().fit_transform(y).toarray().squeeze()
        return y
