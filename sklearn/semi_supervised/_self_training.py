import warnings
from numbers import Integral, Real

import numpy as np

from ..base import BaseEstimator, MetaEstimatorMixin, _fit_context, clone
from ..utils import Bunch, safe_mask
from ..utils._param_validation import HasMethods, Interval, StrOptions
from ..utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _raise_for_params,
    _routing_enabled,
    process_routing,
)
from ..utils.metaestimators import available_if
from ..utils.validation import check_is_fitted

__all__ = ["SelfTrainingClassifier"]

# Authors: Oliver Rausch   <rauscho@ethz.ch>
#          Patrice Becker  <beckerp@ethz.ch>
# License: BSD 3 clause


def _estimator_has(attr):
    """Check if we can delegate a method to the underlying estimator.

    First, we check the fitted `base_estimator_` if available, otherwise we check
    the unfitted `base_estimator`. We raise the original `AttributeError` if
    `attr` does not exist. This function is used together with `available_if`.
    """

    def check(self):
        if hasattr(self, "base_estimator_"):
            getattr(self.base_estimator_, attr)
        else:
            getattr(self.base_estimator, attr)

        return True

    return check


class SelfTrainingClassifier(MetaEstimatorMixin, BaseEstimator):
    """Self-training classifier.

    This :term:`metaestimator` allows a given supervised classifier to function as a
    semi-supervised classifier, allowing it to learn from unlabeled data. It
    does this by iteratively predicting pseudo-labels for the unlabeled data
    and adding them to the training set.

    The classifier will continue iterating until either max_iter is reached, or
    no pseudo-labels were added to the training set in the previous iteration.

    Read more in the :ref:`User Guide <self_training>`.

    Parameters
    ----------
    base_estimator : estimator object
        An estimator object implementing `fit` and `predict_proba`.
        Invoking the `fit` method will fit a clone of the passed estimator,
        which will be stored in the `base_estimator_` attribute.

    threshold : float, default=0.75
        The decision threshold for use with `criterion='threshold'`.
        Should be in [0, 1). When using the `'threshold'` criterion, a
        :ref:`well calibrated classifier <calibration>` should be used.

    criterion : {'threshold', 'k_best'}, default='threshold'
        The selection criterion used to select which labels to add to the
        training set. If `'threshold'`, pseudo-labels with prediction
        probabilities above `threshold` are added to the dataset. If `'k_best'`,
        the `k_best` pseudo-labels with highest prediction probabilities are
        added to the dataset. When using the 'threshold' criterion, a
        :ref:`well calibrated classifier <calibration>` should be used.

    k_best : int, default=10
        The amount of samples to add in each iteration. Only used when
        `criterion='k_best'`.

    max_iter : int or None, default=10
        Maximum number of iterations allowed. Should be greater than or equal
        to 0. If it is `None`, the classifier will continue to predict labels
        until no new pseudo-labels are added, or all unlabeled samples have
        been labeled.

    verbose : bool, default=False
        Enable verbose output.

    Attributes
    ----------
    base_estimator_ : estimator object
        The fitted estimator.

    classes_ : ndarray or list of ndarray of shape (n_classes,)
        Class labels for each output. (Taken from the trained
        `base_estimator_`).

    transduction_ : ndarray of shape (n_samples,)
        The labels used for the final fit of the classifier, including
        pseudo-labels added during fit.

    labeled_iter_ : ndarray of shape (n_samples,)
        The iteration in which each sample was labeled. When a sample has
        iteration 0, the sample was already labeled in the original dataset.
        When a sample has iteration -1, the sample was not labeled in any
        iteration.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_iter_ : int
        The number of rounds of self-training, that is the number of times the
        base estimator is fitted on relabeled variants of the training set.

    termination_condition_ : {'max_iter', 'no_change', 'all_labeled'}
        The reason that fitting was stopped.

        - `'max_iter'`: `n_iter_` reached `max_iter`.
        - `'no_change'`: no new labels were predicted.
        - `'all_labeled'`: all unlabeled samples were labeled before `max_iter`
          was reached.

    See Also
    --------
    LabelPropagation : Label propagation classifier.
    LabelSpreading : Label spreading model for semi-supervised learning.

    References
    ----------
    :doi:`David Yarowsky. 1995. Unsupervised word sense disambiguation rivaling
    supervised methods. In Proceedings of the 33rd annual meeting on
    Association for Computational Linguistics (ACL '95). Association for
    Computational Linguistics, Stroudsburg, PA, USA, 189-196.
    <10.3115/981658.981684>`

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import datasets
    >>> from sklearn.semi_supervised import SelfTrainingClassifier
    >>> from sklearn.svm import SVC
    >>> rng = np.random.RandomState(42)
    >>> iris = datasets.load_iris()
    >>> random_unlabeled_points = rng.rand(iris.target.shape[0]) < 0.3
    >>> iris.target[random_unlabeled_points] = -1
    >>> svc = SVC(probability=True, gamma="auto")
    >>> self_training_model = SelfTrainingClassifier(svc)
    >>> self_training_model.fit(iris.data, iris.target)
    SelfTrainingClassifier(...)
    """

    _estimator_type = "classifier"

    _parameter_constraints: dict = {
        # We don't require `predic_proba` here to allow passing a meta-estimator
        # that only exposes `predict_proba` after fitting.
        "base_estimator": [HasMethods(["fit"])],
        "threshold": [Interval(Real, 0.0, 1.0, closed="left")],
        "criterion": [StrOptions({"threshold", "k_best"})],
        "k_best": [Interval(Integral, 1, None, closed="left")],
        "max_iter": [Interval(Integral, 0, None, closed="left"), None],
        "verbose": ["verbose"],
    }

    def __init__(
        self,
        base_estimator,
        threshold=0.75,
        criterion="threshold",
        k_best=10,
        max_iter=10,
        verbose=False,
    ):
        self.base_estimator = base_estimator
        self.threshold = threshold
        self.criterion = criterion
        self.k_best = k_best
        self.max_iter = max_iter
        self.verbose = verbose

    @_fit_context(
        # SelfTrainingClassifier.base_estimator is not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y, **params):
        """
        Fit self-training classifier using `X`, `y` as training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Array representing the data.

        y : {array-like, sparse matrix} of shape (n_samples,)
            Array representing the labels. Unlabeled samples should have the
            label -1.

        **params : dict
            Parameters to pass to the underlying estimators.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        _raise_for_params(params, self, "fit")

        # we need row slicing support for sparse matrices, but costly finiteness check
        # can be delegated to the base estimator.
        X, y = self._validate_data(
            X, y, accept_sparse=["csr", "csc", "lil", "dok"], force_all_finite=False
        )

        self.base_estimator_ = clone(self.base_estimator)

        if y.dtype.kind in ["U", "S"]:
            raise ValueError(
                "y has dtype string. If you wish to predict on "
                "string targets, use dtype object, and use -1"
                " as the label for unlabeled samples."
            )

        has_label = y != -1

        if np.all(has_label):
            warnings.warn("y contains no unlabeled samples", UserWarning)

        if self.criterion == "k_best" and (
            self.k_best > X.shape[0] - np.sum(has_label)
        ):
            warnings.warn(
                (
                    "k_best is larger than the amount of unlabeled "
                    "samples. All unlabeled samples will be labeled in "
                    "the first iteration"
                ),
                UserWarning,
            )

        if _routing_enabled():
            routed_params = process_routing(self, "fit", **params)
        else:
            routed_params = Bunch()
            routed_params.base_estimator = Bunch(fit={})

        self.transduction_ = np.copy(y)
        self.labeled_iter_ = np.full_like(y, -1)
        self.labeled_iter_[has_label] = 0

        self.n_iter_ = 0

        while not np.all(has_label) and (
            self.max_iter is None or self.n_iter_ < self.max_iter
        ):
            self.n_iter_ += 1
            self.base_estimator_.fit(
                X[safe_mask(X, has_label)],
                self.transduction_[has_label],
                **routed_params.base_estimator.fit,
            )

            # Predict on the unlabeled samples
            prob = self.base_estimator_.predict_proba(X[safe_mask(X, ~has_label)])
            pred = self.base_estimator_.classes_[np.argmax(prob, axis=1)]
            max_proba = np.max(prob, axis=1)

            # Select new labeled samples
            if self.criterion == "threshold":
                selected = max_proba > self.threshold
            else:
                n_to_select = min(self.k_best, max_proba.shape[0])
                if n_to_select == max_proba.shape[0]:
                    selected = np.ones_like(max_proba, dtype=bool)
                else:
                    # NB these are indices, not a mask
                    selected = np.argpartition(-max_proba, n_to_select)[:n_to_select]

            # Map selected indices into original array
            selected_full = np.nonzero(~has_label)[0][selected]

            # Add newly labeled confident predictions to the dataset
            self.transduction_[selected_full] = pred[selected]
            has_label[selected_full] = True
            self.labeled_iter_[selected_full] = self.n_iter_

            if selected_full.shape[0] == 0:
                # no changed labels
                self.termination_condition_ = "no_change"
                break

            if self.verbose:
                print(
                    f"End of iteration {self.n_iter_},"
                    f" added {selected_full.shape[0]} new labels."
                )

        if self.n_iter_ == self.max_iter:
            self.termination_condition_ = "max_iter"
        if np.all(has_label):
            self.termination_condition_ = "all_labeled"

        self.base_estimator_.fit(
            X[safe_mask(X, has_label)],
            self.transduction_[has_label],
            **routed_params.base_estimator.fit,
        )
        self.classes_ = self.base_estimator_.classes_
        return self

    @available_if(_estimator_has("predict"))
    def predict(self, X, **params):
        """Predict the classes of `X`.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Array representing the data.

        **params : dict of str -> object
            Parameters to pass to the underlying estimator's predict method.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        y : ndarray of shape (n_samples,)
            Array with predicted labels.
        """
        check_is_fitted(self)
        _raise_for_params(params, self, "predict")

        if _routing_enabled():
            # metadata routing is enabled.
            routed_params = process_routing(self, "predict", **params)
        else:
            routed_params = Bunch()
            routed_params.base_estimator = Bunch(predict={})

        X = self._validate_data(
            X,
            accept_sparse=True,
            force_all_finite=False,
            reset=False,
        )
        return self.base_estimator_.predict(X, **routed_params.base_estimator.predict)

    @available_if(_estimator_has("predict_proba"))
    def predict_proba(self, X, **params):
        """Predict probability for each possible outcome.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Array representing the data.

        **params : dict of str -> object
            Parameters to pass to the underlying estimator's predict_proba method.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        y : ndarray of shape (n_samples, n_features)
            Array with prediction probabilities.
        """
        check_is_fitted(self)
        _raise_for_params(params, self, "predict_proba")

        if _routing_enabled():
            # metadata routing is enabled.
            routed_params = process_routing(self, "predict_proba", **params)
        else:
            routed_params = Bunch()
            routed_params.base_estimator = Bunch(predict_proba={})

        X = self._validate_data(
            X,
            accept_sparse=True,
            force_all_finite=False,
            reset=False,
        )
        return self.base_estimator_.predict_proba(
            X, **routed_params.base_estimator.predict_proba
        )

    @available_if(_estimator_has("decision_function"))
    def decision_function(self, X, **params):
        """Call decision function of the `base_estimator`.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Array representing the data.

        **params : dict of str -> object
            Parameters to pass to the underlying estimator's decision_function method.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        y : ndarray of shape (n_samples, n_features)
            Result of the decision function of the `base_estimator`.
        """
        check_is_fitted(self)
        _raise_for_params(params, self, "decision_function")

        if _routing_enabled():
            # metadata routing is enabled.
            routed_params = process_routing(self, "decision_function", **params)
        else:
            routed_params = Bunch()
            routed_params.base_estimator = Bunch(decision_function={})

        X = self._validate_data(
            X,
            accept_sparse=True,
            force_all_finite=False,
            reset=False,
        )
        return self.base_estimator_.decision_function(
            X, **routed_params.base_estimator.decision_function
        )

    @available_if(_estimator_has("predict_log_proba"))
    def predict_log_proba(self, X, **params):
        """Predict log probability for each possible outcome.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Array representing the data.

        **params : dict of str -> object
            Parameters to pass to the underlying estimator's predict_log_proba method.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        y : ndarray of shape (n_samples, n_features)
            Array with log prediction probabilities.
        """
        check_is_fitted(self)
        _raise_for_params(params, self, "predict_log_proba")

        if _routing_enabled():
            # metadata routing is enabled.
            routed_params = process_routing(self, "predict_log_proba", **params)
        else:
            routed_params = Bunch()
            routed_params.base_estimator = Bunch(predict_log_proba={})

        X = self._validate_data(
            X,
            accept_sparse=True,
            force_all_finite=False,
            reset=False,
        )
        return self.base_estimator_.predict_log_proba(
            X, **routed_params.base_estimator.predict_log_proba
        )

    @available_if(_estimator_has("score"))
    def score(self, X, y, **params):
        """Call score on the `base_estimator`.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Array representing the data.

        y : array-like of shape (n_samples,)
            Array representing the labels.

        **params : dict of str -> object
            Parameters requested and accepted by base estimator's score function.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        score : float
            Result of calling score on the `base_estimator`.
        """
        check_is_fitted(self)
        _raise_for_params(params, self, "score")

        if _routing_enabled():
            # metadata routing is enabled.
            routed_params = process_routing(self, "score", **params)
        else:
            routed_params = Bunch()
            routed_params.base_estimator = Bunch(score={})

        X = self._validate_data(
            X,
            accept_sparse=True,
            force_all_finite=False,
            reset=False,
        )
        return self.base_estimator_.score(X, y, **routed_params.base_estimator.score)

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.5

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self.__class__.__name__)
        router.add(
            base_estimator=self.base_estimator,
            method_mapping=(
                MethodMapping()
                .add(callee="fit", caller="fit")
                .add(callee="predict", caller="predict")
                .add(callee="predict_proba", caller="predict_proba")
                .add(callee="decision_function", caller="decision_function")
                .add(callee="predict_log_proba", caller="predict_log_proba")
                .add(callee="score", caller="score")
            ),
        )
        return router
