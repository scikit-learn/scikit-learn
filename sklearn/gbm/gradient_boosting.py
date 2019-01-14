"""
Gradient Boosting decision trees for classification and regression.
"""
from abc import ABC, abstractmethod

import numpy as np
from time import time
from sklearn.base import BaseEstimator, RegressorMixin, ClassifierMixin
from sklearn.utils import check_X_y, check_random_state, check_array
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.multiclass import check_classification_targets
from sklearn.metrics import check_scoring
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from ._gradient_boosting import _update_raw_predictions
from .types import Y_DTYPE, X_DTYPE

from .binning import BinMapper
from .grower import TreeGrower
from .loss import _LOSSES


class BaseGradientBoostingMachine(BaseEstimator, ABC):
    """Base class for gradient boosting estimators."""

    @abstractmethod
    def __init__(self, loss, learning_rate, max_iter, max_leaf_nodes,
                 max_depth, min_samples_leaf, l2_regularization, max_bins,
                 scoring, validation_split, n_iter_no_change, tol, verbose,
                 random_state):
        self.loss = loss
        self.learning_rate = learning_rate
        self.max_iter = max_iter
        self.max_leaf_nodes = max_leaf_nodes
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.l2_regularization = l2_regularization
        self.max_bins = max_bins
        self.n_iter_no_change = n_iter_no_change
        self.validation_split = validation_split
        self.scoring = scoring
        self.tol = tol
        self.verbose = verbose
        self.random_state = random_state

    def _validate_parameters(self):
        """Validate parameters passed to __init__.

        The parameters that are directly passed to the grower are checked in
        TreeGrower."""

        if self.loss not in self._VALID_LOSSES:
            raise ValueError(
                "Loss {} is not supported for {}. Accepted losses"
                "are {}.".format(self.loss, self.__class__.__name__,
                                 ', '.join(self._VALID_LOSSES)))

        if self.learning_rate <= 0:
            raise ValueError(f'learning_rate={self.learning_rate} must '
                             f'be strictly positive')
        if self.max_iter < 1:
            raise ValueError(f'max_iter={self.max_iter} must '
                             f'not be smaller than 1.')
        if self.n_iter_no_change is not None and self.n_iter_no_change < 0:
            raise ValueError(f'n_iter_no_change={self.n_iter_no_change} '
                             f'must be positive.')
        if self.validation_split is not None and self.validation_split <= 0:
            raise ValueError(f'validation_split={self.validation_split} '
                             f'must be strictly positive, or None.')
        if self.tol is not None and self.tol < 0:
            raise ValueError(f'tol={self.tol} '
                             f'must not be smaller than 0.')

    def fit(self, X, y):
        """Fit the gradient boosting model.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples.

        y : array-like, shape=(n_samples,)
            Target values.

        Returns
        -------
        self : object
        """

        fit_start_time = time()
        acc_find_split_time = 0.  # time spent finding the best splits
        acc_apply_split_time = 0.  # time spent splitting nodes
        # time spent predicting X for gradient and hessians update
        acc_prediction_time = 0.
        # TODO: add support for mixed-typed (numerical + categorical) data
        # TODO: add support for missing data
        # TODO: add support for pre-binned data (pass-through)?
        X, y = check_X_y(X, y, dtype=[X_DTYPE])
        y = self._encode_y(y)
        rng = check_random_state(self.random_state)

        self._validate_parameters()
        self.n_features_ = X.shape[1]  # used for validation in predict()

        if self.verbose:
            print(f"Binning {X.nbytes / 1e9:.3f} GB of data: ", end="",
                  flush=True)
        tic = time()
        self.bin_mapper_ = BinMapper(max_bins=self.max_bins, random_state=rng)
        X_binned = self.bin_mapper_.fit_transform(X)
        toc = time()
        if self.verbose:
            duration = toc - tic
            troughput = X.nbytes / duration
            print(f"{duration:.3f} s ({troughput / 1e6:.3f} MB/s)")

        self.loss_ = self._get_loss()

        self.do_early_stopping_ = (self.n_iter_no_change is not None and
                                   self.n_iter_no_change > 0)

        if self.do_early_stopping_ and self.validation_split is not None:
            # stratify for classification
            stratify = y if hasattr(self.loss_, 'predict_proba') else None

            X_binned_train, X_binned_val, y_train, y_val = train_test_split(
                X_binned, y, test_size=self.validation_split,
                stratify=stratify, random_state=rng)
            if X_binned_train.size == 0 or X_binned_val.size == 0:
                raise ValueError(
                    f'Not enough data (n_samples={X_binned.shape[0]}) to '
                    f'perform early stopping with validation_split='
                    f'{self.validation_split}. Use more training data or '
                    f'adjust validation_split.'
                )
            # Predicting is faster of C-contiguous arrays, training is faster
            # on Fortran arrays.
            X_binned_val = np.ascontiguousarray(X_binned_val)
            X_binned_train = np.asfortranarray(X_binned_train)
        else:
            X_binned_train, y_train = X_binned, y
            X_binned_val, y_val = None, None

        # Subsample the training set for score-based monitoring.
        if self.do_early_stopping_:
            subsample_size = 10000
            indices = np.arange(X_binned_train.shape[0])
            if X_binned_train.shape[0] > subsample_size:
                indices = rng.choice(indices, subsample_size)
            X_binned_small_train = X_binned_train[indices]
            y_small_train = y_train[indices]
            # Predicting is faster of C-contiguous arrays.
            X_binned_small_train = np.ascontiguousarray(X_binned_small_train)

        if self.verbose:
            print("Fitting gradient boosted rounds:")

        n_samples = X_binned_train.shape[0]
        self.baseline_prediction_ = self.loss_.get_baseline_prediction(
            y_train, self.n_trees_per_iteration_)
        # raw_predictions are the accumulated values predicted by the trees
        # for the training data.
        raw_predictions = np.zeros(
            shape=(n_samples, self.n_trees_per_iteration_),
            dtype=self.baseline_prediction_.dtype
        )
        raw_predictions += self.baseline_prediction_

        # gradients and hessians are 1D arrays of size
        # n_samples * n_trees_per_iteration
        gradients, hessians = self.loss_.init_gradients_and_hessians(
            n_samples=n_samples,
            prediction_dim=self.n_trees_per_iteration_
        )

        # predictors_ is a matrix of TreePredictor objects with shape
        # (n_iter_, n_trees_per_iteration)
        self.predictors_ = predictors = []

        # scorer_ is a callable with signature (est, X, y) and calls
        # est.predict() or est.predict_proba() depending on its nature.
        self.scorer_ = check_scoring(self, self.scoring)
        self.train_scores_ = []
        self.validation_scores_ = []
        if self.do_early_stopping_:
            # Add predictions of the initial model (before the first tree)
            self.train_scores_.append(
                self._get_scores(X_binned_train, y_train))

            if self.validation_split is not None:
                self.validation_scores_.append(
                    self._get_scores(X_binned_val, y_val))

        for iteration in range(self.max_iter):

            if self.verbose:
                iteration_start_time = time()
                print(f"[{iteration + 1}/{self.max_iter}] ", end='',
                      flush=True)

            # Update gradients and hessians, inplace
            self.loss_.update_gradients_and_hessians(gradients, hessians,
                                                     y_train, raw_predictions)

            predictors.append([])

            # Build `n_trees_per_iteration` trees.
            for k, (gradients_at_k, hessians_at_k) in enumerate(zip(
                    np.array_split(gradients, self.n_trees_per_iteration_),
                    np.array_split(hessians, self.n_trees_per_iteration_))):
                # the xxxx_at_k arrays are **views** on the original arrays.
                # Note that for binary classif and regressions,
                # n_trees_per_iteration is 1 and xxxx_at_k is equivalent to the
                # whole array.

                grower = TreeGrower(
                    X_binned_train, gradients_at_k, hessians_at_k,
                    max_bins=self.max_bins,
                    n_bins_per_feature=self.bin_mapper_.n_bins_per_feature_,
                    max_leaf_nodes=self.max_leaf_nodes,
                    max_depth=self.max_depth,
                    min_samples_leaf=self.min_samples_leaf,
                    l2_regularization=self.l2_regularization,
                    shrinkage=self.learning_rate)
                grower.grow()

                acc_apply_split_time += grower.total_apply_split_time
                acc_find_split_time += grower.total_find_split_time

                predictor = grower.make_predictor(
                    bin_thresholds=self.bin_mapper_.bin_thresholds_)
                predictors[-1].append(predictor)

                tic_pred = time()
                _update_raw_predictions(raw_predictions[:, k], grower)
                toc_pred = time()
                acc_prediction_time += toc_pred - tic_pred

            should_early_stop = False
            if self.do_early_stopping_:
                should_early_stop = self._check_early_stopping(
                    X_binned_small_train, y_small_train,
                    X_binned_val, y_val)

            if self.verbose:
                self._print_iteration_stats(iteration_start_time)

            if should_early_stop:
                break

        if self.verbose:
            duration = time() - fit_start_time
            n_total_leaves = sum(
                predictor.get_n_leaf_nodes()
                for predictors_at_ith_iteration in self.predictors_
                for predictor in predictors_at_ith_iteration)
            n_predictors = sum(
                len(predictors_at_ith_iteration)
                for predictors_at_ith_iteration in self.predictors_)
            print(f"Fit {n_predictors} trees in {duration:.3f} s, "
                  f"({n_total_leaves} total leaves)")
            print(f"{'Time spent finding best splits:':<32} "
                  f"{acc_find_split_time:.3f}s")
            print(f"{'Time spent applying splits:':<32} "
                  f"{acc_apply_split_time:.3f}s")
            print(f"{'Time spent predicting:':<32} "
                  f"{acc_prediction_time:.3f}s")

        self.train_scores_ = np.asarray(self.train_scores_)
        self.validation_scores_ = np.asarray(self.validation_scores_)
        return self

    def _check_early_stopping(self, X_binned_train, y_train,
                              X_binned_val, y_val):
        """Check if fitting should be early-stopped.

        Scores are computed on validation data or on training data.
        """

        self.train_scores_.append(
            self._get_scores(X_binned_train, y_train))

        if self.validation_split is not None:
            self.validation_scores_.append(
                self._get_scores(X_binned_val, y_val))
            return self._should_stop(self.validation_scores_)

        return self._should_stop(self.train_scores_)

    def _should_stop(self, scores):
        """
        Return True (do early stopping) if the last n scores aren't better
        than the (n-1)th-to-last score, up to some tolerance.
        """
        reference_position = self.n_iter_no_change + 1
        if len(scores) < reference_position:
            return False

        # A higher score is always better. Higher tol means that it will be
        # harder for subsequent iteration to be considered an improvement upon
        # the reference score, and therefore it is more likely to early stop
        # because of the lack of significant improvement.
        tol = 0 if self.tol is None else self.tol
        reference_score = scores[-reference_position] + tol
        recent_scores = scores[-reference_position + 1:]
        recent_improvements = [score > reference_score
                               for score in recent_scores]
        return not any(recent_improvements)

    def _get_scores(self, X, y):
        """Compute scores on data X with target y.

        Scores are either computed with a scorer if scoring parameter is not
        None, else with the loss. As higher is always better, we return
        -loss_value.
        """
        if self.scoring is not None:
            return self.scorer_(self, X, y)

        # Else, use loss
        raw_predictions = self._raw_predict(X)
        return -self.loss_(y, raw_predictions)

    def _print_iteration_stats(self, iteration_start_time):
        """Print info about the current fitting iteration."""
        log_msg = ''

        predictors_of_ith_iteration = [
            predictors_list for predictors_list in self.predictors_[-1]
            if predictors_list
        ]
        n_trees = len(predictors_of_ith_iteration)
        max_depth = max(predictor.get_max_depth()
                        for predictor in predictors_of_ith_iteration)
        n_leaves = sum(predictor.get_n_leaf_nodes()
                       for predictor in predictors_of_ith_iteration)

        if n_trees == 1:
            log_msg += (f"{n_trees} tree, {n_leaves} leaves, ")
        else:
            log_msg += (f"{n_trees} trees, {n_leaves} leaves ")
            log_msg += (f"({int(n_leaves / n_trees)} on avg), ")

        log_msg += f"max depth = {max_depth}, "

        if self.do_early_stopping_:
            log_msg += f"{self.scoring} train: {self.train_scores_[-1]:.5f}, "
            if self.validation_split is not None:
                log_msg += (f"{self.scoring} val: "
                            f"{self.validation_scores_[-1]:.5f}, ")

        iteration_time = time() - iteration_start_time
        log_msg += f"in {iteration_time:0.3f}s"

        print(log_msg)

    def _raw_predict(self, X):
        """Return the sum of the leaves values over all predictors.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples. If ``X.dtype == np.uint8``, the data is assumed
            to be pre-binned.

        Returns
        -------
        raw_predictions : array, shape (n_samples * n_trees_per_iteration,)
            The raw predicted values.
        """
        X = check_array(X, dtype=X_DTYPE)
        check_is_fitted(self, 'predictors_')
        if X.shape[1] != self.n_features_:
            raise ValueError(
                f'X has {X.shape[1]} features but this estimator was '
                f'trained with {self.n_features_} features.'
            )
        n_samples = X.shape[0]
        raw_predictions = np.zeros(
            shape=(n_samples, self.n_trees_per_iteration_),
            dtype=self.baseline_prediction_.dtype
        )
        raw_predictions += self.baseline_prediction_
        for predictors_of_ith_iteration in self.predictors_:
            for k, predictor in enumerate(predictors_of_ith_iteration):
                raw_predictions[:, k] += predictor.predict(X)

        return raw_predictions

    @abstractmethod
    def _get_loss(self):
        pass

    @abstractmethod
    def _encode_y(self, y=None):
        pass

    @property
    def n_iter_(self):
        check_is_fitted(self, 'predictors_')
        return len(self.predictors_)


class GradientBoostingRegressor(BaseGradientBoostingMachine, RegressorMixin):
    """Scikit-learn compatible Gradient Boosting Tree for regression.

    Parameters
    ----------
    loss : {'least_squares'}, optional(default='least_squares')
        The loss function to use in the boosting process.
    learning_rate : float, optional(default=0.1)
        The learning rate, also known as *shrinkage*. This is used as a
        multiplicative factor for the leaves values. Use ``1`` for no
        shrinkage.
    max_iter : int, optional(default=100)
        The maximum number of iterations of the boosting process, i.e. the
        maximum number of trees.
    max_leaf_nodes : int or None, optional(default=None)
        The maximum number of leaves for each tree. If None, there is no
        maximum limit.
    max_depth : int or None, optional(default=None)
        The maximum depth of each tree. The depth of a tree is the number of
        nodes to go from the root to the deepest leaf.
    min_samples_leaf : int, optional(default=20)
        The minimum number of samples per leaf.
    l2_regularization : float, optional(default=0)
        The L2 regularization parameter. Use 0 for no regularization.
    max_bins : int, optional(default=256)
        The maximum number of bins to use. Before training, each feature of
        the input array ``X`` is binned into at most ``max_bins`` bins, which
        allows for a much faster training stage. Features with a small
        number of unique values may use less than ``max_bins`` bins. Must be no
        larger than 256.
    scoring : str or callable or None, \
        optional (default=None)
        Scoring parameter to use for early stopping (see sklearn.metrics for
        available options). If None, early stopping is check w.r.t the loss
        value.
    validation_split : int or float or None, optional(default=0.1)
        Proportion (or absolute size) of training data to set aside as
        validation data for early stopping. If None, early stopping is done on
        the training data.
    n_iter_no_change : int or None, optional (default=5)
        Used to determine when to "early stop". The fitting process is
        stopped when none of the last ``n_iter_no_change`` scores are better
        than the ``n_iter_no_change - 1``th-to-last one, up to some
        tolerance. If None or 0, no early-stopping is done.
    tol : float or None optional (default=1e-7)
        The absolute tolerance to use when comparing scores. The higher the
        tolerance, the more likely we are to early stop: higher tolerance
        means that it will be harder for subsequent iterations to be
        considered an improvement upon the reference score.
    verbose: int, optional (default=0)
        The verbosity level. If not zero, print some information about the
        fitting process.
    random_state : int, np.random.RandomStateInstance or None, \
        optional (default=None)
        Pseudo-random number generator to control the subsampling in the
        binning process, and the train/validation data split if early stopping
        is enabled. See
        `scikit-learn glossary
        <https://scikit-learn.org/stable/glossary.html#term-random-state>`_.


    Examples
    --------
    >>> from sklearn.datasets import load_boston
    >>> from pygbm import GradientBoostingRegressor
    >>> X, y = load_boston(return_X_y=True)
    >>> est = GradientBoostingRegressor().fit(X, y)
    >>> est.score(X, y)
    0.92...
    """

    _VALID_LOSSES = ('least_squares',)

    def __init__(self, loss='least_squares', learning_rate=0.1,
                 max_iter=100, max_leaf_nodes=31, max_depth=None,
                 min_samples_leaf=20, l2_regularization=0., max_bins=256,
                 scoring=None, validation_split=0.1, n_iter_no_change=5,
                 tol=1e-7, verbose=0, random_state=None):
        super(GradientBoostingRegressor, self).__init__(
            loss=loss, learning_rate=learning_rate, max_iter=max_iter,
            max_leaf_nodes=max_leaf_nodes, max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            l2_regularization=l2_regularization, max_bins=max_bins,
            scoring=scoring, validation_split=validation_split,
            n_iter_no_change=n_iter_no_change, tol=tol, verbose=verbose,
            random_state=random_state)

    def predict(self, X):
        """Predict values for X.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples. If ``X.dtype == np.uint8``, the data is assumed
            to be pre-binned.

        Returns
        -------
        y : array, shape (n_samples,)
            The predicted values.
        """
        # Return raw predictions after converting shape
        # (n_samples, 1) to (n_samples,)
        return self._raw_predict(X).ravel()

    def _encode_y(self, y):
        # Just convert y to float32
        self.n_trees_per_iteration_ = 1
        y = y.astype(Y_DTYPE, copy=False)
        return y

    def _get_loss(self):
        return _LOSSES[self.loss]()


class GradientBoostingClassifier(BaseGradientBoostingMachine, ClassifierMixin):
    """Scikit-learn compatible Gradient Boosting Tree for classification.

    Parameters
    ----------
    loss : {'auto', 'binary_crossentropy', 'categorical_crossentropy'}, \
        optional(default='auto')
        The loss function to use in the boosting process. 'binary_crossentropy'
        (also known as logistic loss) is used for binary classification and
        generalizes to 'categorical_crossentropy' for multiclass
        classification. 'auto' will automatically choose either loss depending
        on the nature of the problem.
    learning_rate : float, optional(default=1)
        The learning rate, also known as *shrinkage*. This is used as a
        multiplicative factor for the leaves values. Use ``1`` for no
        shrinkage.
    max_iter : int, optional(default=100)
        The maximum number of iterations of the boosting process, i.e. the
        maximum number of trees for binary classification. For multiclass
        classification, `n_classes` trees per iteration are built.
    max_leaf_nodes : int or None, optional(default=None)
        The maximum number of leaves for each tree. If None, there is no
        maximum limit.
    max_depth : int or None, optional(default=None)
        The maximum depth of each tree. The depth of a tree is the number of
        nodes to go from the root to the deepest leaf.
    min_samples_leaf : int, optional(default=20)
        The minimum number of samples per leaf.
    l2_regularization : float, optional(default=0)
        The L2 regularization parameter. Use 0 for no regularization.
    max_bins : int, optional(default=256)
        The maximum number of bins to use. Before training, each feature of
        the input array ``X`` is binned into at most ``max_bins`` bins, which
        allows for a much faster training stage. Features with a small
        number of unique values may use less than ``max_bins`` bins. Must be no
        larger than 256.
    scoring : str or callable or None, optional (default=None)
        Scoring parameter to use for early stopping (see sklearn.metrics for
        available options). If None, early stopping is check w.r.t the loss
        value.
    validation_split : int or float or None, optional(default=0.1)
        Proportion (or absolute size) of training data to set aside as
        validation data for early stopping. If None, early stopping is done on
        the training data.
    n_iter_no_change : int or None, optional (default=5)
        Used to determine when to "early stop". The fitting process is
        stopped when none of the last ``n_iter_no_change`` scores are better
        than the ``n_iter_no_change - 1``th-to-last one, up to some
        tolerance. If None or 0, no early-stopping is done.
    tol : float or None optional (default=1e-7)
        The absolute tolerance to use when comparing scores. The higher the
        tolerance, the more likely we are to early stop: higher tolerance
        means that it will be harder for subsequent iterations to be
        considered an improvement upon the reference score.
    verbose: int, optional(default=0)
        The verbosity level. If not zero, print some information about the
        fitting process.
    random_state : int, np.random.RandomStateInstance or None, \
        optional(default=None)
        Pseudo-random number generator to control the subsampling in the
        binning process, and the train/validation data split if early stopping
        is enabled. See `scikit-learn glossary
        <https://scikit-learn.org/stable/glossary.html#term-random-state>`_.

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from pygbm import GradientBoostingClassifier
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = GradientBoostingClassifier().fit(X, y)
    >>> clf.score(X, y)
    0.97...
    """

    _VALID_LOSSES = ('binary_crossentropy', 'categorical_crossentropy',
                     'auto')

    def __init__(self, loss='auto', learning_rate=0.1, max_iter=100,
                 max_leaf_nodes=31, max_depth=None, min_samples_leaf=20,
                 l2_regularization=0., max_bins=256, scoring=None,
                 validation_split=0.1, n_iter_no_change=5, tol=1e-7,
                 verbose=0, random_state=None):
        super(GradientBoostingClassifier, self).__init__(
            loss=loss, learning_rate=learning_rate, max_iter=max_iter,
            max_leaf_nodes=max_leaf_nodes, max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            l2_regularization=l2_regularization, max_bins=max_bins,
            scoring=scoring, validation_split=validation_split,
            n_iter_no_change=n_iter_no_change, tol=tol, verbose=verbose,
            random_state=random_state)

    def predict(self, X):
        """Predict classes for X.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples. If ``X.dtype == np.uint8``, the data is assumed
            to be pre-binned.

        Returns
        -------
        y : array, shape (n_samples,)
            The predicted classes.
        """
        # This could be done in parallel
        encoded_classes = np.argmax(self.predict_proba(X), axis=1)
        return self.classes_[encoded_classes]

    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples. If ``X.dtype == np.uint8``, the data is assumed
            to be pre-binned.

        Returns
        -------
        p : array, shape (n_samples, n_classes)
            The class probabilities of the input samples.
        """
        raw_predictions = self._raw_predict(X)
        return self.loss_.predict_proba(raw_predictions)

    def _encode_y(self, y):
        # encode classes into 0 ... n_classes - 1 and sets attributes classes_
        # and n_trees_per_iteration_
        check_classification_targets(y)

        label_encoder = LabelEncoder()
        encoded_y = label_encoder.fit_transform(y)
        self.classes_ = label_encoder.classes_
        n_classes = self.classes_.shape[0]
        # only 1 tree for binary classification. For multiclass classification,
        # we build 1 tree per class.
        self.n_trees_per_iteration_ = 1 if n_classes <= 2 else n_classes
        encoded_y = encoded_y.astype(Y_DTYPE, copy=False)
        return encoded_y

    def _get_loss(self):
        if self.loss == 'auto':
            if self.n_trees_per_iteration_ == 1:
                return _LOSSES['binary_crossentropy']()
            else:
                return _LOSSES['categorical_crossentropy']()

        return _LOSSES[self.loss]()
