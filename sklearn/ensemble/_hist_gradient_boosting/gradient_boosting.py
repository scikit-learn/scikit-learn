"""Fast Gradient Boosting decision trees for classification and regression."""
# Author: Nicolas Hug

from abc import ABC, abstractmethod

import numpy as np
from timeit import default_timer as time
from sklearn.base import BaseEstimator, RegressorMixin, ClassifierMixin
from sklearn.utils import check_X_y, check_random_state, check_array
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.multiclass import check_classification_targets
from sklearn.metrics import check_scoring
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from ._gradient_boosting import _update_raw_predictions
from .types import Y_DTYPE, X_DTYPE, X_BINNED_DTYPE

from .binning import _BinMapper
from .grower import TreeGrower
from .loss import _LOSSES


class BaseHistGradientBoosting(BaseEstimator, ABC):
    """Base class for histogram-based gradient boosting estimators."""

    @abstractmethod
    def __init__(self, loss, learning_rate, max_iter, max_leaf_nodes,
                 max_depth, min_samples_leaf, l2_regularization, max_bins,
                 scoring, validation_fraction, n_iter_no_change, tol, verbose,
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
        self.validation_fraction = validation_fraction
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
                "Loss {} is not supported for {}. Accepted losses: "
                "{}.".format(self.loss, self.__class__.__name__,
                             ', '.join(self._VALID_LOSSES)))

        if self.learning_rate <= 0:
            raise ValueError('learning_rate={} must '
                             'be strictly positive'.format(self.learning_rate))
        if self.max_iter < 1:
            raise ValueError('max_iter={} must not be smaller '
                             'than 1.'.format(self.max_iter))
        if self.n_iter_no_change is not None and self.n_iter_no_change < 0:
            raise ValueError('n_iter_no_change={} must be '
                             'positive.'.format(self.n_iter_no_change))
        if (self.validation_fraction is not None and
                self.validation_fraction <= 0):
            raise ValueError(
                'validation_fraction={} must be strictly '
                'positive, or None.'.format(self.validation_fraction))
        if self.tol is not None and self.tol < 0:
            raise ValueError('tol={} '
                             'must not be smaller than 0.'.format(self.tol))

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
        acc_compute_hist_time = 0.  # time spent computing histograms
        # time spent predicting X for gradient and hessians update
        acc_prediction_time = 0.
        X, y = check_X_y(X, y, dtype=[X_DTYPE])
        y = self._encode_y(y)
        rng = check_random_state(self.random_state)

        self._validate_parameters()
        self.n_features_ = X.shape[1]  # used for validation in predict()

        # we need this stateful variable to tell raw_predict() that it was
        # called from fit() (this current method), and that the data it has
        # received is pre-binned.
        # predicting is faster on pre-binned data, so we want early stopping
        # predictions to be made on pre-binned data. Unfortunately the scorer_
        # can only call predict() or predict_proba(), not raw_predict(), and
        # there's no way to tell the scorer that it needs to predict binned
        # data.
        self._in_fit = True


        self.loss_ = self._get_loss()

        self.do_early_stopping_ = (self.n_iter_no_change is not None and
                                   self.n_iter_no_change > 0)

        # create validation data if needed
        self._use_validation_data = self.validation_fraction is not None
        if self.do_early_stopping_ and self._use_validation_data:
            # stratify for classification
            stratify = y if hasattr(self.loss_, 'predict_proba') else None

            X_train, X_val, y_train, y_val = train_test_split(
                X, y, test_size=self.validation_fraction, stratify=stratify,
                random_state=rng)
        else:
            X_train, y_train = X, y
            X_val, y_val = None, None

        # Bin the data
        self.bin_mapper_ = _BinMapper(max_bins=self.max_bins, random_state=rng)
        X_binned_train = self._bin_data(X_train, rng, is_training_data=True)
        if X_val is not None:
            X_binned_val = self._bin_data(X_val, rng, is_training_data=False)
        else:
            X_binned_val = None

        if self.verbose:
            print("Fitting gradient boosted rounds:")

        # initialize raw_predictions: those are the accumulated values
        # predicted by the trees for the training data. raw_predictions has
        # shape (n_trees_per_iteration, n_samples) where
        # n_trees_per_iterations is n_classes in multiclass classification,
        # else 1.
        n_samples = X_binned_train.shape[0]
        self._baseline_prediction = self.loss_.get_baseline_prediction(
            y_train, self.n_trees_per_iteration_
        )
        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=self._baseline_prediction.dtype
        )
        raw_predictions += self._baseline_prediction

        # initialize gradients and hessians (empty arrays).
        # shape = (n_trees_per_iteration, n_samples).
        gradients, hessians = self.loss_.init_gradients_and_hessians(
            n_samples=n_samples,
            prediction_dim=self.n_trees_per_iteration_
        )

        # predictors is a matrix (list of lists) of TreePredictor objects
        # with shape (n_iter_, n_trees_per_iteration)
        self._predictors = predictors = []

        # Initialize structures and attributes related to early stopping
        self.scorer_ = None  # set if scoring != loss
        raw_predictions_val = None  # set if scoring == loss and use val
        self.train_score_ = []
        self.validation_score_ = []
        if self.do_early_stopping_:
            # populate train_score and validation_score with the predictions
            # of the initial model (before the first tree)

            if self.scoring == 'loss':
                # we're going to compute scoring w.r.t the loss. As losses
                # take raw predictions as input (unlike the scorers), we can
                # optimize a bit and avoid repeating computing the predictions
                # of the previous trees. We'll re-use raw_predictions (as it's
                # needed for training anyway) for evaluating the training
                # loss, and create raw_predictions_val for storing the
                # raw predictions of the validation data.

                if self._use_validation_data:
                    raw_predictions_val = np.zeros(
                        shape=(self.n_trees_per_iteration_,
                               X_binned_val.shape[0]),
                        dtype=self._baseline_prediction.dtype
                    )

                    raw_predictions_val += self._baseline_prediction

                self._check_early_stopping_loss(raw_predictions, y_train,
                                                raw_predictions_val, y_val)
            else:
                self.scorer_ = check_scoring(self, self.scoring)
                # scorer_ is a callable with signature (est, X, y) and calls
                # est.predict() or est.predict_proba() depending on its nature.
                # Unfortunately, each call to scorer_() will compute
                # the predictions of all the trees. So we use a subset of the
                # training set to compute train scores.
                subsample_size = 10000  # should we expose this parameter?
                indices = np.arange(X_binned_train.shape[0])
                if X_binned_train.shape[0] > subsample_size:
                    # TODO: not critical but stratify using resample()
                    indices = rng.choice(indices, subsample_size,
                                         replace=False)
                X_binned_small_train = X_binned_train[indices]
                y_small_train = y_train[indices]
                # Predicting is faster on C-contiguous arrays.
                X_binned_small_train = np.ascontiguousarray(
                    X_binned_small_train)

                self._check_early_stopping_scorer(
                    X_binned_small_train, y_small_train,
                    X_binned_val, y_val,
                )

        for iteration in range(self.max_iter):

            if self.verbose:
                iteration_start_time = time()
                print("[{}/{}] ".format(iteration + 1, self.max_iter),
                      end='', flush=True)

            # Update gradients and hessians, inplace
            self.loss_.update_gradients_and_hessians(gradients, hessians,
                                                     y_train, raw_predictions)

            # Append a list since there may be more than 1 predictor per iter
            predictors.append([])

            # Build `n_trees_per_iteration` trees.
            for k in range(self.n_trees_per_iteration_):

                grower = TreeGrower(
                    X_binned_train, gradients[k, :], hessians[k, :],
                    max_bins=self.max_bins,
                    actual_n_bins=self.bin_mapper_.actual_n_bins_,
                    max_leaf_nodes=self.max_leaf_nodes,
                    max_depth=self.max_depth,
                    min_samples_leaf=self.min_samples_leaf,
                    l2_regularization=self.l2_regularization,
                    shrinkage=self.learning_rate)
                grower.grow()

                acc_apply_split_time += grower.total_apply_split_time
                acc_find_split_time += grower.total_find_split_time
                acc_compute_hist_time += grower.total_compute_hist_time

                predictor = grower.make_predictor(
                    bin_thresholds=self.bin_mapper_.bin_thresholds_
                )
                predictors[-1].append(predictor)

                # Update raw_predictions with the predictions of the newly
                # created tree.
                tic_pred = time()
                _update_raw_predictions(raw_predictions[k, :], grower)
                toc_pred = time()
                acc_prediction_time += toc_pred - tic_pred

            should_early_stop = False
            if self.do_early_stopping_:
                if self.scoring == 'loss':
                    # Update raw_predictions_val with the newest tree(s)
                    if self._use_validation_data:
                        for k, pred in enumerate(self._predictors[-1]):
                            raw_predictions_val[k, :] += (
                                pred.predict_binned(X_binned_val))

                    should_early_stop = self._check_early_stopping_loss(
                        raw_predictions, y_train,
                        raw_predictions_val, y_val
                    )

                else:
                    should_early_stop = self._check_early_stopping_scorer(
                        X_binned_small_train, y_small_train,
                        X_binned_val, y_val,
                    )

            if self.verbose:
                self._print_iteration_stats(iteration_start_time)

            # maybe we could also early stop if all the trees are stumps?
            if should_early_stop:
                break

        if self.verbose:
            duration = time() - fit_start_time
            n_total_leaves = sum(
                predictor.get_n_leaf_nodes()
                for predictors_at_ith_iteration in self._predictors
                for predictor in predictors_at_ith_iteration
            )
            n_predictors = sum(
                len(predictors_at_ith_iteration)
                for predictors_at_ith_iteration in self._predictors)
            print("Fit {} trees in {:.3f} s, ({} total leaves)".format(
                n_predictors, duration, n_total_leaves))
            print("{:<32} {:.3f}s".format('Time spent computing histograms:',
                                          acc_compute_hist_time))
            print("{:<32} {:.3f}s".format('Time spent finding best splits:',
                                          acc_find_split_time))
            print("{:<32} {:.3f}s".format('Time spent applying splits:',
                                          acc_apply_split_time))
            print("{:<32} {:.3f}s".format('Time spent predicting:',
                                          acc_prediction_time))

        self.train_score_ = np.asarray(self.train_score_)
        self.validation_score_ = np.asarray(self.validation_score_)
        del self._in_fit  # hard delete so we're sure it can't be used anymore
        return self

    def _check_early_stopping_scorer(self, X_binned_small_train, y_small_train,
                                     X_binned_val, y_val):
        """Check if fitting should be early-stopped based on scorer.

        Scores are computed on validation data or on training data.
        """

        self.train_score_.append(
            self.scorer_(self, X_binned_small_train, y_small_train)
        )

        if self._use_validation_data:
            self.validation_score_.append(
                self.scorer_(self, X_binned_val, y_val)
            )
            return self._should_stop(self.validation_score_)
        else:
            return self._should_stop(self.train_score_)

    def _check_early_stopping_loss(self,
                                   raw_predictions,
                                   y_train,
                                   raw_predictions_val,
                                   y_val):
        """Check if fitting should be early-stopped based on loss.

        Scores are computed on validation data or on training data.
        """

        self.train_score_.append(
            -self.loss_(y_train, raw_predictions)
        )

        if self._use_validation_data:
            self.validation_score_.append(
                -self.loss_(y_val, raw_predictions_val)
            )
            return self._should_stop(self.validation_score_)
        else:
            return self._should_stop(self.train_score_)

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

    def _bin_data(self, X, rng, is_training_data):
        """Bin data X.

        If is_training_data, then set the bin_mapper_ attribute.
        Else, the binned data is converted to a C-contiguous array.
        """

        description = 'training' if is_training_data else 'validation'
        if self.verbose:
            print("Binning {:.3f} GB of {} data: ".format(
                X.nbytes / 1e9, description), end="", flush=True)
        tic = time()
        if is_training_data:
            X_binned = self.bin_mapper_.fit_transform(X)  # F-aligned array
        else:
            X_binned = self.bin_mapper_.transform(X)  # F-aligned array
            # We convert the array to C-contiguous since predicting is faster
            # with this layout (training is faster on F-arrays though)
            X_binned = np.ascontiguousarray(X_binned)
        toc = time()
        if self.verbose:
            duration = toc - tic
            print("{:.3f} s".format(duration))

        return X_binned

    def _print_iteration_stats(self, iteration_start_time):
        """Print info about the current fitting iteration."""
        log_msg = ''

        predictors_of_ith_iteration = [
            predictors_list for predictors_list in self._predictors[-1]
            if predictors_list
        ]
        n_trees = len(predictors_of_ith_iteration)
        max_depth = max(predictor.get_max_depth()
                        for predictor in predictors_of_ith_iteration)
        n_leaves = sum(predictor.get_n_leaf_nodes()
                       for predictor in predictors_of_ith_iteration)

        if n_trees == 1:
            log_msg += ("{} tree, {} leaves, ".format(n_trees, n_leaves))
        else:
            log_msg += ("{} trees, {} leaves ".format(n_trees, n_leaves))
            log_msg += ("({} on avg), ".format(int(n_leaves / n_trees)))

        log_msg += "max depth = {}, ".format(max_depth)

        if self.do_early_stopping_:
            if self.scoring == 'loss':
                factor = -1  # score_ arrays contain the negative loss
                name = 'loss'
            else:
                factor = 1
                name = 'score'
            log_msg += "train {}: {:.5f}, ".format(name, factor *
                                                   self.train_score_[-1])
            if self._use_validation_data:
                log_msg += "val {}: {:.5f}, ".format(
                    name, factor * self.validation_score_[-1])

        iteration_time = time() - iteration_start_time
        log_msg += "in {:0.3f}s".format(iteration_time)

        print(log_msg)

    def _raw_predict(self, X):
        """Return the sum of the leaves values over all predictors.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples.

        Returns
        -------
        raw_predictions : array, shape (n_samples * n_trees_per_iteration,)
            The raw predicted values.
        """
        X = check_array(X, dtype=[X_DTYPE, X_BINNED_DTYPE])
        check_is_fitted(self, '_predictors')
        if X.shape[1] != self.n_features_:
            raise ValueError(
                'X has {} features but this estimator was trained with '
                '{} features.'.format(X.shape[1], self.n_features_)
            )
        is_binned = getattr(self, '_in_fit', False)
        n_samples = X.shape[0]
        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=self._baseline_prediction.dtype
        )
        raw_predictions += self._baseline_prediction
        for predictors_of_ith_iteration in self._predictors:
            for k, predictor in enumerate(predictors_of_ith_iteration):
                predict = (predictor.predict_binned if is_binned
                           else predictor.predict)
                raw_predictions[k, :] += predict(X)

        return raw_predictions

    @abstractmethod
    def _get_loss(self):
        pass

    @abstractmethod
    def _encode_y(self, y=None):
        pass

    @property
    def n_iter_(self):
        check_is_fitted(self, '_predictors')
        return len(self._predictors)


class HistGradientBoostingRegressor(BaseHistGradientBoosting, RegressorMixin):
    """Histogram-based Gradient Boosting Regression Tree.

    This estimator is much faster than
    :class:`GradientBoostingRegressor<sklearn.ensemble.GradientBoostingRegressor>`
    for big datasets (n_samples >= 10 000). The input data ``X`` is pre-binned
    into integer-valued bins, which considerably reduces the number of
    splitting points to consider, and allows the algorithm to leverage
    integer-based data structures. For small sample sizes,
    :class:`GradientBoostingRegressor<sklearn.ensemble.GradientBoostingRegressor>`
    might be preferred since binning may lead to split points that are too
    approximate in this setting.

    This implementation is inspired by
    `LightGBM <https://github.com/Microsoft/LightGBM>`_.

    .. note::

      This estimator is still **experimental** for now: the predictions
      and the API might change without any deprecation cycle. To use it,
      you need to explicitly import ``enable_hist_gradient_boosting``::

        >>> # explicitly require this experimental feature
        >>> from sklearn.experimental import enable_hist_gradient_boosting  # noqa
        >>> # now you can import normally from ensemble
        >>> from sklearn.ensemble import HistGradientBoostingClassifier


    Parameters
    ----------
    loss : {'least_squares'}, optional (default='least_squares')
        The loss function to use in the boosting process. Note that the
        "least squares" loss actually implements an "half least squares loss"
        to simplify the computation of the gradient.
    learning_rate : float, optional (default=0.1)
        The learning rate, also known as *shrinkage*. This is used as a
        multiplicative factor for the leaves values. Use ``1`` for no
        shrinkage.
    max_iter : int, optional (default=100)
        The maximum number of iterations of the boosting process, i.e. the
        maximum number of trees.
    max_leaf_nodes : int or None, optional (default=31)
        The maximum number of leaves for each tree. Must be strictly greater
        than 1. If None, there is no maximum limit.
    max_depth : int or None, optional (default=None)
        The maximum depth of each tree. The depth of a tree is the number of
        nodes to go from the root to the deepest leaf. Must be strictly greater
        than 1. Depth isn't constrained by default.
    min_samples_leaf : int, optional (default=20)
        The minimum number of samples per leaf. For small datasets with less
        than a few hundred samples, it is recommended to lower this value
        since only very shallow trees would be built.
    l2_regularization : float, optional (default=0)
        The L2 regularization parameter. Use ``0`` for no regularization
        (default).
    max_bins : int, optional (default=256)
        The maximum number of bins to use. Before training, each feature of
        the input array ``X`` is binned into at most ``max_bins`` bins, which
        allows for a much faster training stage. Features with a small
        number of unique values may use less than ``max_bins`` bins. Must be no
        larger than 256.
    scoring : str or callable or None, optional (default=None)
        Scoring parameter to use for early stopping. It can be a single
        string (see :ref:`scoring_parameter`) or a callable (see
        :ref:`scoring`). If None, the estimator's default scorer is used. If
        ``scoring='loss'``, early stopping is checked w.r.t the loss value.
        Only used if ``n_iter_no_change`` is not None.
    validation_fraction : int or float or None, optional (default=0.1)
        Proportion (or absolute size) of training data to set aside as
        validation data for early stopping. If None, early stopping is done on
        the training data. Only used if ``n_iter_no_change`` is not None.
    n_iter_no_change : int or None, optional (default=None)
        Used to determine when to "early stop". The fitting process is
        stopped when none of the last ``n_iter_no_change`` scores are better
        than the ``n_iter_no_change - 1``th-to-last one, up to some
        tolerance. If None or 0, no early-stopping is done.
    tol : float or None, optional (default=1e-7)
        The absolute tolerance to use when comparing scores during early
        stopping. The higher the tolerance, the more likely we are to early
        stop: higher tolerance means that it will be harder for subsequent
        iterations to be considered an improvement upon the reference score.
    verbose: int, optional (default=0)
        The verbosity level. If not zero, print some information about the
        fitting process.
    random_state : int, np.random.RandomStateInstance or None, \
        optional (default=None)
        Pseudo-random number generator to control the subsampling in the
        binning process, and the train/validation data split if early stopping
        is enabled. See :term:`random_state`.

    Attributes
    ----------
    n_iter_ : int
        The number of iterations as selected by early stopping (if
        n_iter_no_change is not None). Otherwise it corresponds to max_iter.
    n_trees_per_iteration_ : int
        The number of tree that are built at each iteration. For regressors,
        this is always 1.
    train_score_ : ndarray, shape (max_iter + 1,)
        The scores at each iteration on the training data. The first entry
        is the score of the ensemble before the first iteration. Scores are
        computed according to the ``scoring`` parameter. If ``scoring`` is
        not 'loss', scores are computed on a subset of at most 10 000
        samples. Empty if no early stopping.
    validation_score_ : ndarray, shape (max_iter + 1,)
        The scores at each iteration on the held-out validation data. The
        first entry is the score of the ensemble before the first iteration.
        Scores are computed according to the ``scoring`` parameter. Empty if
        no early stopping or if ``validation_fraction`` is None.

    Examples
    --------
    >>> # To use this experimental feature, we need to explicitly ask for it:
    >>> from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    >>> from sklearn.ensemble import HistGradientBoostingRegressor
    >>> from sklearn.datasets import load_boston
    >>> X, y = load_boston(return_X_y=True)
    >>> est = HistGradientBoostingRegressor().fit(X, y)
    >>> est.score(X, y)
    0.98...
    """

    _VALID_LOSSES = ('least_squares',)

    def __init__(self, loss='least_squares', learning_rate=0.1,
                 max_iter=100, max_leaf_nodes=31, max_depth=None,
                 min_samples_leaf=20, l2_regularization=0., max_bins=256,
                 scoring=None, validation_fraction=0.1, n_iter_no_change=None,
                 tol=1e-7, verbose=0, random_state=None):
        super(HistGradientBoostingRegressor, self).__init__(
            loss=loss, learning_rate=learning_rate, max_iter=max_iter,
            max_leaf_nodes=max_leaf_nodes, max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            l2_regularization=l2_regularization, max_bins=max_bins,
            scoring=scoring, validation_fraction=validation_fraction,
            n_iter_no_change=n_iter_no_change, tol=tol, verbose=verbose,
            random_state=random_state)

    def predict(self, X):
        """Predict values for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The predicted values.
        """
        # Return raw predictions after converting shape
        # (n_samples, 1) to (n_samples,)
        return self._raw_predict(X).ravel()

    def _encode_y(self, y):
        # Just convert y to the expected dtype
        self.n_trees_per_iteration_ = 1
        y = y.astype(Y_DTYPE, copy=False)
        return y

    def _get_loss(self):
        return _LOSSES[self.loss]()


class HistGradientBoostingClassifier(BaseHistGradientBoosting,
                                     ClassifierMixin):
    """Histogram-based Gradient Boosting Classification Tree.

    This estimator is much faster than
    :class:`GradientBoostingClassifier<sklearn.ensemble.GradientBoostingClassifier>`
    for big datasets (n_samples >= 10 000). The input data ``X`` is pre-binned
    into integer-valued bins, which considerably reduces the number of
    splitting points to consider, and allows the algorithm to leverage
    integer-based data structures. For small sample sizes,
    :class:`GradientBoostingClassifier<sklearn.ensemble.GradientBoostingClassifier>`
    might be preferred since binning may lead to split points that are too
    approximate in this setting.

    This implementation is inspired by
    `LightGBM <https://github.com/Microsoft/LightGBM>`_.

    .. note::

      This estimator is still **experimental** for now: the predictions
      and the API might change without any deprecation cycle. To use it,
      you need to explicitly import ``enable_hist_gradient_boosting``::

        >>> # explicitly require this experimental feature
        >>> from sklearn.experimental import enable_hist_gradient_boosting  # noqa
        >>> # now you can import normally from ensemble
        >>> from sklearn.ensemble import HistGradientBoostingClassifier

    Parameters
    ----------
    loss : {'auto', 'binary_crossentropy', 'categorical_crossentropy'}, \
            optional (default='auto')
        The loss function to use in the boosting process. 'binary_crossentropy'
        (also known as logistic loss) is used for binary classification and
        generalizes to 'categorical_crossentropy' for multiclass
        classification. 'auto' will automatically choose either loss depending
        on the nature of the problem.
    learning_rate : float, optional (default=1)
        The learning rate, also known as *shrinkage*. This is used as a
        multiplicative factor for the leaves values. Use ``1`` for no
        shrinkage.
    max_iter : int, optional (default=100)
        The maximum number of iterations of the boosting process, i.e. the
        maximum number of trees for binary classification. For multiclass
        classification, `n_classes` trees per iteration are built.
    max_leaf_nodes : int or None, optional (default=31)
        The maximum number of leaves for each tree. Must be strictly greater
        than 1. If None, there is no maximum limit.
    max_depth : int or None, optional (default=None)
        The maximum depth of each tree. The depth of a tree is the number of
        nodes to go from the root to the deepest leaf. Must be strictly greater
        than 1. Depth isn't constrained by default.
    min_samples_leaf : int, optional (default=20)
        The minimum number of samples per leaf. For small datasets with less
        than a few hundred samples, it is recommended to lower this value
        since only very shallow trees would be built.
    l2_regularization : float, optional (default=0)
        The L2 regularization parameter. Use 0 for no regularization.
    max_bins : int, optional (default=256)
        The maximum number of bins to use. Before training, each feature of
        the input array ``X`` is binned into at most ``max_bins`` bins, which
        allows for a much faster training stage. Features with a small
        number of unique values may use less than ``max_bins`` bins. Must be no
        larger than 256.
    scoring : str or callable or None, optional (default=None)
        Scoring parameter to use for early stopping. It can be a single
        string (see :ref:`scoring_parameter`) or a callable (see
        :ref:`scoring`). If None, the estimator's default scorer
        is used. If ``scoring='loss'``, early stopping is checked
        w.r.t the loss value. Only used if ``n_iter_no_change`` is not None.
    validation_fraction : int or float or None, optional (default=0.1)
        Proportion (or absolute size) of training data to set aside as
        validation data for early stopping. If None, early stopping is done on
        the training data.
    n_iter_no_change : int or None, optional (default=None)
        Used to determine when to "early stop". The fitting process is
        stopped when none of the last ``n_iter_no_change`` scores are better
        than the ``n_iter_no_change - 1``th-to-last one, up to some
        tolerance. If None or 0, no early-stopping is done.
    tol : float or None, optional (default=1e-7)
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
        is enabled. See :term:`random_state`.

    Attributes
    ----------
    n_iter_ : int
        The number of estimators as selected by early stopping (if
        n_iter_no_change is not None). Otherwise it corresponds to max_iter.
    n_trees_per_iteration_ : int
        The number of tree that are built at each iteration. This is equal to 1
        for binary classification, and to ``n_classes`` for multiclass
        classification.
    train_score_ : ndarray, shape (max_iter + 1,)
        The scores at each iteration on the training data. The first entry
        is the score of the ensemble before the first iteration. Scores are
        computed according to the ``scoring`` parameter. If ``scoring`` is
        not 'loss', scores are computed on a subset of at most 10 000
        samples. Empty if no early stopping.
    validation_score_ : ndarray, shape (max_iter + 1,)
        The scores at each iteration on the held-out validation data. The
        first entry is the score of the ensemble before the first iteration.
        Scores are computed according to the ``scoring`` parameter. Empty if
        no early stopping or if ``validation_fraction`` is None.

    Examples
    --------
    >>> # To use this experimental feature, we need to explicitly ask for it:
    >>> from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    >>> from sklearn.ensemble import HistGradientBoostingRegressor
    >>> from sklearn.datasets import load_iris
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = HistGradientBoostingClassifier().fit(X, y)
    >>> clf.score(X, y)
    1.0
    """

    _VALID_LOSSES = ('binary_crossentropy', 'categorical_crossentropy',
                     'auto')

    def __init__(self, loss='auto', learning_rate=0.1, max_iter=100,
                 max_leaf_nodes=31, max_depth=None, min_samples_leaf=20,
                 l2_regularization=0., max_bins=256, scoring=None,
                 validation_fraction=0.1, n_iter_no_change=None, tol=1e-7,
                 verbose=0, random_state=None):
        super(HistGradientBoostingClassifier, self).__init__(
            loss=loss, learning_rate=learning_rate, max_iter=max_iter,
            max_leaf_nodes=max_leaf_nodes, max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            l2_regularization=l2_regularization, max_bins=max_bins,
            scoring=scoring, validation_fraction=validation_fraction,
            n_iter_no_change=n_iter_no_change, tol=tol, verbose=verbose,
            random_state=random_state)

    def predict(self, X):
        """Predict classes for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The predicted classes.
        """
        # TODO: This could be done in parallel
        encoded_classes = np.argmax(self.predict_proba(X), axis=1)
        return self.classes_[encoded_classes]

    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        p : ndarray, shape (n_samples, n_classes)
            The class probabilities of the input samples.
        """
        raw_predictions = self._raw_predict(X)
        return self.loss_.predict_proba(raw_predictions)

    def decision_function(self, X):
        """Compute the decision function of X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        decision : ndarray, shape (n_samples,) or \
                (n_samples, n_trees_per_iteration)
            The raw predicted values (i.e. the sum of the trees leaves) for
            each sample. n_trees_per_iteration is equal to the number of
            classes in multiclass classification.
        """
        decision = self._raw_predict(X)
        if decision.shape[0] == 1:
            decision = decision.ravel()
        return decision.T

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
