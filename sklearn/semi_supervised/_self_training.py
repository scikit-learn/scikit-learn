import numpy as np

from ..base import MetaEstimatorMixin, clone
from ..utils.validation import check_X_y, check_array, check_is_fitted
from ..utils.metaestimators import if_delegate_has_method
from ..utils import safe_mask
from ..exceptions import ConvergenceWarning
from ..utils.metaestimators import _BaseComposition
import warnings

__all__ = ["SelfTrainingClassifier"]

# Authors: Oliver Rausch   <rauscho@ethz.ch>
#          Patrice Becker  <beckerp@ethz.ch>
# License: BSD 3 clause


def _validate_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if not hasattr(estimator, "predict_proba"):
        name = type(estimator).__name__
        msg = "base_classifier ({}) should implement predict_proba!"
        msg = msg.format(name)
        raise ValueError(msg)


class SelfTrainingClassifier(MetaEstimatorMixin, _BaseComposition):
    """Self-training classifier

    This class allows a given supervised classifier to function as a
    semi-supervised classifier, allowing it to learn from unlabeled data. It
    does this by iteratively predicting labels for the unlabeled data.

    Read more in the :ref:`User Guide <self_training>`.

    Parameters
    ----------
    base_classifier : estimator object
        An estimator object implementing ``fit`` and ``predict_proba``.
        Invoking the ``fit`` method will fit a clone of the passed estimator,
        which will be stored in the ``base_classifier_`` attribute.

    threshold : float, optional (default=0.75)
        The decision threshold. If the ``base_classifier`` makes a prediction
        with a ``predict_proba`` above this threshold, it will be added to the
        labeled dataset.
        Should be in [0, 1).

    max_iter : int or ``None``, optional (default=10)
        Maximum number of iterations allowed. Should be greater than or equal
        to 0. If it is ``None``, the classifier will continue to predict labels
        until all unlabeled samples have been labeled. In this case, be aware
        that the fit may never terminate.

    n_iter_no_change: int, optional (default=3)
        If not `None`, early stopping will be applied. In this case, the
        classifier will count the amount of new labels in each iteration. If
        no new labels are added to the training set for `n_iter_no_change`
        iterations, the classifier will stop fitting early.

    Attributes
    ----------
    base_classifier_ : estimator object
        The fitted estimator.

    y_labels_ : array, shape=(n_samples,)
        The labels used for the final fit of the classifier.

    y_labeled_iter_ : array, shape=(n_samples,)
        The iteration in which each sample was labeled. When a sample has
        iteration 0, the sample was already labeled in the original dataset.
        When a sample has iteration -1, the sample was not labeled in any
        iteration.

    n_iter_ : int
        The amount of iterations the classifier used during fitting.

    termination_condition_: string {'max_iter',
                                    'early_stopping',
                                    'all_labeled'}
        The reason that fitting was stopped.

        - 'max_iter': `n_iter_` reached `max_iter`.
        - 'early_stopping': no new labels were predicted for `n_iter_no_change`
          iterations.
        - 'all_labeled': all unlabeled samples were labeled before `max_iter`
          was reached.

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
    ... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    SelfTrainingClassifier(...)

    References
    ----------
    David Yarowsky. 1995. Unsupervised word sense disambiguation rivaling
    supervised methods. In Proceedings of the 33rd annual meeting on
    Association for Computational Linguistics (ACL '95). Association for
    Computational Linguistics, Stroudsburg, PA, USA, 189-196. DOI:
    https://doi.org/10.3115/981658.981684
    """
    _estimator_type = "classifier"

    def __init__(self,
                 base_classifier,
                 threshold=0.75,
                 max_iter=10,
                 n_iter_no_change=3):
        self.base_classifier = base_classifier
        self.threshold = threshold
        self.max_iter = max_iter
        self.n_iter_no_change = n_iter_no_change

    def fit(self, X, y):
        """
        Fits this ``SelfTrainingClassifier`` to a dataset.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        y : array-like, shape = (n_samples,)
            array representing the labels. Unlabeled samples should have the
            label -1.

        Returns
        -------
        self : object
            returns an instance of self.
        """
        X, y = check_X_y(X, y)

        if self.base_classifier is None:
            raise ValueError("base_classifier cannot be None!")

        self.base_classifier_ = clone(self.base_classifier)

        if self.max_iter is not None and self.max_iter < 0:
            msg = "max_iter must be >= 0 or None, got {}".format(self.max_iter)
            raise ValueError(msg)

        if not 0 <= self.threshold < 1:
            msg = "threshold must be in [0,1), got {}".format(self.threshold)
            raise ValueError(msg)

        if self.n_iter_no_change is not None and not 0 < self.n_iter_no_change:
            msg = "n_iter_no_change must be > 0, got {}"
            raise ValueError(msg.format(self.n_iter_no_change))

        if self.max_iter is not None and self.n_iter_no_change is not None \
                and self.n_iter_no_change >= self.max_iter:
            msg = "n_iter_no_change >= max_iter means " \
                "early stopping is ineffective"
            warnings.warn(msg, UserWarning)

        # XXX: y != -1 somehow doesn't work when y is an array of only strings
        # (it returns a scalar instead of an array). Is there a better way?
        has_label = np.vectorize(lambda x: x != -1)(y)

        if np.all(has_label):
            warnings.warn("y contains no unlabeled samples", UserWarning)

        self.y_labels_ = np.copy(y)
        self.y_labeled_iter_ = np.full_like(y, -1)
        self.y_labeled_iter_[has_label] = 0

        self.n_iter_ = 0
        patience = self.n_iter_no_change

        while not np.all(has_label) and (self.max_iter is None or
                                         self.n_iter_ < self.max_iter):
            self.n_iter_ += 1
            self.base_classifier_.fit(
                X[safe_mask(X, has_label)],
                self.y_labels_[safe_mask(self.y_labels_, has_label)])

            if self.n_iter_ == 1:
                _validate_estimator(self.base_classifier)

            # Predict on the unlabeled samples
            prob = self.base_classifier_.predict_proba(
                X[safe_mask(X, ~has_label)])
            pred = self.base_classifier_.classes_[np.argmax(prob, axis=1)]
            max_proba = np.max(prob, axis=1)

            # Select samples where confidence is above the threshold
            confident_labels_mask = max_proba > self.threshold

            new_labels_idx = np.flatnonzero(~has_label)[confident_labels_mask]

            # Add newly labeled confident predictions to the dataset
            self.y_labels_[new_labels_idx] = pred[confident_labels_mask]
            has_label[new_labels_idx] = True
            self.y_labeled_iter_[new_labels_idx] = self.n_iter_

            if self.n_iter_no_change is not None and (
                    new_labels_idx.shape[0] == 0):
                # no changed labels
                patience = patience - 1
                if patience == 0:
                    self.termination_condition_ = "early_stopping"
                    break
            elif self.n_iter_no_change is not None:
                # we changed some labels => reset patience
                patience = self.n_iter_no_change

        if self.n_iter_ == self.max_iter:
            self.termination_condition_ = "max_iter"
        if np.all(has_label):
            self.termination_condition_ = "all_labeled"
        if self.n_iter_ == self.max_iter and self.n_iter_no_change:
            warnings.warn("Maximum number of iterations reached before "
                          "early stopping. Consider increasing max_iter to "
                          "improve the fit.",
                          ConvergenceWarning)

        self.base_classifier_.fit(
            X[safe_mask(X, has_label)],
            self.y_labels_[safe_mask(self.y_labels_, has_label)])
        self.classes_ = self.base_classifier_.classes_
        return self

    @if_delegate_has_method(delegate='base_classifier')
    def predict(self, X):
        """Predict the classes of X.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples,)
            array with predicted labels
        """
        check_is_fitted(self, 'y_labeled_iter_')
        X = check_array(X)
        return self.base_classifier_.predict(X)

    def predict_proba(self, X):
        """Predict probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples, n_features)
            array with prediction probabilities
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.predict_proba(X)

    @if_delegate_has_method(delegate='base_classifier')
    def decision_function(self, X):
        """Calls decision function of the base_classifier

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples, n_features)
            result of the decision function of the base_classifier
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.decision_function(X)

    @if_delegate_has_method(delegate='base_classifier')
    def predict_log_proba(self, X):
        """Predict log probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples, n_features)
            array with log prediction probabilities
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.predict_log_proba(X)

    @if_delegate_has_method(delegate='base_classifier')
    def score(self, X, y):
        """Calls score on the base_classifier

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        y : array-like, shape=(n_samples,)
            array representing the labels

        Returns
        -------
        score : float
            result of calling score on the base_classifier
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.score(X, y)
