import numpy as np
from ..base import BaseEstimator, ClassifierMixin
from ..base import clone, is_classifier
from ..model_selection import cross_val_predict
from ..preprocessing import LabelEncoder
from ..externals.joblib import Parallel, delayed
from ..utils.metaestimators import if_delegate_has_method
from ..utils.validation import has_fit_parameter, check_is_fitted
from ..externals import six


def _fit_and_predict(estimator, method, cv, X, y, sample_weight):
    """Private function used to fit an estimator within a job."""
    if sample_weight is not None:
        estimator.fit(X, y, sample_weight)
        scores = cross_val_predict(
            estimator, X, y, method=method, cv=cv,
            fit_params={'sample_weight': sample_weight})
    else:
        estimator.fit(X, y)
        scores = cross_val_predict(
            estimator, X, y, method=method, cv=cv)
    return estimator, scores


class StackingClassifier(BaseEstimator, ClassifierMixin):
    """ Stacking classifier for combining unfitted estimators

    The cross-validated predictions of unfitted estimators will be used as
    inputs to a meta-estimator for making a final prediction

    Parameters
    ----------
    estimators : list of (string, estimator) tuples
        Invoking the ``fit`` method on the ``StackingClassifier`` will fit
        clones of those original estimators that will be stored in the class
        attribute `self.estimators_`.

    meta_estimator : estimator
        The meta-estimator to combine the predictions of each individual
        estimator

    method : string, optional, default='predict_proba'
        Specifies which method of the estimators will be called to generate
        inputs to the meta_estimator

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy to generate inputs
        to the meta_estimator
        Possible inputs for cv are:
          - None, to use the default 3-fold cross validation,
          - integer, to specify the number of folds in a `(Stratified)KFold`,
          - An object to be used as a cross-validation generator.
          - An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_jobs : int, optional (default=1)
        The number of jobs to run in parallel for ``fit``.
        If -1, then the number of jobs is set to the number of cores.

    Attributes
    ----------
    estimators_ : list of classifiers
        The collection of fitted sub-estimators.

    meta_estimators_ : classifier
        Fitted meta-estimator

    classes_ : array-like, shape = [n_predictions]
        The classes labels.
    """

    def __init__(self, estimators, meta_estimator, method='predict_proba',
                 cv=None, n_jobs=1):
        self.estimators = estimators
        self.named_estimators = dict(estimators)
        self.meta_estimator = meta_estimator
        self.cv = cv
        self.method = method
        self.n_jobs = n_jobs

    def fit(self, X, y, sample_weight=None):
        """ Fit the estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted.
            Note that this is supported only if all underlying estimators
            support sample weights.

        Returns
        -------
        self : object
        """
        if isinstance(y, np.ndarray) and len(y.shape) > 1 and y.shape[1] > 1:
            raise NotImplementedError('Multilabel and multi-output'
                                      ' classification is not supported.')

        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError('Invalid `estimators` attribute, `estimators`'
                                 ' should be a list of (string, estimator)'
                                 ' tuples')

        if not is_classifier(self.meta_estimator):
            raise AttributeError('Invalid `meta_estimator` attribute, '
                                 '`meta_estimator` should be a classifier')

        for name, step in self.estimators:
            if not hasattr(step, self.method):
                raise ValueError('Underlying estimator {0} does not support '
                                 '{1}.'.format(name, self.method))

        if sample_weight is not None:
            for name, step in self.estimators:
                if not has_fit_parameter(step, 'sample_weight'):
                    raise ValueError('Underlying estimator {0} does not'
                                     'support sample weights.'.format(name))

        self.le_ = LabelEncoder()
        self.le_.fit(y)
        self.classes_ = self.le_.classes_

        transformed_y = self.le_.transform(y)

        prediction_blocks = Parallel(n_jobs=self.n_jobs)(
                delayed(_fit_and_predict)(clone(clf), self.method, self.cv, X,
                                          transformed_y, sample_weight)
                for _, clf in self.estimators)

        scores = self._form_meta_inputs([p for _, p in prediction_blocks])
        self.estimators_ = [est for est, _ in prediction_blocks]
        self.meta_estimator_ = clone(self.meta_estimator)
        self.meta_estimator_.fit(scores, transformed_y,
                                 sample_weight=sample_weight)
        return self

    def _form_meta_inputs(self, predicted):
        """ Concatenate estimator predictions to form input to meta-estimator
        """
        if self.method == 'predict_log_proba':
            # Replace inf
            predicted = [np.nan_to_num(s) for s in predicted]

        if self.method == 'predict_proba':
            # Remove first column to avoid multicollinearity since sum of
            # probabilities equals 1
            predicted = [s[:, 1:] for s in predicted]

        return np.concatenate([s.reshape(-1, 1) if s.ndim == 1 else s
                               for s in predicted], axis=1)

    def _est_predict(self, X):
        """ Generate input to meta_estimator from predictions of estimators
        """
        return self._form_meta_inputs([getattr(clf, self.method)(X)
                                       for clf in self.estimators_])

    @if_delegate_has_method(delegate=('meta_estimator_', 'meta_estimator'))
    def predict(self, X):
        """Predict class labels for X.

        Only available if fitted and ``meta_estimator`` supports ``predict``.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Predicted class label per sample.
        """
        check_is_fitted(self, ['estimators_', 'meta_estimator_'])
        predicted = self.meta_estimator_.predict(self._est_predict(X))
        return self.le_.inverse_transform(predicted)

    @if_delegate_has_method(delegate=('meta_estimator_', 'meta_estimator'))
    def predict_proba(self, X):
        """ Return probability estimates for the test vector X.

        Only available if ``meta_estimator`` supports ``predict_proba``.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array-like, shape = [n_samples, n_classes]
            Returns the probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute `classes_`.
        """
        check_is_fitted(self, ['estimators_', 'meta_estimator_'])
        return self.meta_estimator_.predict_proba(self._est_predict(X))

    @if_delegate_has_method(delegate=('best_estimator_', 'estimator'))
    def predict_log_proba(self, X):
        """ Return log-probability estimates for the test vector X.

        Only available if ``meta_estimator`` supports ``predict_log_proba``.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array-like, shape = [n_samples, n_classes]
            Returns the log-probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute `classes_`.
        """
        check_is_fitted(self, ['estimators_', 'meta_estimator_'])
        return self.meta_estimator_.predict_log_proba(self._est_predict(X))

    @if_delegate_has_method(delegate=('best_estimator_', 'estimator'))
    def decision_function(self, X):
        """Predict confidence scores for samples.

        Only available if ``meta_estimator`` supports ``predict_log_proba``.
        The confidence score for a sample is the signed distance of that
        sample to the hyperplane.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.

        Returns
        -------
        array, shape=(n_samples,) if n_classes == 2 else (n_samples, n_classes)
            Confidence scores per (sample, class) combination. In the binary
            case, confidence score for self.classes_[1] where >0 means this
            class would be predicted.
        """
        check_is_fitted(self, ['estimators_', 'meta_estimator_'])
        return self.meta_estimator_.decision_function(self._est_predict(X))

    def get_params(self, deep=True):
        """Return estimator parameter names for GridSearch support
        """
        if not deep:
            return super(StackingClassifier, self).get_params(deep=False)
        else:
            out = super(StackingClassifier, self).get_params(deep=False)
            out.update(self.named_estimators.copy())
            out.update({'meta_estimator': self.meta_estimator}.copy())
            for name, step in six.iteritems(self.named_estimators):
                for key, value in six.iteritems(step.get_params(deep=True)):
                    out['{0}__{1}'.format(name, key)] = value
            for key, value in six.iteritems(
                    self.meta_estimator.get_params(deep=True)):
                out['{0}__{1}'.format('meta_estimator', key)] = value
            return out
