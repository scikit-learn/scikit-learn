from ..base import clone
from ..base import ClassifierMixin, RegressorMixin, TransformerMixin
from ..base import is_classifier, is_regressor
from ..base import MetaEstimatorMixin

from sklearn.externals.joblib import Parallel, delayed

from sklearn.linear_model import LogisticRegression

from ..model_selection import cross_val_predict
from ..model_selection._split import check_cv

from ..utils import check_random_state
from ..utils.metaestimators import _BaseComposition
from ..utils.validation import has_fit_parameter


class StackingClassifier(_BaseComposition, MetaEstimatorMixin, ClassifierMixin,
                         TransformerMixin):

    def __init__(self, estimators, classifier=None, cv=3,
                 method_estimators='auto', method_classifier='auto',
                 n_jobs=1, random_state=None):
        self.estimators = estimators
        self.classifier = classifier
        self.cv = cv
        self.method_estimators = method_estimators
        self.method_classifier = method_classifier
        self.n_jobs = n_jobs
        self.random_state = random_state

    @property
    def named_estimators(self):
        return dict(self.estimators)

    @staticmethod
    def _method_name(estimator, method):
        if method == 'auto':
            if getattr(estimator, 'predict_proba', None):
                return 'predict_proba'
            elif getattr(estimator, 'decision_function', None):
                return 'decision_function'
            else:
                return 'predict'
        else:
            return method

    def fit(self, X, y, sample_weight=None):

        random_state = check_random_state(self.random_state)

        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError('Invalid `estimators` attribute, `estimators`'
                                 ' should be a list of (string, estimator)'
                                 ' tuples')

        if self.classifier is None:
            self.classifier_ = LogisticRegression(random_state=random_state)
        else:
            if not is_classifier(self.classifier):
                raise AttributeError('`classifier` attribute should be a '
                                     'classifier.')

        if sample_weight is not None:
            for name, step in self.estimators:
                if not has_fit_parameter(step, 'sample_weight'):
                    raise ValueError('Underlying estimator \'%s\' does not'
                                     ' support sample weights.' % name)

        names, estimators_ = zip(*self.estimators)
        self._validate_names(names)

        method_estimators = ([self.method_estimators] * len(estimators_)
                             if self.method_estimators == 'auto'
                             else self.method_estimators)

        self._method_estimators = [
            self._method_name(est, meth)
            for est, meth in zip(estimators_, method_estimators)]

        cv = check_cv(self.cv)

        return self
