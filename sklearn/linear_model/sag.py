import numpy as np
import scipy.sparse as sp

from abc import ABCMeta
import warnings

from .base import LinearClassifierMixin, LinearModel, SparseCoefMixin
from ..base import RegressorMixin, BaseEstimator
from ..utils import check_X_y, compute_class_weight, check_random_state
from ..utils import ConvergenceWarning
from ..utils.seq_dataset import ArrayDataset, CSRDataset
from ..externals import six
from ..externals.joblib import Parallel, delayed
from .sag_fast import Log, SquaredLoss
from .sag_fast import sag_sparse, get_auto_eta

MAX_INT = np.iinfo(np.int32).max

"""For sparse data intercept updates are scaled by this decay factor to avoid
intercept oscillation."""
SPARSE_INTERCEPT_DECAY = 0.01


# taken from http://stackoverflow.com/questions/1816958
# useful for passing instance methods to Parallel
def multiprocess_method(instance, name, args=()):
    "indirect caller for instance methods and multiprocessing"
    return getattr(instance, name)(*args)


# The inspiration for SAG comes from:
# "Minimizing Finite Sums with the Stochastic Average Gradient" by
# Mark Schmidt, Nicolas Le Roux, Francis Bach. 2013. <hal-00860051>
#
# https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf
class BaseSAG(six.with_metaclass(ABCMeta, SparseCoefMixin)):
    def __init__(self, alpha=0.0001, fit_intercept=True, max_iter=1000,
                 tol=0.001, verbose=0,
                 random_state=None, eta0='auto', warm_start=False):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        self.eta0 = eta0
        self.random_state = random_state
        self.warm_start = warm_start

        self._validate_params()

        self.coef_ = None
        self.intercept_ = None

        self.num_seen_ = None
        self.seen_ = None
        self.sum_gradient_ = None
        self.gradient_memory_ = None
        self.intercept_sum_gradient_ = None

    def _validate_params(self):
        if not isinstance(self.max_iter, int):
            raise ValueError("max_iter must be an integer")
        if self.max_iter < 1:
            raise ValueError("max_iter must be greater than 0")

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None, sum_gradient_init=None,
             gradient_memory_init=None, seen_init=None, num_seen_init=None,
             intercept_sum_gradient_init=None,
             weight_pos=1.0, weight_neg=1.0):

        n_samples, n_features = X.shape[0], X.shape[1]

        # initialize all parameters if there is no init
        if sample_weight is None:
            sample_weight = np.ones(n_samples, dtype=np.float64, order='C')

        if intercept_init is None:
            intercept_init = 0.0

        if intercept_sum_gradient_init is None:
            intercept_sum_gradient_init = 0.0

        if coef_init is None:
            coef_init = np.zeros(n_features, dtype=np.float64, order='C')

        if sum_gradient_init is None:
            sum_gradient_init = np.zeros(n_features, dtype=np.float64,
                                         order='C')

        if gradient_memory_init is None:
            gradient_memory_init = np.zeros(n_samples, dtype=np.float64,
                                            order='C')

        if seen_init is None:
            seen_init = np.zeros(n_samples, dtype=np.int32, order='C')

        if num_seen_init is None:
            num_seen_init = 0

        random_state = check_random_state(self.random_state)

        # check which type of Sequential Dataset is needed
        if sp.issparse(X):
            dataset = CSRDataset(X.data, X.indptr, X.indices,
                                 y, sample_weight,
                                 seed=random_state.randint(MAX_INT))
            intercept_decay = SPARSE_INTERCEPT_DECAY
        else:
            dataset = ArrayDataset(X, y, sample_weight,
                                   seed=random_state.randint(MAX_INT))
            intercept_decay = 1.0

        # set the eta0 if needed, 'auto' is 1 / 4L where L is the max sum of
        # squares for over all samples
        if self.eta0 == 'auto':
            step_size = get_auto_eta(dataset, self.alpha, n_samples,
                                     self.loss_function, self.fit_intercept)
        else:
            step_size = self.eta0

        intercept_, num_seen, max_iter_reached, intercept_sum_gradient = \
            sag_sparse(dataset, coef_init.ravel(),
                       intercept_init, n_samples,
                       n_features, self.tol,
                       self.max_iter,
                       self.loss_function,
                       step_size, self.alpha,
                       sum_gradient_init.ravel(),
                       gradient_memory_init.ravel(),
                       seen_init.ravel(),
                       num_seen_init, weight_pos,
                       weight_neg,
                       self.fit_intercept,
                       intercept_sum_gradient_init,
                       intercept_decay,
                       self.verbose)

        if max_iter_reached:
            warnings.warn("The max_iter was reached which means "
                          "the coef_ did not converge", ConvergenceWarning)

        return (coef_init.reshape(1, -1), intercept_,
                sum_gradient_init.reshape(1, -1),
                gradient_memory_init.reshape(1, -1),
                seen_init.reshape(1, -1),
                num_seen, intercept_sum_gradient)


class SAGClassifier(BaseSAG, LinearClassifierMixin, BaseEstimator):
    """Linear classifiers (SVM, logistic regression, a.o.) with SAG training.

    This estimator implements regularized linear models with stochastic
    average gradient (SAG) learning: the gradient of the loss is estimated
    using a random sample from the dataset. The weights are then updated
    according to the sum of gradients seen thus far divided by the number of
    unique samples seen. The inspiration for SAG comes from "Minimizing Finite
    Sums with the Stochastic Average Gradient" by Mark Schmidt, Nicolas Le
    Roux, and Francis Bach. 2013. <hal-00860051>
    https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf

    IMPORTANT NOTE: SAGClassifier and models from linear_model in general
    depend on columns that are on the same scale. You can make sure that the
    data will be normalized by using sklearn.preprocessing.StandardScaler on
    your data before passing it to the fit method.

    This implementation works with data represented as dense or sparse arrays
    of floating point values for the features. It will fit the data according
    to log loss.

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean norm
    L2.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the regularization term. Defaults to 0.0001

    fit_intercept: bool, optional
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    max_iter: int, optional
        The max number of passes over the training data if the stopping
        criterea is not reached. Defaults to 1000.

    tol: double, optional
        The stopping criterea for the weights. THe iterations will stop when
        max(change in weights) / max(weights) < tol. Defaults to .001

    random_state: int or numpy.random.RandomState, optional
        The random_state of the pseudo random number generator to use when
        sampling the data.

    verbose: integer, optional
        The verbosity level

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    eta0 : double or "auto"
        The initial learning rate. The default value is 0.001.

    class_weight : dict, {class_label : weight} or "auto" or None, optional
        Preset for the class_weight fit parameter.

        Weights associated with classes. If not given, all classes
        are supposed to have weight one.

        The "auto" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.


    Attributes
    ----------
    coef_ : array, shape (1, n_features) if n_classes == 2 else (n_classes,
    n_features)
        Weights assigned to the features.

    intercept_ : array, shape (1,) if n_classes == 2 else (n_classes,)
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> clf = linear_model.SAGClassifier()
    >>> clf.fit(X, Y)
    ... #doctest: +NORMALIZE_WHITESPACE
    SAGClassifier(alpha=0.0001, class_weight=None,
                  eta0='auto', fit_intercept=True,
                  max_iter=1000, n_jobs=1, random_state=None,
                  tol=0.001, verbose=0, warm_start=False)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    SGDClassifier, LinearSVC, LogisticRegression, Perceptron

    """
    def __init__(self, alpha=0.0001,
                 fit_intercept=True, max_iter=1000, tol=0.001, verbose=0,
                 n_jobs=1, random_state=None,
                 eta0='auto', class_weight=None, warm_start=False):
        self.n_jobs = n_jobs
        self.class_weight = class_weight
        self.loss_function = Log()
        super(SAGClassifier, self).__init__(alpha=alpha,
                                            fit_intercept=fit_intercept,
                                            max_iter=max_iter,
                                            verbose=verbose,
                                            random_state=random_state,
                                            tol=tol,
                                            eta0=eta0,
                                            warm_start=warm_start)

    """Fit linear model with Stochastic Average Gradient.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : numpy array, shape (n_samples,)
            Target values

        sample_weight : array-like, shape (n_samples,), optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : returns an instance of self.
        """
    def fit(self, X, y, sample_weight=None):
        X, y = check_X_y(X, y, "csr", copy=False, order='C',
                         dtype=np.float64)
        n_samples, n_features = X.shape[0], X.shape[1]

        self.classes_ = np.unique(y)
        self.expanded_class_weight_ = compute_class_weight(self.class_weight,
                                                           self.classes_, y)

        if self.classes_.shape[0] <= 1:
            # there is only one class
            raise ValueError("The number of class labels must be "
                             "greater than one.")
        elif self.classes_.shape[0] == 2:
            # binary classifier
            (coef, intercept, sum_gradient, gradient_memory,
             seen, num_seen, intercept_sum_gradient) = \
                self._fit_target_class(X, y, self.classes_[1], sample_weight)
        else:
            # multiclass classifier
            coef = []
            intercept = []
            sum_gradient = []
            gradient_memory = []
            seen = []
            num_seen = []
            intercept_sum_gradient = []

            # perform a fit for all classes, one verse all
            results = Parallel(n_jobs=self.n_jobs,
                               backend="threading",
                               verbose=self.verbose)(
                # we have to use a call to multiprocess_method instead of the
                # plain instance method because pickle will not work on
                # instance methods in python 2.6 and 2.7
                delayed(multiprocess_method)(self, "_fit_target_class",
                                             (X, y, cl, sample_weight))
                for cl in self.classes_)

            # append results to the correct array
            for (coef_cl, intercept_cl, sum_gradient_cl, gradient_memory_cl,
                 seen_cl, num_seen_cl, intercept_sum_gradient_cl) in results:
                coef.append(coef_cl)
                intercept.append(intercept_cl)
                sum_gradient.append(sum_gradient_cl)
                gradient_memory.append(gradient_memory_cl)
                seen.append(seen_cl)
                num_seen.append(num_seen_cl)
                intercept_sum_gradient.append(intercept_sum_gradient_cl)

            # stack all arrays to transform into np arrays
            coef = np.vstack(coef)
            intercept = np.array(intercept)
            sum_gradient = np.vstack(sum_gradient)
            gradient_memory = np.vstack(gradient_memory)
            seen = np.vstack(seen)
            num_seen = np.array(num_seen)
            intercept_sum_gradient = np.array(intercept_sum_gradient)

        self.coef_ = coef
        self.intercept_ = intercept
        self.sum_gradient_ = sum_gradient
        self.gradient_memory_ = gradient_memory
        self.seen_ = seen
        self.num_seen_ = num_seen
        self.intercept_sum_gradient_ = intercept_sum_gradient

        return self

    def _fit_target_class(self, X, y, target_class, sample_weight=None):
        coef_init = None
        intercept_init = None
        sum_gradient_init = None
        gradient_memory_init = None
        seen_init = None
        num_seen_init = None
        intercept_sum_gradient_init = None

        if self.classes_.shape[0] == 2:
            if self.warm_start:
                # init parameters for binary classifier
                coef_init = self.coef_
                intercept_init = self.intercept_
                sum_gradient_init = self.sum_gradient_
                gradient_memory_init = self.gradient_memory_
                seen_init = self.seen_
                num_seen_init = self.num_seen_
                intercept_sum_gradient_init = \
                    self.intercept_sum_gradient_

            weight_pos = self.expanded_class_weight_[1]
            weight_neg = self.expanded_class_weight_[0]
        else:
            class_index = np.where(self.classes_ == target_class)[0][0]
            if self.warm_start:
                # init parameters for multi-class classifier
                if self.coef_ is not None:
                    coef_init = self.coef_[class_index]
                if self.intercept_ is not None:
                    intercept_init = self.intercept_[class_index]
                if self.sum_gradient_ is not None:
                    sum_gradient_init = self.sum_gradient_[class_index]
                if self.gradient_memory_ is not None:
                    gradient_memory_init = self.gradient_memory_[class_index]
                if self.seen_ is not None:
                    seen_init = self.seen_[class_index]
                if self.num_seen_ is not None:
                    num_seen_init = self.num_seen_[class_index]
                if self.intercept_sum_gradient_ is not None:
                    intercept_sum_gradient_init = \
                        self.intercept_sum_gradient_[class_index]

            weight_pos = self.expanded_class_weight_[class_index]
            weight_neg = 1.0

        n_samples, n_features = X.shape[0], X.shape[1]

        y_encoded = np.ones(n_samples)
        y_encoded[y != target_class] = -1.0

        return super(SAGClassifier, self).\
            _fit(X, y_encoded,
                 coef_init, intercept_init,
                 sample_weight,
                 sum_gradient_init,
                 gradient_memory_init,
                 seen_init, num_seen_init,
                 intercept_sum_gradient_init,
                 weight_pos, weight_neg)


class SAGRegressor(BaseSAG, LinearModel, RegressorMixin,
                   BaseEstimator):
    """Linear model fitted by minimizing a regularized empirical loss with SAG

    SAG stands for Stochastic Average Gradient: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a constant learning rate. The inspiration for SAG comes from "Minimizing
    Finite Sums with the Stochastic Average Gradient" by Mark Schmidt,
    Nicolas Le Roux, and Francis Bach. 2013. <hal-00860051>
    https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf

    IMPORTANT NOTE: SAGRegressor and models from linear_model in general depend
    on columns that are on the same scale. You can make sure that the data will
    be normalized by using sklearn.preprocessing.StandardScaler on your data
    before passing it to the fit method.

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using the squared euclidean norm
    L2.

    This implementation works with data represented as dense or sparse numpy
    arrays of floating point values for the features.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the regularization term. Defaults to 0.0001

    fit_intercept: bool, optional
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    max_iter: int, optional
        The max number of passes over the training data if the stopping
        criterea is not reached. Defaults to 1000.

    tol: double, optional
        The stopping criterea for the weights. THe iterations will stop when
        max(change in weights) / max(weights) < tol. Defaults to .001

    random_state: int or numpy.random.RandomState, optional
        The random_state of the pseudo random number generator to use when
        sampling the data.

    verbose: integer, optional
        The verbosity level.

    eta0 : double or "auto"
        The initial learning rate [default 0.01].

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Weights asigned to the features.

    intercept_ : array, shape (1,)
        The intercept term.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = linear_model.SAGRegressor()
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    SAGRegressor(alpha=0.0001, eta0='auto',
                 fit_intercept=True, max_iter=1000, random_state=None,
                 tol=0.001, verbose=0, warm_start=False)

    See also
    --------
    SGDRegressor, Ridge, ElasticNet, Lasso, SVR

    """
    def __init__(self, alpha=0.0001, fit_intercept=True, max_iter=1000,
                 tol=0.001, verbose=0, random_state=None, eta0='auto',
                 warm_start=False):

        self.loss_function = SquaredLoss()
        super(SAGRegressor, self).__init__(alpha=alpha,
                                           fit_intercept=fit_intercept,
                                           max_iter=max_iter,
                                           verbose=verbose,
                                           random_state=random_state,
                                           tol=tol,
                                           eta0=eta0,
                                           warm_start=warm_start)

    """Fit linear model with Stochastic Average Gradient.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : numpy array, shape (n_samples,)
            Target values

        sample_weight : array-like, shape (n_samples,), optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : returns an instance of self.
        """
    def fit(self, X, y, sample_weight=None):
        X, y = check_X_y(X, y, "csr", copy=False, order='C', dtype=np.float64)
        y = y.astype(np.float64)

        coef_init = None
        intercept_init = None
        sum_gradient_init = None
        gradient_memory_init = None
        seen_init = None
        num_seen_init = None
        intercept_sum_gradient_init = None

        if self.warm_start:
            coef_init = self.coef_
            intercept_init = self.intercept_
            sum_gradient_init = self.sum_gradient_
            gradient_memory_init = self.gradient_memory_
            seen_init = self.seen_
            num_seen_init = self.num_seen_
            intercept_sum_gradient_init = self.intercept_sum_gradient_

        (self.coef_, self.intercept_, self.sum_gradient_,
         self.gradient_memory_, self.seen_, self.num_seen_,
         self.intercept_sum_gradient_) = \
            super(SAGRegressor, self)._fit(X, y, coef_init,
                                           intercept_init,
                                           sample_weight,
                                           sum_gradient_init,
                                           gradient_memory_init,
                                           seen_init, num_seen_init,
                                           intercept_sum_gradient_init)

        return self
