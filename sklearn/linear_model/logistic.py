# Authors: Fabian Pedregosa
#          Alexandre Gramfort
# License: 3-clause BSD

"""
Logistic Regression
"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Fabian Pedregosa <f@bianp.net>
#         Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>

import numbers
>>>>>>> Implementation of logistic_regression_path.
import numpy as np
from scipy import optimize, sparse

from .base import LinearClassifierMixin, SparseCoefMixin, BaseEstimator
from ..feature_selection.from_model import _LearntSelectorMixin
from ..preprocessing import LabelEncoder
from ..svm.base import BaseLibLinear
from ..utils import as_float_array
from ..externals.joblib import Parallel, delayed
from ..cross_validation import check_cv
from ..utils.optimize import newton_cg
from ..externals import six
from ..metrics import SCORERS


# .. some helper functions for logistic_regression_path ..
def _phi(t, copy=True):
    # helper function: return 1. / (1 + np.exp(-t))
    if copy:
        t = np.copy(t)
    t *= -1.
    t = np.exp(t, t)
    t += 1
    t = np.reciprocal(t, t)
    return t


def _logistic_loss_and_grad(w, X, y, alpha):
    # the logistic loss and its gradient
    z = X.dot(w)
    yz = y * z
    out = np.empty_like(yz)
    idx = yz > 0
    out[idx] = np.log(1 + np.exp(-yz[idx]))
    out[~idx] = (-yz[~idx] + np.log(1 + np.exp(yz[~idx])))
    out = out.sum() + .5 * alpha * w.dot(w)

    z = _phi(yz, copy=False)
    z0 = (z - 1) * y
    grad = X.T.dot(z0) + alpha * w
    return out, grad


def _logistic_loss(w, X, y, alpha):
    # the logistic loss and
    z = X.dot(w)
    yz = y * z
    out = np.empty_like(yz)
    idx = yz > 0
    out[idx] = np.log(1 + np.exp(-yz[idx]))
    out[~idx] = (-yz[~idx] + np.log(1 + np.exp(yz[~idx])))
    out = out.sum() + .5 * alpha * w.dot(w)
    #print 'Loss %r' % out
    return out


def _logistic_loss_grad_hess(w, X, y, alpha):
    # the logistic loss, its gradient, and the matvec application of the
    # Hessian
    z = X.dot(w)
    yz = y * z
    out = np.empty_like(yz)
    idx = yz > 0
    out[idx] = np.log(1 + np.exp(-yz[idx]))
    out[~idx] = (-yz[~idx] + np.log(1 + np.exp(yz[~idx])))
    out = out.sum() + .5 * alpha * w.dot(w)

    z = _phi(yz, copy=False)
    z0 = (z - 1) * y
    grad = X.T.dot(z0) + alpha * w

    # The mat-vec product of the Hessian
    d = z * (1 - z)
    d = np.sqrt(d, out=d)
    if sparse.issparse(X):
        dX = sparse.dia_matrix((d, 0), shape=(d.size, d.size)).dot(X)
    else:
        # Precompute as much as possible
        dX = d[:, np.newaxis] * X

    def Hs(s):
        ret = dX.T.dot(dX.dot(s))
        ret += alpha * s
        return ret
    #print 'Loss/grad/hess %r, %r' % (out, grad.dot(grad))
    return out, grad, Hs


def _logistic_loss_and_grad_intercept(w_c, X, y, alpha):
    w = w_c[:-1]
    c = w_c[-1]

    z = X.dot(w)
    z += c
    yz = y * z
    out = np.empty_like(yz)
    idx = yz > 0
    out[idx] = np.log(1 + np.exp(-yz[idx]))
    out[~idx] = (-yz[~idx] + np.log(1 + np.exp(yz[~idx])))
    out = out.sum() + .5 * alpha * w.dot(w)

    z = _phi(yz, copy=False)
    z0 = (z - 1) * y
    grad = np.empty_like(w_c)
    grad[:-1] = X.T.dot(z0) + alpha * w
    grad[-1] = z0.sum()
    return out, grad


def _logistic_loss_intercept(w_c, X, y, alpha):
    w = w_c[:-1]
    c = w_c[-1]

    z = X.dot(w)
    z += c
    yz = y * z
    out = np.empty_like(yz)
    idx = yz > 0
    out[idx] = np.log(1 + np.exp(-yz[idx]))
    out[~idx] = (-yz[~idx] + np.log(1 + np.exp(yz[~idx])))
    out = out.sum() + .5 * alpha * w.dot(w)

    #print 'Loss %r' % out
    return out


def _logistic_loss_grad_hess_intercept(w_c, X, y, alpha):
    w = w_c[:-1]
    c = w_c[-1]

    z = X.dot(w)
    z += c
    yz = y * z
    out = np.empty_like(yz)
    idx = yz > 0
    out[idx] = np.log(1 + np.exp(-yz[idx]))
    out[~idx] = (-yz[~idx] + np.log(1 + np.exp(yz[~idx])))
    out = out.sum() + .5 * alpha * w.dot(w)

    z = _phi(yz, copy=False)
    z0 = (z - 1) * y
    grad = np.empty_like(w_c)
    grad[:-1] = X.T.dot(z0) + alpha * w
    z0_sum = z0.sum()
    grad[-1] = z0_sum
    # The mat-vec product of the Hessian
    d = z * (1 - z)
    d = np.sqrt(d, out=d)
    if sparse.issparse(X):
        dX = sparse.dia_matrix((d, 0), shape=(d.size, d.size)).dot(X)
    else:
        # Precompute as much as possible
        dX = d[:, np.newaxis] * X
    def Hs(s):
        ret = np.empty_like(s)
        ret[:-1] = dX.T.dot(dX.dot(s[:-1]))
        ret[:-1] += alpha * s[:-1]
        # XXX: I am not sure that this last line of the Hessian is right
        # Without the intercept the Hessian is right, though
        ret[-1] = z0_sum * s[-1]
        return ret

    return out, grad, Hs

def logistic_regression_path(X, y, Cs=10, fit_intercept=True,
                             max_iter=100, gtol=1e-4, verbose=0,
                             solver='liblinear', callback=None,
                             coef=None):
    """
    Compute a Logistic Regression model for a list of regularization
    parameters.

    This is an implementation that uses the result of the previous model
    to speed up computations along the set of solutions, making it faster
    than sequentially calling LogisticRegression for the different parameters.

    Parameters
    ----------
    X : array-like or sparse matrix, shape (n_samples, n_features)
        Input data

    y : array-like, shape (n_samples,)
        Input data, target values

    Cs : array-like or integer of shape (n_cs,)
        List of values for the regularization parameter or integer specifying
        the number of regularization parameters that should be used. In this
        case, the parameters will be chosen in a logarithmic scale between
        1e-4 and 1e4.

    fit_intercept : boolean
        Whether to fit an intercept for the model. In this case the shape of
        the returned array is (n_cs, n_features + 1).

    max_iter : integer
        Maximum number of iterations for the solver.

    gtol : float
        Stopping criterion. The iteration will stop when
        ``max{|g_i | i = 1, ..., n} <= gtol``
        where ``g_i`` is the i-th component of the gradient. Only used
        by the methods 'lbfgs' and 'trust-ncg'

    verbose: int
        Print convergence message if True.

    solver : {'lbfgs', 'newton-cg', 'liblinear'}
        Numerical solver to use.

    callback : callable
        Function to be called before and after the fit of each regularization
        parameter. Must have the signature callback(w, X, y, alpha).

    coef: array-lime, shape (n_features,)
        Initialization value for coefficients of logistic regression.

    Returns
    -------
    coefs: array of shape (n_cs, n_features) or (n_cs, n_features + 1)
        List of coefficients for the Logistic Regression model. If
        fit_intercept is set to True then the seconds dimension will be
        n_features + 1, where the last item represents the intercept.


    Notes
    -----
    You might get slighly different results with the solver trust-ncg than
    with the others since this uses LIBLINEAR penalizes the intercept.
    """
    if isinstance(Cs, numbers.Integral):
        Cs = np.logspace(-4, 4, Cs)
    Cs = np.sort(Cs)
    y = np.sign(y - np.asarray(y).mean())
    X = as_float_array(X, copy=False)
    if not (np.unique(y).size == 2):
        raise NotImplementedError('logistic_regression_path is currently only '
                                  'implemented for the binary class case')
    if fit_intercept:
        w0 = np.zeros(X.shape[1] + 1)
        func = _logistic_loss_and_grad_intercept
    else:
        w0 = np.zeros(X.shape[1])
        func = _logistic_loss_and_grad

    if coef is not None:
        # it must work both giving the bias term and not
        if not coef.size in (X.shape[1], w0.size):
            raise ValueError('Initialization coef is not of correct shape')
        w0[:coef.size] = coef
    coefs = list()

    for C in Cs:
        if callback is not None:
            callback(w0, X, y, 1. / C)
        if solver == 'lbfgs':
            out = optimize.fmin_l_bfgs_b(
                func, w0, fprime=None,
                args=(X, y, 1. / C),
                iprint=verbose > 0, pgtol=gtol, maxiter=max_iter)
            w0 = out[0]
        elif solver == 'newton-cg':
            if fit_intercept:
                func_grad_hess = _logistic_loss_grad_hess_intercept
                func = _logistic_loss_intercept
            else:
                func_grad_hess = _logistic_loss_grad_hess
                func = _logistic_loss

            w0 = newton_cg(func_grad_hess, func, w0, args=(X, y, 1./C),
                        maxiter=max_iter)
        elif solver == 'liblinear':
            lr = LogisticRegression(C=C, fit_intercept=fit_intercept, tol=gtol)
            lr.fit(X, y)
            if fit_intercept:
                w0 = np.concatenate([lr.coef_.ravel(), lr.intercept_])
            else:
                w0 = lr.coef_.ravel()
        else:
            raise ValueError("solver must be one of {'liblinear', 'lbfgs', "
                             "'newton-cg'}, got '%s' instead" % solver)
        if callback is not None:
            callback(w0, X, y, 1. / C)
        coefs.append(w0)
    return coefs, Cs


# helper function for LogisticCV
def _log_reg_scoring_path(X, y, train, test, Cs=10, scoring=None,
                          fit_intercept=False,
                          max_iter=100, gtol=1e-4,
                          tol=1e-4, verbose=0, method='liblinear'):
    log_reg = LogisticRegression(fit_intercept=fit_intercept)
    log_reg._enc = LabelEncoder()
    log_reg._enc.fit_transform([-1, 1])

    coefs, Cs = logistic_regression_path(X[train], y[train], Cs=Cs,
                                         fit_intercept=fit_intercept,
                                         solver=method,
                                         max_iter=max_iter,
                                         gtol=gtol, verbose=verbose)
    scores = list()
    X_test = X[test]
    y_test = y[test]
    if isinstance(scoring, six.string_types):
        scoring = SCORERS[scoring]
    for w in coefs:
        if fit_intercept:
            log_reg.coef_ = w[np.newaxis, :-1]
            log_reg.intercept_ = w[-1]
        else:
            log_reg.coef_ = w[np.newaxis, :]
            log_reg.intercept_ = 0.
        if scoring is None:
            scores.append(log_reg.score(X_test, y_test))
        else:
            scores.append(scoring(log_reg, X_test, y_test))
    return coefs, Cs, np.array(scores)


class LogisticRegression(BaseLibLinear, LinearClassifierMixin,
                         _LearntSelectorMixin, SparseCoefMixin):
    """Logistic Regression (aka logit, MaxEnt) classifier.

    In the multiclass case, the training algorithm uses a one-vs.-all (OvA)
    scheme, rather than the "true" multinomial LR.

    This class implements L1 and L2 regularized logistic regression using the
    `liblinear` library. It can handle both dense and sparse input. Use
    C-ordered arrays or CSR matrices containing 64-bit floats for optimal
    performance; any other input format will be converted (and copied).

    Parameters
    ----------
    penalty : string, 'l1' or 'l2'
        Used to specify the norm used in the penalization.

    dual : boolean
        Dual or primal formulation. Dual formulation is only
        implemented for l2 penalty. Prefer dual=False when
        n_samples > n_features.

    C : float, optional (default=1.0)
        Inverse of regularization strength; must be a positive float.
        Like in support vector machines, smaller values specify stronger
        regularization.

    fit_intercept : bool, default: True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added the decision function.

    intercept_scaling : float, default: 1
        when self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased

    class_weight : {dict, 'auto'}, optional
        Over-/undersamples the samples of each class according to the given
        weights. If not given, all classes are supposed to have weight one.
        The 'auto' mode selects weights inversely proportional to class
        frequencies in the training set.

    random_state: int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    tol: float, optional
        Tolerance for stopping criteria.

    Attributes
    ----------
    `coef_` : array, shape = [n_classes, n_features]
        Coefficient of the features in the decision function.

    `intercept_` : array, shape = [n_classes]
        Intercept (a.k.a. bias) added to the decision function.
        If `fit_intercept` is set to False, the intercept is set to zero.

    See also
    --------
    SGDClassifier: incrementally trained logistic regression (when given
        the parameter ``loss="log"``).
    sklearn.svm.LinearSVC: learns SVM models using the same algorithm.

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    References:

    LIBLINEAR -- A Library for Large Linear Classification
        http://www.csie.ntu.edu.tw/~cjlin/liblinear/

    Hsiang-Fu Yu, Fang-Lan Huang, Chih-Jen Lin (2011). Dual coordinate descent
        methods for logistic regression and maximum entropy models.
        Machine Learning 85(1-2):41-75.
        http://www.csie.ntu.edu.tw/~cjlin/papers/maxent_dual.pdf
    """

    def __init__(self, penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None):

        super(LogisticRegression, self).__init__(
            penalty=penalty, dual=dual, loss='lr', tol=tol, C=C,
            fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
            class_weight=class_weight, random_state=random_state)

    def predict_proba(self, X):
        """Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in ``self.classes_``.
        """
        return self._predict_proba_lr(X)

    def predict_log_proba(self, X):
        """Log of probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the log-probability of the sample for each class in the
            model, where classes are ordered as they are in ``self.classes_``.
        """
        return np.log(self.predict_proba(X))


class LogisticRegressionCV(BaseEstimator, LinearClassifierMixin,
                           _LearntSelectorMixin):
    """Logistic Regression (aka logit, MaxEnt) classifier.

    This class implements L2 regularized logistic regression using and
    LBFGS optimizer.

    Parameters
    ----------
    Cs: list of floats, integer
        Each of the values in Cs describes the inverse of regularization
        strength and must be a positive float.
        Like in support vector machines, smaller values specify stronger
        regularization.

    fit_intercept: bool, default: True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added the decision function.

    max_iter: integer, optional
        Maximum number of iterations of the optimization algorithm.

    tol: float, optional
        Tolerance for stopping criteria.

    scoring: callabale
        Scoring function to use as cross-validation criteria.

    solver: {'newton-cg', 'lbfgs', 'liblinear'}
        Algorithm to use in the optimization problem.

    Attributes
    ----------
    `coef_` : array, shape = [n_classes-1, n_features]
        Coefficient of the features in the decision function.

        `coef_` is readonly property derived from `raw_coef_` that \
        follows the internal memory layout of liblinear.

    `intercept_` : array, shape = [n_classes-1]
        Intercept (a.k.a. bias) added to the decision function.
        It is available only when parameter intercept is set to True.

    See also
    --------
    LogisticRegression

    """

    def __init__(self, Cs=10, fit_intercept=True, cv=None, scoring=None,
                 solver='newton-cg', tol=1e-4, gtol=1e-4, max_iter=100,
                 n_jobs=1, verbose=False):
        self.Cs = Cs
        self.fit_intercept = fit_intercept
        self.cv = cv
        self.scoring = scoring
        self.tol = tol
        self.gtol = gtol
        self.max_iter = max_iter
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.solver = solver

    def fit(self, X, y):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target vector relative to X

        Returns
        -------
        self : object
            Returns self.
        """
        self._enc = LabelEncoder()
        X = as_float_array(X, copy=False)
        y = self._enc.fit_transform(y)
        if len(self.classes_) != 2:
            raise ValueError("LogisticRegressionCV works only on 2 "
                             "class problems. Please use "
                             "OneVsOneClassifier or OneVsRestClassifier")

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n"
                             "X has %s samples, but y has %s." %
                             (X.shape[0], y.shape[0]))

        # Transform to [-1, 1] classes, as y is [0, 1]
        y *= 2
        y -= 1

        # init cross-validation generator
        cv = check_cv(self.cv, X, y, classifier=True)
        folds = list(cv)

        fold_coefs_ = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                            delayed(_log_reg_scoring_path)(X, y, train, test,
                                        Cs=self.Cs,
                                        fit_intercept=self.fit_intercept,
                                        method=self.solver,
                                        max_iter=self.max_iter,
                                        gtol=self.gtol, tol=self.tol,
                                        verbose=max(0, self.verbose - 1),
                                        scoring=self.scoring,
                                    )
                              for train, test in folds
                            )
        coefs_paths, Cs, scores = zip(*fold_coefs_)
        self.Cs_ = Cs[0]
        self.coefs_paths_ = coefs_paths
        self.scores_ = np.array(scores)
        best_index = self.scores_.sum(axis=0).argmax()
        self.C_ = self.Cs_[best_index]
        coef_init = np.mean([c[best_index] for c in coefs_paths], axis=0)
        w = logistic_regression_path(
            X, y, Cs=[self.C_], fit_intercept=self.fit_intercept,
            coef=coef_init, solver=self.solver, max_iter=self.max_iter,
            gtol=self.gtol, verbose=max(0, self.verbose - 1))
        w = w[0][0][:, np.newaxis].T
        if self.fit_intercept:
            self.coef_ = w[:, :-1]
            self.intercept_ = w[:, -1]
        else:
            self.coef_ = w
            self.intercept_ = 0.
        return self

    @property
    def classes_(self):
        return self._enc.classes_

