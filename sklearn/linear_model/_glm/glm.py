"""
Generalized Linear Models with Exponential Dispersion Family
"""

# Author: Christian Lorentzen <lorentzen.ch@googlemail.com>
# some parts and tricks stolen from other sklearn files.
# License: BSD 3 clause

import numbers

import numpy as np
import scipy.optimize

from ...base import BaseEstimator, RegressorMixin
from ...utils import check_array, check_X_y
from ...utils.optimize import _check_optimize_result
from ...utils.validation import check_is_fitted, _check_sample_weight
from .distribution import (
        ExponentialDispersionModel,
        TweedieDistribution,
        EDM_DISTRIBUTIONS
)
from .link import (
        Link,
        IdentityLink,
        LogLink,
)


class GeneralizedLinearRegressor(BaseEstimator, RegressorMixin):
    """Regression via a Generalized Linear Model (GLM) with penalties.

    GLMs based on a reproductive Exponential Dispersion Model (EDM) aim at
    fitting and predicting the mean of the target y as mu=h(X*w). Therefore,
    the fit minimizes the following objective function with L2
    priors as regularizer::

            1/(2*sum(s)) * deviance(y, h(X*w); s)
            + 1/2 * alpha * |w|_2

    with inverse link function h and s=sample_weight.
    The parameter ``alpha`` corresponds to the lambda parameter in glmnet.

    Read more in the :ref:`User Guide <Generalized_linear_regression>`.

    Parameters
    ----------
    alpha : float, optional (default=1)
        Constant that multiplies the penalty terms and thus determines the
        regularization strength.
        See the notes for the exact mathematical meaning of this
        parameter.``alpha = 0`` is equivalent to unpenalized GLMs. In this
        case, the design matrix X must have full column rank
        (no collinearities).

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    family : {'normal', 'poisson', 'gamma', 'inverse-gaussian'} \
            or an instance of class ExponentialDispersionModel, \
            optional(default='normal')
        The distributional assumption of the GLM, i.e. which distribution from
        the EDM, specifies the loss function to be minimized.

    link : {'auto', 'identity', 'log'} or an instance of class Link, \
            optional (default='auto')
        The link function of the GLM, i.e. mapping from linear predictor
        (X*coef) to expectation (mu). Option 'auto' sets the link depending on
        the chosen family as follows:

        - 'identity' for family 'normal'

        - 'log' for families 'poisson', 'gamma', 'inverse-gaussian'

    solver : 'lbfgs', optional (default='lbfgs')
        Algorithm to use in the optimization problem:

        'lbfgs'
            Calls scipy's L-BFGS-B optimizer.

    max_iter : int, optional (default=100)
        The maximal number of iterations for solver algorithms.

    tol : float, optional (default=1e-4)
        Stopping criterion. For the lbfgs solver,
        the iteration will stop when ``max{|g_i|, i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative) of
        the objective function.

    warm_start : boolean, optional (default=False)
        If set to ``True``, reuse the solution of the previous call to ``fit``
        as initialization for ``coef_`` and ``intercept_``.

    copy_X : boolean, optional, (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    check_input : boolean, optional (default=True)
        Allow to bypass several checks on input: y values in range of family,
        sample_weight non-negative.
        Don't use this parameter unless you know what you do.

    verbose : int, optional (default=0)
        For the lbfgs solver set verbose to any positive number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the linear predictor (X*coef_+intercept_) in
        the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    n_iter_ : int
        Actual number of iterations used in solver.

    Notes
    -----
    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments to be :math:`E[Y_i]=\\mu_i=h((Xw)_i)` and
    :math:`Var[Y_i]=\\frac{\\phi}{s_i} v(\\mu_i)`. The unit variance function
    :math:`v(\\mu_i)` is a property of and given by the specific EDM, see
    :ref:`User Guide <Generalized_linear_regression>`.

    The parameters :math:`w` (`coef_` and `intercept_`) are estimated by
    minimizing the deviance plus penalty term, which is equivalent to
    (penalized) maximum likelihood estimation.

    For alpha > 0, the feature matrix X should be standardized in order to
    penalize features equally strong. Call
    :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``.

    If the target y is a ratio, appropriate sample weights s should be
    provided.
    As an example, consider Poisson distributed counts z (integers) and
    weights s=exposure (time, money, persons years, ...). Then you fit
    y = z/s, i.e. ``GeneralizedLinearModel(family='poisson').fit(X, y,
    sample_weight=s)``. The weights are necessary for the right (finite
    sample) mean.
    Consider :math:`\\bar{y} = \\frac{\\sum_i s_i y_i}{\\sum_i s_i}`,
    in this case one might say that y has a 'scaled' Poisson distributions.
    The same holds for other distributions.

    References
    ----------
    .. McCullagh, Peter; Nelder, John (1989). Generalized Linear Models,
       Second Edition. Boca Raton: Chapman and Hall/CRC. ISBN 0-412-31760-5.

    .. Jørgensen, B. (1992). The theory of exponential dispersion models
       and analysis of deviance. Monografias de matemática, no. 51.  See also
       `Exponential dispersion model.
       <https://en.wikipedia.org/wiki/Exponential_dispersion_model>`_
    """
    def __init__(self, alpha=1.0,
                 fit_intercept=True, family='normal', link='auto',
                 solver='lbfgs', max_iter=100, tol=1e-4, warm_start=False,
                 copy_X=True, check_input=True, verbose=0):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.family = family
        self.link = link
        self.solver = solver
        self.max_iter = max_iter
        self.tol = tol
        self.warm_start = warm_start
        self.copy_X = copy_X
        self.check_input = check_input
        self.verbose = verbose

    def fit(self, X, y, sample_weight=None):
        """Fit a Generalized Linear Model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : {None, array-like}, shape (n_samples,),\
                optional (default=None)
            Individual weights w_i for each sample. Note that for an
            Exponential Dispersion Model (EDM), one has
            Var[Y_i]=phi/w_i * v(mu).
            If Y_i ~ EDM(mu, phi/w_i), then
            sum(w*Y)/sum(w) ~ EDM(mu, phi/sum(w)), i.e. the mean of y is a
            weighted average with weights=sample_weight.

        Returns
        -------
        self : returns an instance of self.
        """
        if isinstance(self.family, ExponentialDispersionModel):
            self._family_instance = self.family
        elif self.family in EDM_DISTRIBUTIONS:
            self._family_instance = EDM_DISTRIBUTIONS[self.family]()
        else:
            raise ValueError(
                "The family must be an instance of class"
                " ExponentialDispersionModel or an element of"
                " ['normal', 'poisson', 'gamma', 'inverse-gaussian']"
                "; got (family={0})".format(self.family))

        # Guarantee that self._link_instance is set to an instance of
        # class Link
        if isinstance(self.link, Link):
            self._link_instance = self.link
        else:
            if self.link == 'auto':
                if isinstance(self._family_instance, TweedieDistribution):
                    if self._family_instance.power <= 0:
                        self._link_instance = IdentityLink()
                    if self._family_instance.power >= 1:
                        self._link_instance = LogLink()
                else:
                    raise ValueError("No default link known for the "
                                     "specified distribution family. Please "
                                     "set link manually, i.e. not to 'auto'; "
                                     "got (link='auto', family={}"
                                     .format(self.family))
            elif self.link == 'identity':
                self._link_instance = IdentityLink()
            elif self.link == 'log':
                self._link_instance = LogLink()
            else:
                raise ValueError(
                    "The link must be an instance of class Link or "
                    "an element of ['auto', 'identity', 'log']; "
                    "got (link={0})".format(self.link))

        if not isinstance(self.alpha, numbers.Number) or self.alpha < 0:
            raise ValueError("Penalty term must be a non-negative number;"
                             " got (alpha={0})".format(self.alpha))
        if not isinstance(self.fit_intercept, bool):
            raise ValueError("The argument fit_intercept must be bool;"
                             " got {0}".format(self.fit_intercept))
        if self.solver not in ['lbfgs']:
            raise ValueError("GeneralizedLinearRegressor supports only solvers"
                             "'lbfgs'; got {0}".format(self.solver))
        solver = self.solver
        if (not isinstance(self.max_iter, numbers.Integral)
                or self.max_iter <= 0):
            raise ValueError("Maximum number of iteration must be a positive "
                             "integer;"
                             " got (max_iter={0!r})".format(self.max_iter))
        if not isinstance(self.tol, numbers.Number) or self.tol <= 0:
            raise ValueError("Tolerance for stopping criteria must be "
                             "positive; got (tol={0!r})".format(self.tol))
        if not isinstance(self.warm_start, bool):
            raise ValueError("The argument warm_start must be bool;"
                             " got {0}".format(self.warm_start))
        if not isinstance(self.copy_X, bool):
            raise ValueError("The argument copy_X must be bool;"
                             " got {0}".format(self.copy_X))
        if not isinstance(self.check_input, bool):
            raise ValueError("The argument check_input must be bool; got "
                             "(check_input={0})".format(self.check_input))

        family = self._family_instance
        link = self._link_instance

        _dtype = [np.float64, np.float32]
        _stype = ['csc', 'csr']
        X, y = check_X_y(X, y, accept_sparse=_stype,
                         dtype=_dtype, y_numeric=True, multi_output=False,
                         copy=self.copy_X)
        y = np.asarray(y, dtype=np.float64)

        weights = _check_sample_weight(sample_weight, X)

        n_samples, n_features = X.shape

        if self.check_input:
            if not np.all(family.in_y_range(y)):
                raise ValueError("Some value(s) of y are out of the valid "
                                 "range for family {0}"
                                 .format(family.__class__.__name__))
            # TODO: if alpha=0 check that X is not rank deficient

        # rescaling of sample_weight
        #
        # IMPORTANT NOTE: Since we want to minimize
        # 1/(2*sum(sample_weight)) * deviance + L2,
        # deviance = sum(sample_weight * unit_deviance),
        # we rescale weights such that sum(weights) = 1 and this becomes
        # 1/2*deviance + L2 with deviance=sum(weights * unit_deviance)
        weights_sum = np.sum(weights)
        weights = weights/weights_sum

        # initialization of coef = (intercept_, coef)
        # Note: The dispersion parameter phi does not enter the estimation
        #       of mu_i=E[y_i].

        if self.warm_start and hasattr(self, 'coef_'):
            if self.fit_intercept:
                coef = np.concatenate((np.array([self.intercept_]),
                                       self.coef_))
            else:
                coef = self.coef_
        else:
            if self.fit_intercept:
                coef = np.zeros(n_features+1)
                coef[0] = link.link(np.average(y, weights=weights))
            else:
                coef = np.zeros(n_features)

        # algorithms for optimization

        if solver == 'lbfgs':
            def func(coef, X, y, weights, alpha, family, link):
                mu, devp = \
                    family._mu_deviance_derivative(coef, X, y, weights, link)
                dev = family.deviance(y, mu, weights)
                intercept = (coef.size == X.shape[1] + 1)
                idx = 1 if intercept else 0  # offset if coef[0] is intercept
                L2 = alpha * coef[idx:]
                obj = 0.5 * dev + 0.5 * (coef[idx:] @ L2)
                objp = 0.5 * devp
                objp[idx:] += L2
                return obj, objp

            args = (X, y, weights, self.alpha, family, link)

            opt_res = scipy.optimize.minimize(
                func, coef, method="L-BFGS-B", jac=True,
                options={
                    "maxiter": self.max_iter,
                    "iprint": (self.verbose > 0) - 1,
                    "gtol": self.tol,
                    "ftol": 1e3*np.finfo(float).eps,
                },
                args=args)
            self.n_iter_ = _check_optimize_result("lbfgs", opt_res)
            coef = opt_res.x

        if self.fit_intercept:
            self.intercept_ = coef[0]
            self.coef_ = coef[1:]
        else:
            # set intercept to zero as the other linear models do
            self.intercept_ = 0.
            self.coef_ = coef

        return self

    def _linear_predictor(self, X):
        """Compute the linear_predictor = X*coef_ + intercept_.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        C : array, shape (n_samples,)
            Returns predicted values of linear predictor.
        """
        check_is_fitted(self, "coef_")
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'],
                        dtype='numeric', ensure_2d=True,
                        allow_nd=False)
        return X @ self.coef_ + self.intercept_

    def predict(self, X):
        """Predict using GLM with feature matrix X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y_pred : array, shape (n_samples,)
            Returns predicted values.
        """
        # check_array is done in _linear_predictor
        eta = self._linear_predictor(X)
        mu = self._link_instance.inverse(eta)
        return mu

    def score(self, X, y, sample_weight=None):
        """Compute D^2, the percentage of deviance explained.

        D^2 is a generalization of the coefficient of determination R^2.
        R^2 uses squared error and D^2 deviance. Note that those two are equal
        for ``family='normal'``.

        D^2 is defined as
        :math:`D^2 = 1-\\frac{D(y_{true},y_{pred})}{D_{null}}`,
        :math:`D_{null}` is the null deviance, i.e. the deviance of a model
        with intercept alone, which corresponds to :math:`y_{pred} = \\bar{y}`.
        The mean :math:`\\bar{y}` is averaged by sample_weight.
        Best possible score is 1.0 and it can be negative (because the model
        can be arbitrarily worse).

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Test samples.

        y : array-like, shape (n_samples,)
            True values of target.

        sample_weight : {None, array-like}, shape (n_samples,), optional \
                (default=None)
            Sample weights.

        Returns
        -------
        score : float
            D^2 of self.predict(X) w.r.t. y.
        """
        # Note, default score defined in RegressorMixin is R^2 score.
        # TODO: make D^2 a score function in module metrics (and thereby get
        #       input validation and so on)
        weights = _check_sample_weight(sample_weight, X)
        mu = self.predict(X)
        dev = self._family_instance.deviance(y, mu, weights=weights)
        y_mean = np.average(y, weights=weights)
        dev_null = self._family_instance.deviance(y, y_mean, weights=weights)
        return 1. - dev / dev_null

    def _more_tags(self):
        return {"requires_positive_y": True}


class PoissonRegressor(GeneralizedLinearRegressor):
    """Regression with the response variable y following a Poisson distribution

    GLMs based on a reproductive Exponential Dispersion Model (EDM) aim at
    fitting and predicting the mean of the target y as mu=h(X*w).
    The fit minimizes the following objective function with L2 regularization::

            1/(2*sum(s)) * deviance(y, h(X*w); s) + 1/2 * alpha * ||w||_2^2

    with inverse link function h and s=sample_weight. Note that for
    ``sample_weight=None``, one has s_i=1 and sum(s)=n_samples).

    Read more in the :ref:`User Guide <Generalized_linear_regression>`.

    Parameters
    ----------
    alpha : float, optional (default=1)
        Constant that multiplies the penalty terms and thus determines the
        regularization strength.
        See the notes for the exact mathematical meaning of this
        parameter.``alpha = 0`` is equivalent to unpenalized GLMs. In this
        case, the design matrix X must have full column rank
        (no collinearities).

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    solver : {'lbfgs'}, optional (default='lbfgs')
        Algorithm to use in the optimization problem:

        'lbfgs'
            Calls scipy's L-BFGS-B optimizer.

    max_iter : int, optional (default=100)
        The maximal number of iterations for solver algorithms.

    tol : float, optional (default=1e-4)
        Stopping criterion. For the lbfgs solver,
        the iteration will stop when ``max{|g_i|, i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative) of
        the objective function.

    warm_start : boolean, optional (default=False)
        If set to ``True``, reuse the solution of the previous call to ``fit``
        as initialization for ``coef_`` and ``intercept_`` .

    copy_X : boolean, optional, (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    verbose : int, optional (default=0)
        For the lbfgs solver set verbose to any positive number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the linear predictor (X*coef_+intercept_) in
        the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    n_iter_ : int
        Actual number of iterations used in solver.

    Notes
    -----
    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments to be :math:`E[Y_i]=\\mu_i=h((Xw)_i)` and
    :math:`Var[Y_i]=\\frac{\\phi}{s_i} v(\\mu_i)`. The unit variance function
    :math:`v(\\mu_i)` is a property of and given by the specific EDM, see
    :ref:`User Guide <Generalized_linear_regression>`.

    The parameters :math:`w` (`coef_` and `intercept_`) are estimated by
    minimizing the deviance plus penalty term, which is equivalent to
    (penalized) maximum likelihood estimation.

    For alpha > 0, the feature matrix X should be standardized in order to
    penalize features equally strong.

    If the target y is a ratio, appropriate sample weights s should be
    provided.
    As an example, consider Poisson distributed counts z (integers) and
    weights s=exposure (time, money, persons years, ...). Then you fit
    y = z/s, i.e. ``PoissonRegressor().fit(X, y, sample_weight=s)``.
    The weights are necessary for the right (finite sample) mean.
    Consider :math:`\\bar{y} = \\frac{\\sum_i s_i y_i}{\\sum_i s_i}`,
    in this case one might say that y has a 'scaled' Poisson distributions.

    References
    ----------
    .. McCullagh, Peter; Nelder, John (1989). Generalized Linear Models,
       Second Edition. Boca Raton: Chapman and Hall/CRC. ISBN 0-412-31760-5.

    .. Jørgensen, B. (1992). The theory of exponential dispersion models
       and analysis of deviance. Monografias de matemática, no. 51.  See also
       `Exponential dispersion model.
       <https://en.wikipedia.org/wiki/Exponential_dispersion_model>`_
    """
    def __init__(self, alpha=1.0, fit_intercept=True, link='log',
                 solver='lbfgs', max_iter=100, tol=1e-4, warm_start=False,
                 copy_X=True, check_input=True, verbose=0):

        super().__init__(alpha=alpha, fit_intercept=fit_intercept,
                         family="poisson", link=link,
                         solver=solver, max_iter=max_iter, tol=tol,
                         warm_start=warm_start, copy_X=copy_X, verbose=verbose)

    @property
    def family(self):
        return "poisson"

    @family.setter
    def family(self, value):
        if value != "poisson":
            raise ValueError("PoissonRegressor.family must be 'poisson'!")


class GammaRegressor(GeneralizedLinearRegressor):
    """Regression with the response variable y following a Gamma distribution

    GLMs based on a reproductive Exponential Dispersion Model (EDM) aim at
    fitting and predicting the mean of the target y as mu=h(X*w).
    The fit minimizes the following objective function with L2 regularization::

            1/(2*sum(s)) * deviance(y, h(X*w); s) + 1/2 * alpha * ||w||_2^2

    with inverse link function h and s=sample_weight. Note that for
    ``sample_weight=None``, one has s_i=1 and sum(s)=n_samples).

    Read more in the :ref:`User Guide <Generalized_linear_regression>`.

    Parameters
    ----------
    alpha : float, optional (default=1)
        Constant that multiplies the penalty terms and thus determines the
        regularization strength.
        See the notes for the exact mathematical meaning of this
        parameter.``alpha = 0`` is equivalent to unpenalized GLMs. In this
        case, the design matrix X must have full column rank
        (no collinearities).

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    solver : {'lbfgs'}, optional (default='lbfgs')
        Algorithm to use in the optimization problem:

        'lbfgs'
            Calls scipy's L-BFGS-B optimizer.

    max_iter : int, optional (default=100)
        The maximal number of iterations for solver algorithms.

    tol : float, optional (default=1e-4)
        Stopping criterion. For the lbfgs solver,
        the iteration will stop when ``max{|g_i|, i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative) of
        the objective function.

    warm_start : boolean, optional (default=False)
        If set to ``True``, reuse the solution of the previous call to ``fit``
        as initialization for ``coef_`` and ``intercept_`` .

    copy_X : boolean, optional, (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    verbose : int, optional (default=0)
        For the lbfgs solver set verbose to any positive number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the linear predictor (X*coef_+intercept_) in
        the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    n_iter_ : int
        Actual number of iterations used in solver.

    Notes
    -----
    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments to be :math:`E[Y_i]=\\mu_i=h((Xw)_i)` and
    :math:`Var[Y_i]=\\frac{\\phi}{s_i} v(\\mu_i)`. The unit variance function
    :math:`v(\\mu_i)` is a property of and given by the specific EDM, see
    :ref:`User Guide <Generalized_linear_regression>`.

    The parameters :math:`w` (`coef_` and `intercept_`) are estimated by
    minimizing the deviance plus penalty term, which is equivalent to
    (penalized) maximum likelihood estimation.


    References
    ----------
    .. McCullagh, Peter; Nelder, John (1989). Generalized Linear Models,
       Second Edition. Boca Raton: Chapman and Hall/CRC. ISBN 0-412-31760-5.

    .. Jørgensen, B. (1992). The theory of exponential dispersion models
       and analysis of deviance. Monografias de matemática, no. 51.  See also
       `Exponential dispersion model.
       <https://en.wikipedia.org/wiki/Exponential_dispersion_model>`_
    """
    def __init__(self, alpha=1.0, fit_intercept=True, link='log',
                 solver='lbfgs', max_iter=100, tol=1e-4, warm_start=False,
                 copy_X=True, check_input=True, verbose=0):

        super().__init__(alpha=alpha, fit_intercept=fit_intercept,
                         family="gamma", link=link,
                         solver=solver, max_iter=max_iter, tol=tol,
                         warm_start=warm_start, copy_X=copy_X, verbose=verbose)

    @property
    def family(self):
        return "gamma"

    @family.setter
    def family(self, value):
        if value != "gamma":
            raise ValueError("GammaRegressor.family must be 'gamma'!")


class TweedieRegressor(GeneralizedLinearRegressor):
    """Regression with the response variable y following a Tweedie distribution

    GLMs based on a reproductive Exponential Dispersion Model (EDM) aim at
    fitting and predicting the mean of the target y as mu=h(X*w).
    The fit minimizes the following objective function with L2 regularization::

            1/(2*sum(s)) * deviance(y, h(X*w); s) + 1/2 * alpha * ||w||_2^2

    with inverse link function h and s=sample_weight. Note that for
    ``sample_weight=None``, one has s_i=1 and sum(s)=n_samples).

    Read more in the :ref:`User Guide <Generalized_linear_regression>`.

    Parameters
    ----------
    power : float (default=0)
            The variance power: :math:`v(\\mu) = \\mu^{power}`.
            For ``0<power<1``, no distribution exists.

            Special cases are:

            ===== ================
            Power Distribution
            ===== ================
            0     Normal
            1     Poisson
            (1,2) Compound Poisson
            2     Gamma
            3     Inverse Gaussian

    alpha : float, optional (default=1)
        Constant that multiplies the penalty terms and thus determines the
        regularization strength.
        See the notes for the exact mathematical meaning of this
        parameter.``alpha = 0`` is equivalent to unpenalized GLMs. In this
        case, the design matrix X must have full column rank
        (no collinearities).

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    solver : {'lbfgs'}, optional (default='lbfgs')
        Algorithm to use in the optimization problem:

        'lbfgs'
            Calls scipy's L-BFGS-B optimizer.

    max_iter : int, optional (default=100)
        The maximal number of iterations for solver algorithms.

    tol : float, optional (default=1e-4)
        Stopping criterion. For the lbfgs solver,
        the iteration will stop when ``max{|g_i|, i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative) of
        the objective function.

    warm_start : boolean, optional (default=False)
        If set to ``True``, reuse the solution of the previous call to ``fit``
        as initialization for ``coef_`` and ``intercept_`` .

    copy_X : boolean, optional, (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    verbose : int, optional (default=0)
        For the lbfgs solver set verbose to any positive number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the linear predictor (X*coef_+intercept_) in
        the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    n_iter_ : int
        Actual number of iterations used in solver.

    Notes
    -----
    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments to be :math:`E[Y_i]=\\mu_i=h((Xw)_i)` and
    :math:`Var[Y_i]=\\frac{\\phi}{s_i} v(\\mu_i)`. The unit variance function
    :math:`v(\\mu_i)` is a property of and given by the specific EDM, see
    :ref:`User Guide <Generalized_linear_regression>`.

    The parameters :math:`w` (`coef_` and `intercept_`) are estimated by
    minimizing the deviance plus penalty term, which is equivalent to
    (penalized) maximum likelihood estimation.


    References
    ----------
    .. McCullagh, Peter; Nelder, John (1989). Generalized Linear Models,
       Second Edition. Boca Raton: Chapman and Hall/CRC. ISBN 0-412-31760-5.

    .. Jørgensen, B. (1992). The theory of exponential dispersion models
       and analysis of deviance. Monografias de matemática, no. 51.  See also
       `Exponential dispersion model.
       <https://en.wikipedia.org/wiki/Exponential_dispersion_model>`_
    """
    def __init__(self, power=0.0, alpha=1.0, fit_intercept=True, link='log',
                 solver='lbfgs', max_iter=100, tol=1e-4, warm_start=False,
                 copy_X=True, check_input=True, verbose=0):

        super().__init__(alpha=alpha, fit_intercept=fit_intercept,
                         family=TweedieDistribution(power=power), link=link,
                         solver=solver, max_iter=max_iter, tol=tol,
                         warm_start=warm_start, copy_X=copy_X, verbose=verbose)

    @property
    def family(self):
        dist = TweedieDistribution(power=self.power)
        # TODO: make the returned object immutable
        return dist

    @family.setter
    def family(self, value):
        if isinstance(value, TweedieDistribution):
            self.power = value.power
        else:
            raise TypeError("TweedieRegressor.family must be of type "
                            "TweedieDistribution!")
