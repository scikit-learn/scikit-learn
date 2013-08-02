from __future__ import print_function
#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# Licence: BSD 3 clause

import numpy as np
from scipy import linalg, optimize, rand, sparse, stats

from ..base import BaseEstimator, RegressorMixin, ClassifierMixin
from ..metrics.pairwise import manhattan_distances
from ..utils import array2d, check_random_state
from . import regression_models as regression
from . import correlation_models as correlation
from ..utils.extmath import logsumexp

MACHINE_EPSILON = np.finfo(np.double).eps
if hasattr(linalg, 'solve_triangular'):
    # only in scipy since 0.9
    solve_triangular = linalg.solve_triangular
else:
    # slower, but works
    def solve_triangular(x, y, lower=True):
        return linalg.solve(x, y)

def inv_triangular(A, lower = False):
    """
    Compute the inverse of triangular matrix
    
    Parameters
    ----------
    
    A - Triangular matrix
    lower - if A is lower-triangular
    
    Returns
    -------
    A_inv - Inverse of A (equivalant to linalg.inv(A))
    
    """
    
    A_inv = None
    # see if available use scipy.linalg.lapack.clapack.dtrtri instead
    if hasattr(linalg, 'lapack'):
        if hasattr(linalg.lapack, 'clapack'):
            dtrtri = linalg.lapack.clapack.dtrtri
            A = np.double(A)
            A_inv, info = dtrtri(A, lower = lower)
            
    if A_inv is None:
        A_inv = solve_triangular(A, np.eye(len(A)), lower = lower)
    
    return A_inv

def l1_cross_distances(X):
    """
    Computes the nonzero componentwise L1 cross-distances between the vectors
    in X.

    Parameters
    ----------

    X: array_like
        An array with shape (n_samples, n_features)

    Returns
    -------

    D: array with shape (n_samples * (n_samples - 1) / 2, n_features)
        The array of componentwise L1 cross-distances.

    ij: arrays with shape (n_samples * (n_samples - 1) / 2, 2)
        The indices i and j of the vectors in X associated to the cross-
        distances in D: D[k] = np.abs(X[ij[k, 0]] - Y[ij[k, 1]]).
    """
    X = array2d(X)
    n_samples, n_features = X.shape
    n_nonzero_cross_dist = n_samples * (n_samples - 1) / 2
    ij = np.zeros((n_nonzero_cross_dist, 2), dtype=np.int)
    D = np.zeros((n_nonzero_cross_dist, n_features))
    ll_1 = 0
    for k in range(n_samples - 1):
        ll_0 = ll_1
        ll_1 = ll_0 + n_samples - k - 1
        ij[ll_0:ll_1, 0] = k
        ij[ll_0:ll_1, 1] = np.arange(k + 1, n_samples)
        D[ll_0:ll_1] = np.abs(X[k] - X[(k + 1):n_samples])

    return D, ij.astype(np.int)


class GaussianProcess(BaseEstimator, RegressorMixin):
    """The Gaussian Process model class.

    Parameters
    ----------
    regr : string or callable, optional
        A regression function returning an array of outputs of the linear
        regression functional basis. The number of observations n_samples
        should be greater than the size p of this basis.
        Default assumes a simple constant regression trend.
        Available built-in regression models are::

            'constant', 'linear', 'quadratic'

    corr : string or callable, optional
        A stationary autocorrelation function returning the autocorrelation
        between two points x and x'.
        Default assumes a squared-exponential autocorrelation model.
        Built-in correlation models are::

            'absolute_exponential', 'squared_exponential',
            'generalized_exponential', 'cubic', 'linear'

    beta0 : double array_like, optional
        The regression weight vector to perform Ordinary Kriging (OK).
        Default assumes Universal Kriging (UK) so that the vector beta of
        regression weights is estimated using the maximum likelihood
        principle.

    storage_mode : string, optional
        A string specifying whether the Cholesky decomposition of the
        correlation matrix should be stored in the class (storage_mode =
        'full') or not (storage_mode = 'light').
        Default assumes storage_mode = 'full', so that the
        Cholesky decomposition of the correlation matrix is stored.
        This might be a useful parameter when one is not interested in the
        MSE and only plan to estimate the BLUP, for which the correlation
        matrix is not required.

    verbose : boolean, optional
        A boolean specifying the verbose level.
        Default is verbose = False.

    theta0 : double array_like, optional
        An array with shape (n_features, ) or (1, ).
        The parameters in the autocorrelation model.
        If thetaL and thetaU are also specified, theta0 is considered as
        the starting point for the maximum likelihood estimation of the
        best set of parameters.
        Default assumes isotropic autocorrelation model with theta0 = 1e-1.

    thetaL : double array_like, optional
        An array with shape matching theta0's.
        Lower bound on the autocorrelation parameters for maximum
        likelihood estimation.
        Default is None, so that it skips maximum likelihood estimation and
        it uses theta0.

    thetaU : double array_like, optional
        An array with shape matching theta0's.
        Upper bound on the autocorrelation parameters for maximum
        likelihood estimation.
        Default is None, so that it skips maximum likelihood estimation and
        it uses theta0.

    normalize : boolean, optional
        Input X and observations y are centered and reduced wrt
        means and standard deviations estimated from the n_samples
        observations provided.
        Default is normalize = True so that data is normalized to ease
        maximum likelihood estimation.

    nugget : double or ndarray, optional
        Introduce a nugget effect to allow smooth predictions from noisy
        data.  If nugget is an ndarray, it must be the same length as the
        number of data points used for the fit.
        The nugget is added to the diagonal of the assumed training covariance;
        in this way it acts as a Tikhonov regularization in the problem.  In
        the special case of the squared exponential correlation function, the
        nugget mathematically represents the variance of the input values.
        Default assumes a nugget close to machine precision for the sake of
        robustness (nugget = 10. * MACHINE_EPSILON).

    optimizer : string, optional
        A string specifying the optimization algorithm to be used.
        Default uses 'fmin_cobyla' algorithm from scipy.optimize.
        Available optimizers are::

            'fmin_cobyla', 'Welch'

        'Welch' optimizer is dued to Welch et al., see reference [WBSWM1992]_.
        It consists in iterating over several one-dimensional optimizations
        instead of running one single multi-dimensional optimization.

    random_start : int, optional
        The number of times the Maximum Likelihood Estimation should be
        performed from a random starting point.
        The first MLE always uses the specified starting point (theta0),
        the next starting points are picked at random according to an
        exponential distribution (log-uniform on [thetaL, thetaU]).
        Default does not use random starting point (random_start = 1).

    random_state: integer or numpy.RandomState, optional
        The generator used to shuffle the sequence of coordinates of theta in
        the Welch optimizer. If an integer is given, it fixes the seed.
        Defaults to the global numpy random number generator.


    Attributes
    ----------
    `theta_`: array
        Specified theta OR the best set of autocorrelation parameters (the \
        sought maximizer of the reduced likelihood function).

    `reduced_likelihood_function_value_`: array
        The optimal reduced likelihood function value.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.gaussian_process import GaussianProcess
    >>> X = np.array([[1., 3., 5., 6., 7., 8.]]).T
    >>> y = (X * np.sin(X)).ravel()
    >>> gp = GaussianProcess(theta0=0.1, thetaL=.001, thetaU=1.)
    >>> gp.fit(X, y)                                      # doctest: +ELLIPSIS
    GaussianProcess(beta0=None...
            ...

    Notes
    -----
    The presentation implementation is based on a translation of the DACE
    Matlab toolbox, see reference [NLNS2002]_.

    References
    ----------

    .. [NLNS2002] `H.B. Nielsen, S.N. Lophaven, H. B. Nielsen and J.
        Sondergaard.  DACE - A MATLAB Kriging Toolbox.` (2002)
        http://www2.imm.dtu.dk/~hbn/dace/dace.pdf

    .. [WBSWM1992] `W.J. Welch, R.J. Buck, J. Sacks, H.P. Wynn, T.J. Mitchell,
        and M.D.  Morris (1992). Screening, predicting, and computer
        experiments.  Technometrics, 34(1) 15--25.`
        http://www.jstor.org/pss/1269548
    """

    _regression_types = {
        'constant': regression.constant,
        'linear': regression.linear,
        'quadratic': regression.quadratic}

    _correlation_types = {
        'absolute_exponential': correlation.absolute_exponential,
        'squared_exponential': correlation.squared_exponential,
        'generalized_exponential': correlation.generalized_exponential,
        'cubic': correlation.cubic,
        'linear': correlation.linear}

    _optimizer_types = [
        'fmin_cobyla',
        'Welch']

    def __init__(self, regr='constant', corr='squared_exponential', beta0=None,
                 storage_mode='full', verbose=False, theta0=1e-1,
                 thetaL=None, thetaU=None, optimizer='fmin_cobyla',
                 random_start=1, normalize=True,
                 nugget=10. * MACHINE_EPSILON, random_state=None):

        self.regr = regr
        self.corr = corr
        self.beta0 = beta0
        self.storage_mode = storage_mode
        self.verbose = verbose
        self.theta0 = theta0
        self.thetaL = thetaL
        self.thetaU = thetaU
        self.normalize = normalize
        self.nugget = nugget
        self.optimizer = optimizer
        self.random_start = random_start
        self.random_state = random_state

    def fit(self, X, y):
        """
        The Gaussian Process model fitting method.

        Parameters
        ----------
        X : double array_like
            An array with shape (n_samples, n_features) with the input at which
            observations were made.

        y : double array_like
            An array with shape (n_samples, ) with the observations of the
            scalar output to be predicted.

        Returns
        -------
        gp : self
            A fitted Gaussian Process model object awaiting data to perform
            predictions.
        """
        # Run input checks
        self._check_params()

        self.random_state = check_random_state(self.random_state)

        # Force data to 2D numpy.array
        X = array2d(X)
        y = np.asarray(y).ravel()[:, np.newaxis]

        # Check shapes of DOE & observations
        n_samples_X, n_features = X.shape
        n_samples_y = y.shape[0]

        if n_samples_X != n_samples_y:
            raise ValueError("X and y must have the same number of rows.")
        else:
            n_samples = n_samples_X

        # Run input checks
        self._check_params(n_samples)

        # Normalize data or don't
        if self.normalize:
            X_mean = np.mean(X, axis=0)
            X_std = np.std(X, axis=0)
            y_mean = np.mean(y, axis=0)
            y_std = np.std(y, axis=0)
            X_std[X_std == 0.] = 1.
            y_std[y_std == 0.] = 1.
            # center and scale X if necessary
            X = (X - X_mean) / X_std
            y = (y - y_mean) / y_std
        else:
            X_mean = np.zeros(1)
            X_std = np.ones(1)
            y_mean = np.zeros(1)
            y_std = np.ones(1)

        # Calculate matrix of distances D between samples
        D, ij = l1_cross_distances(X)
        if (np.min(np.sum(D, axis=1)) == 0.
                and self.corr != correlation.pure_nugget):
            raise Exception("Multiple input features cannot have the same"
                            " value.")

        # Regression matrix and parameters
        F = self.regr(X)
        n_samples_F = F.shape[0]
        if F.ndim > 1:
            p = F.shape[1]
        else:
            p = 1
        if n_samples_F != n_samples:
            raise Exception("Number of rows in F and X do not match. Most "
                            "likely something is going wrong with the "
                            "regression model.")
        if p > n_samples_F:
            raise Exception(("Ordinary least squares problem is undetermined "
                             "n_samples=%d must be greater than the "
                             "regression model size p=%d.") % (n_samples, p))
        if self.beta0 is not None:
            if self.beta0.shape[0] != p:
                raise Exception("Shapes of beta0 and F do not match.")

        # Set attributes
        self.X = X
        self.y = y
        self.D = D
        self.ij = ij
        self.F = F
        self.X_mean, self.X_std = X_mean, X_std
        self.y_mean, self.y_std = y_mean, y_std

        # Determine Gaussian Process model parameters
        if self.thetaL is not None and self.thetaU is not None:
            # Maximum Likelihood Estimation of the parameters
            if self.verbose:
                print("Performing Maximum Likelihood Estimation of the "
                      "autocorrelation parameters...")
            self.theta_, self.reduced_likelihood_function_value_, par = \
                self._arg_max_reduced_likelihood_function()
            if np.isinf(self.reduced_likelihood_function_value_):
                raise Exception("Bad parameter region. "
                                "Try increasing upper bound")

        else:
            # Given parameters
            if self.verbose:
                print("Given autocorrelation parameters. "
                      "Computing Gaussian Process model parameters...")
            self.theta_ = self.theta0
            self.reduced_likelihood_function_value_, par = \
                self.reduced_likelihood_function()
            if np.isinf(self.reduced_likelihood_function_value_):
                raise Exception("Bad point. Try increasing theta0.")

        self.beta = par['beta']
        self.gamma = par['gamma']
        self.sigma2 = par['sigma2']
        self.C = par['C']
        self.Ft = par['Ft']
        self.G = par['G']

        if self.storage_mode == 'light':
            # Delete heavy data (it will be computed again if required)
            # (it is required only when MSE is wanted in self.predict)
            if self.verbose:
                print("Light storage mode specified. "
                      "Flushing autocorrelation matrix...")
            self.D = None
            self.ij = None
            self.F = None
            self.C = None
            self.Ft = None
            self.G = None

        return self

    def predict(self, X, eval_MSE=False, batch_size=None):
        """
        This function evaluates the Gaussian Process model at x.

        Parameters
        ----------
        X : array_like
            An array with shape (n_eval, n_features) giving the point(s) at
            which the prediction(s) should be made.

        eval_MSE : boolean, optional
            A boolean specifying whether the Mean Squared Error should be
            evaluated or not.
            Default assumes evalMSE = False and evaluates only the BLUP (mean
            prediction).

        batch_size : integer, optional
            An integer giving the maximum number of points that can be
            evaluated simultaneously (depending on the available memory).
            Default is None so that all given points are evaluated at the same
            time.

        Returns
        -------
        y : array_like
            An array with shape (n_eval, ) with the Best Linear Unbiased
            Prediction at x.

        MSE : array_like, optional (if eval_MSE == True)
            An array with shape (n_eval, ) with the Mean Squared Error at x.
        """

        # Check input shapes
        X = array2d(X)
        n_eval, n_features_X = X.shape
        n_samples, n_features = self.X.shape

        # Run input checks
        self._check_params(n_samples)

        if n_features_X != n_features:
            raise ValueError(("The number of features in X (X.shape[1] = %d) "
                             "should match the sample size used for fit() "
                             "which is %d.") % (n_features_X, n_features))

        if batch_size is None:
            # No memory management
            # (evaluates all given points in a single batch run)

            # Normalize input
            X = (X - self.X_mean) / self.X_std

            # Initialize output
            y = np.zeros(n_eval)
            if eval_MSE:
                MSE = np.zeros(n_eval)

            # Get pairwise componentwise L1-distances to the input training set
            dx = manhattan_distances(X, Y=self.X, sum_over_features=False)
            # Get regression function and correlation
            f = self.regr(X)
            r = self.corr(self.theta_, dx).reshape(n_eval, n_samples)

            # Scaled predictor
            y_ = np.dot(f, self.beta) + np.dot(r, self.gamma)

            # Predictor
            y = (self.y_mean + self.y_std * y_).ravel()

            # Mean Squared Error
            if eval_MSE:
                C = self.C
                if C is None:
                    # Light storage mode (need to recompute C, F, Ft and G)
                    if self.verbose:
                        print("This GaussianProcess used 'light' storage mode "
                              "at instantiation. Need to recompute "
                              "autocorrelation matrix...")
                    reduced_likelihood_function_value, par = \
                        self.reduced_likelihood_function()
                    self.C = par['C']
                    self.Ft = par['Ft']
                    self.G = par['G']

                rt = solve_triangular(self.C, r.T, lower=True)

                if self.beta0 is None:
                    # Universal Kriging
                    u = solve_triangular(self.G.T,
                                         np.dot(self.Ft.T, rt) - f.T)
                else:
                    # Ordinary Kriging
                    u = np.zeros(y.shape)

                MSE = self.sigma2 * (1. - (rt ** 2.).sum(axis=0)
                                     + (u ** 2.).sum(axis=0))

                # Mean Squared Error might be slightly negative depending on
                # machine precision: force to zero!
                MSE[MSE < 0.] = 0.

                return y, MSE

            else:

                return y

        else:
            # Memory management

            if type(batch_size) is not int or batch_size <= 0:
                raise Exception("batch_size must be a positive integer")

            if eval_MSE:

                y, MSE = np.zeros(n_eval), np.zeros(n_eval)
                for k in range(max(1, n_eval / batch_size)):
                    batch_from = k * batch_size
                    batch_to = min([(k + 1) * batch_size + 1, n_eval + 1])
                    y[batch_from:batch_to], MSE[batch_from:batch_to] = \
                        self.predict(X[batch_from:batch_to],
                                     eval_MSE=eval_MSE, batch_size=None)

                return y, MSE

            else:

                y = np.zeros(n_eval)
                for k in range(max(1, n_eval / batch_size)):
                    batch_from = k * batch_size
                    batch_to = min([(k + 1) * batch_size + 1, n_eval + 1])
                    y[batch_from:batch_to] = \
                        self.predict(X[batch_from:batch_to],
                                     eval_MSE=eval_MSE, batch_size=None)

                return y

    def reduced_likelihood_function(self, theta=None):
        """
        This function determines the BLUP parameters and evaluates the reduced
        likelihood function for the given autocorrelation parameters theta.

        Maximizing this function wrt the autocorrelation parameters theta is
        equivalent to maximizing the likelihood of the assumed joint Gaussian
        distribution of the observations y evaluated onto the design of
        experiments X.

        Parameters
        ----------
        theta : array_like, optional
            An array containing the autocorrelation parameters at which the
            Gaussian Process model parameters should be determined.
            Default uses the built-in autocorrelation parameters
            (ie ``theta = self.theta_``).

        Returns
        -------
        reduced_likelihood_function_value : double
            The value of the reduced likelihood function associated to the
            given autocorrelation parameters theta.

        par : dict
            A dictionary containing the requested Gaussian Process model
            parameters:

                sigma2
                        Gaussian Process variance.
                beta
                        Generalized least-squares regression weights for
                        Universal Kriging or given beta0 for Ordinary
                        Kriging.
                gamma
                        Gaussian Process weights.
                C
                        Cholesky decomposition of the correlation matrix [R].
                Ft
                        Solution of the linear equation system : [R] x Ft = F
                G
                        QR decomposition of the matrix Ft.
        """

        if theta is None:
            # Use built-in autocorrelation parameters
            theta = self.theta_

        # Initialize output
        reduced_likelihood_function_value = - np.inf
        par = {}

        # Retrieve data
        n_samples = self.X.shape[0]
        D = self.D
        ij = self.ij
        F = self.F

        if D is None:
            # Light storage mode (need to recompute D, ij and F)
            D, ij = l1_cross_distances(self.X)
            if (np.min(np.sum(D, axis=1)) == 0.
                    and self.corr != correlation.pure_nugget):
                raise Exception("Multiple X are not allowed")
            F = self.regr(self.X)

        # Set up R
        r = self.corr(theta, D)
        R = np.eye(n_samples) * (1. + self.nugget)
        R[ij[:, 0], ij[:, 1]] = r
        R[ij[:, 1], ij[:, 0]] = r

        # Cholesky decomposition of R
        try:
            C = linalg.cholesky(R, lower=True)
        except linalg.LinAlgError:
            return reduced_likelihood_function_value, par

        # Get generalized least squares solution
        Ft = solve_triangular(C, F, lower=True)
        try:
            Q, G = linalg.qr(Ft, econ=True)
        except:
            #/usr/lib/python2.6/dist-packages/scipy/linalg/decomp.py:1177:
            # DeprecationWarning: qr econ argument will be removed after scipy
            # 0.7. The economy transform will then be available through the
            # mode='economic' argument.
            Q, G = linalg.qr(Ft, mode='economic')
            pass

        sv = linalg.svd(G, compute_uv=False)
        rcondG = sv[-1] / sv[0]
        if rcondG < 1e-10:
            # Check F
            sv = linalg.svd(F, compute_uv=False)
            condF = sv[0] / sv[-1]
            if condF > 1e15:
                raise Exception("F is too ill conditioned. Poor combination "
                                "of regression model and observations.")
            else:
                # Ft is too ill conditioned, get out (try different theta)
                return reduced_likelihood_function_value, par

        Yt = solve_triangular(C, self.y, lower=True)
        if self.beta0 is None:
            # Universal Kriging
            beta = solve_triangular(G, np.dot(Q.T, Yt))
        else:
            # Ordinary Kriging
            beta = np.array(self.beta0)

        rho = Yt - np.dot(Ft, beta)
        sigma2 = (rho ** 2.).sum(axis=0) / n_samples
        # The determinant of R is equal to the squared product of the diagonal
        # elements of its Cholesky decomposition C
        detR = (np.diag(C) ** (2. / n_samples)).prod()

        # Compute/Organize output
        reduced_likelihood_function_value = - sigma2.sum() * detR
        par['sigma2'] = sigma2 * self.y_std ** 2.
        par['beta'] = beta
        par['gamma'] = solve_triangular(C.T, rho)
        par['C'] = C
        par['Ft'] = Ft
        par['G'] = G

        return reduced_likelihood_function_value, par

    def _arg_max_reduced_likelihood_function(self):
        """
        This function estimates the autocorrelation parameters theta as the
        maximizer of the reduced likelihood function.
        (Minimization of the opposite reduced likelihood function is used for
        convenience)

        Parameters
        ----------
        self : All parameters are stored in the Gaussian Process model object.

        Returns
        -------
        optimal_theta : array_like
            The best set of autocorrelation parameters (the sought maximizer of
            the reduced likelihood function).

        optimal_reduced_likelihood_function_value : double
            The optimal reduced likelihood function value.

        optimal_par : dict
            The BLUP parameters associated to thetaOpt.
        """

        # Initialize output
        best_optimal_theta = []
        best_optimal_rlf_value = []
        best_optimal_par = []

        if self.verbose:
            print("The chosen optimizer is: " + str(self.optimizer))
            if self.random_start > 1:
                print(str(self.random_start) + " random starts are required.")

        percent_completed = 0.

        # Force optimizer to fmin_cobyla if the model is meant to be isotropic
        if self.optimizer == 'Welch' and self.theta0.size == 1:
            self.optimizer = 'fmin_cobyla'

        if self.optimizer == 'fmin_cobyla':

            def minus_reduced_likelihood_function(log10t):
                return - self.reduced_likelihood_function(
                    theta=10. ** log10t)[0]

            constraints = []
            for i in range(self.theta0.size):
                constraints.append(lambda log10t:
                                   log10t[i] - np.log10(self.thetaL[0, i]))
                constraints.append(lambda log10t:
                                   np.log10(self.thetaU[0, i]) - log10t[i])

            for k in range(self.random_start):

                if k == 0:
                    # Use specified starting point as first guess
                    theta0 = self.theta0
                else:
                    # Generate a random starting point log10-uniformly
                    # distributed between bounds
                    log10theta0 = np.log10(self.thetaL) \
                        + rand(self.theta0.size).reshape(self.theta0.shape) \
                        * np.log10(self.thetaU / self.thetaL)
                    theta0 = 10. ** log10theta0

                # Run Cobyla
                try:
                    log10_optimal_theta = \
                        optimize.fmin_cobyla(minus_reduced_likelihood_function,
                                             np.log10(theta0), constraints,
                                             iprint=0)
                except ValueError as ve:
                    print("Optimization failed. Try increasing the ``nugget``")
                    raise ve

                optimal_theta = 10. ** log10_optimal_theta
                optimal_minus_rlf_value, optimal_par = \
                    self.reduced_likelihood_function(theta=optimal_theta)
                optimal_rlf_value = - optimal_minus_rlf_value

                # Compare the new optimizer to the best previous one
                if k > 0:
                    if optimal_rlf_value > best_optimal_rlf_value:
                        best_optimal_rlf_value = optimal_rlf_value
                        best_optimal_par = optimal_par
                        best_optimal_theta = optimal_theta
                else:
                    best_optimal_rlf_value = optimal_rlf_value
                    best_optimal_par = optimal_par
                    best_optimal_theta = optimal_theta
                if self.verbose and self.random_start > 1:
                    if (20 * k) / self.random_start > percent_completed:
                        percent_completed = (20 * k) / self.random_start
                        print("%s completed" % (5 * percent_completed))

            optimal_rlf_value = best_optimal_rlf_value
            optimal_par = best_optimal_par
            optimal_theta = best_optimal_theta

        elif self.optimizer == 'Welch':

            # Backup of the given atrributes
            theta0, thetaL, thetaU = self.theta0, self.thetaL, self.thetaU
            corr = self.corr
            verbose = self.verbose

            # This will iterate over fmin_cobyla optimizer
            self.optimizer = 'fmin_cobyla'
            self.verbose = False

            # Initialize under isotropy assumption
            if verbose:
                print("Initialize under isotropy assumption...")
            self.theta0 = array2d(self.theta0.min())
            self.thetaL = array2d(self.thetaL.min())
            self.thetaU = array2d(self.thetaU.max())
            theta_iso, optimal_rlf_value_iso, par_iso = \
                self._arg_max_reduced_likelihood_function()
            optimal_theta = theta_iso + np.zeros(theta0.shape)

            # Iterate over all dimensions of theta allowing for anisotropy
            if verbose:
                print("Now improving allowing for anisotropy...")
            for i in self.random_state.permutation(theta0.size):
                if verbose:
                    print("Proceeding along dimension %d..." % (i + 1))
                self.theta0 = array2d(theta_iso)
                self.thetaL = array2d(thetaL[0, i])
                self.thetaU = array2d(thetaU[0, i])

                def corr_cut(t, d):
                    return corr(array2d(np.hstack([optimal_theta[0][0:i],
                                                   t[0],
                                                   optimal_theta[0][(i + 1)::]]
                                                  )), d)

                self.corr = corr_cut
                optimal_theta[0, i], optimal_rlf_value, optimal_par = \
                    self._arg_max_reduced_likelihood_function()

            # Restore the given atrributes
            self.theta0, self.thetaL, self.thetaU = theta0, thetaL, thetaU
            self.corr = corr
            self.optimizer = 'Welch'
            self.verbose = verbose

        else:

            raise NotImplementedError("This optimizer ('%s') is not "
                                      "implemented yet. Please contribute!"
                                      % self.optimizer)

        return optimal_theta, optimal_rlf_value, optimal_par

    def _check_params(self, n_samples=None):

        # Check regression model
        if not callable(self.regr):
            if self.regr in self._regression_types:
                self.regr = self._regression_types[self.regr]
            else:
                raise ValueError("regr should be one of %s or callable, "
                                 "%s was given."
                                 % (self._regression_types.keys(), self.regr))

        # Check regression weights if given (Ordinary Kriging)
        if self.beta0 is not None:
            self.beta0 = array2d(self.beta0)
            if self.beta0.shape[1] != 1:
                # Force to column vector
                self.beta0 = self.beta0.T

        # Check correlation model
        if not callable(self.corr):
            if self.corr in self._correlation_types:
                self.corr = self._correlation_types[self.corr]
            else:
                raise ValueError("corr should be one of %s or callable, "
                                 "%s was given."
                                 % (self._correlation_types.keys(), self.corr))

        # Check storage mode
        if self.storage_mode != 'full' and self.storage_mode != 'light':
            raise ValueError("Storage mode should either be 'full' or "
                             "'light', %s was given." % self.storage_mode)

        # Check correlation parameters
        self.theta0 = array2d(self.theta0)
        lth = self.theta0.size

        if self.thetaL is not None and self.thetaU is not None:
            self.thetaL = array2d(self.thetaL)
            self.thetaU = array2d(self.thetaU)
            if self.thetaL.size != lth or self.thetaU.size != lth:
                raise ValueError("theta0, thetaL and thetaU must have the "
                                 "same length.")
            if np.any(self.thetaL <= 0) or np.any(self.thetaU < self.thetaL):
                raise ValueError("The bounds must satisfy O < thetaL <= "
                                 "thetaU.")

        elif self.thetaL is None and self.thetaU is None:
            if np.any(self.theta0 <= 0):
                raise ValueError("theta0 must be strictly positive.")

        elif self.thetaL is None or self.thetaU is None:
            raise ValueError("thetaL and thetaU should either be both or "
                             "neither specified.")

        # Force verbose type to bool
        self.verbose = bool(self.verbose)

        # Force normalize type to bool
        self.normalize = bool(self.normalize)

        # Check nugget value
        self.nugget = np.asarray(self.nugget)
        if np.any(self.nugget) < 0.:
            raise ValueError("nugget must be positive or zero.")
        if (n_samples is not None
                and self.nugget.shape not in [(), (n_samples,)]):
            raise ValueError("nugget must be either a scalar "
                             "or array of length n_samples.")

        # Check optimizer
        if not self.optimizer in self._optimizer_types:
            raise ValueError("optimizer should be one of %s"
                             % self._optimizer_types)

        # Force random_start type to int
        self.random_start = int(self.random_start)


class GaussianProcessClassifier(BaseEstimator, ClassifierMixin):
    """The Gaussian Process Classifier class.

    Parameters
    ----------
    method : string, optional
        'laplace' - Use laplace approximation (default)
    
    verbose : boolean, optional
        A boolean specifying the verbose level.
        Default is verbose = False.
        
    normalize : boolean, optional
        Input X is centered and reduced wrt
        means and standard deviations estimated from the n_samples
        observations provided.
        Default is normalize = True so that data is normalized to ease
        maximum likelihood estimations.
        
    nugget : double or ndarray, optional
        Introduce a nugget effect to allow smooth predictions from noisy
        data.  If nugget is an ndarray, it must be the same length as the
        number of data points used for the fit.
        The nugget is added to the diagonal of the assumed training covariance;
        in this way it acts as a Tikhonov regularization in the problem.  In
        the special case of the squared exponential correlation function, the
        nugget mathematically represents the variance of the input values.
        Default assumes a nugget close to machine precision for the sake of
        robustness (nugget = 10. * MACHINE_EPSILON).
    
    max_iter : int, optional
        Maximum number of iterations in newton algorithm in MAP
        Default: 0 means 10 * n_samples * n_classes
    
    mc_iter: int, optional
        Number of Monte Carlo (MC) samples to take to estimate test
        Default: None, do enough sampling for 95% confidence interval
        
    tol: float
        Convergence tolerance in MAP algorithm
        Default: 10. * MACHINE_EPSILON

    theta0 : double array_like, optional
        An array with shape (n_classes, n_features, ) or (n_classes,) or (1, ).
        The parameters in the autocorrelation model.
        Default assumes isotropic autocorrelation model with theta0 = 1e-1 for all classes and all features.

    Attributes
    ----------
    `theta_`: array
        Specified theta OR the best set of autocorrelation parameters (the \
        sought maximizer of the reduced likelihood function).
        
    `means_` : array, shape (`n_classes`, `n_samples`)
        Mean parameters for each latent function value.

    `reduced_likelihood_function_value_`: array
        The optimal reduced likelihood function value.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.gaussian_process import GaussianProcessClassifier
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([0, 0, 0, 1, 1, 1])
    >>> gpc = GaussianProcessClassifier()
    >>> gpc.fit(X, y)
    >>> print(gpc.predict([[-0.8, -1]]))

    Notes
    -----
    The presentation implementation is based on the Gaussian 
    Processes for Machine Learning (Rassmussen Williams), 
    see reference [RASWIL06]

    References
    ----------

    .. [RASWIL06] C. E. Rasmussen & C. K. I. Williams, G, `Gaussian Processes for Machine Learning`, 
                the MIT Press, 2006, ISBN 026218253X
        http://www.gaussianprocess.org/gpml/

    """

    _correlation_types = {
        'absolute_exponential': correlation.absolute_exponential,
        'squared_exponential': correlation.squared_exponential,
        'generalized_exponential': correlation.generalized_exponential,
        'cubic': correlation.cubic,
        'linear': correlation.linear}

    def __init__(self, method='laplace', corr='squared_exponential', verbose=False,
                 normalize=True, max_iter=0, mc_iter = None, tol=10 * MACHINE_EPSILON, theta0=1e-1, 
                 nugget=10. * MACHINE_EPSILON):

        self.method = method
        self.verbose = verbose
        self.normalize = normalize
        self.max_iter = max_iter
        self.mc_iter = mc_iter
        self.tol = tol
        self.theta0 = theta0
        self.nugget = nugget 
        self.corr = corr

    def fit(self, X, y):
        """
        The Gaussian Process model fitting method.

        Parameters
        ----------
        X : double array_like
            An array with shape (n_samples, n_features) with the input at which
            observations were made.

        y : double array_like
            An array with shape (n_samples, ) with the class labels

        Returns
        -------
        gpc : self
            A fitted Gaussian Process Classifier model object awaiting data to perform
            predictions.
        """

        # TODO: maximize hyperparameter theta using arg_max_reduced_likelyhood
        # TODO: Implement expectation maximization 'EM'/'EP' method for multiclass
        # TODO: Can use PyMC for multiclass posterior approximation
        # TODO: Implement variational approximation to the posterior
        # TODO: better initialize f_vec
        # TODO: investigate reformulating to use fmin_ncg but without the need to K^-1 
        # TODO: implement numerical integration for prediction to use if mc_iter=0
        
        # Force data to 2D numpy.array
        X = array2d(X)

        # Check shapes of DOE & observations
        n_samples_X, n_features = X.shape
        n_samples_y = y.shape[0]

        if n_samples_X != n_samples_y:
            raise ValueError("X and y must have the same number of rows.")
        else:
            n_samples = n_samples_X

        n_classes = len(np.unique(y))
        
        # Run input checks
        self._check_params(n_classes, n_samples)
        
        # Create empty cache
        self._remove_cache()
        
        # Normalize data or don't
        if self.normalize:
            X_mean = np.mean(X, axis=0)
            X_std = np.std(X, axis=0)
            X_std[X_std == 0.] = 1.
            # center and scale X if necessary
            X = (X - X_mean) / X_std
        else:
            X_mean = np.zeros(1)
            X_std = np.ones(1)
        
        D, ij = l1_cross_distances(X)
        if (np.min(np.sum(D, axis=1)) == 0.
                and self.corr != correlation.pure_nugget):
            raise Exception("Multiple input features cannot have the same"
                            " value.")
        
        self._cache['D'] = D
        self._cache['ij'] = ij
        
        # 0/1 encoding of the labels     
        y_vec = np.zeros((n_classes, n_samples))
        for c in range(n_classes):
            y_vec[c, y == c] = 1 
        y_vec = y_vec.reshape((n_classes * n_samples, 1))
        
        # Keep model parameters
        self.X = X
        self.y = y_vec
        self.X_mean, self.X_std = X_mean, X_std
        self.n_classes = n_classes
        
        # Compute the maximum aposteriori value
        self.theta_ = self.theta0
        self.reduced_likelihood_function_value_, par = self.reduced_likelihood_function()
        
        # Keep some needed computed values
        self.means_ = par['f_vec'].reshape(n_classes, n_samples) # posterios mean
        self.E = par['E']
        self.F_inv = par['F_inv']
        self.pi_vec = par['pi_vec']
        
        # Remove cache
        self._remove_cache()
        
        return self
    
    def predict(self, X):
        """Predict label for data.

        Parameters
        ----------
        X : array-like, shape = [n_eval, n_features]

        Returns
        -------
        C : array, shape = (n_eval,) with predicted class labels
        """
        responsibilities = self.predict_proba(X)
        return responsibilities.argmax(axis=1)
    
    def predict_proba(self, X):
        """
        This function evaluates the Gaussian Process model at X.

        Parameters
        ----------
        X : array_like
            An array with shape (n_eval, n_features) giving the point(s) at
            which the prediction(s) should be made.

        Returns
        -------
        responsibilities : array-like, shape = (n_eval, n_classes)
            Returns the probability of the sample for each class.

        """

        # Check input shapes
        X = array2d(X)
        n_eval, n_features_X = X.shape
        n_samples, n_features = self.X.shape
        n_classes = self.n_classes
        
        # Run input checks
        self._check_params()

        if n_features_X != n_features:
            raise ValueError(("The number of features in X (X.shape[1] = %d) "
                             "should match the sample size used for fit() "
                             "which is %d.") % (n_features_X, n_features))

        # No memory management for now
        # (evaluates all given points in a single batch run)

        # Normalize input
        X = (X - self.X_mean) / self.X_std

        y_mat = self.y.reshape(n_classes, n_samples)
        pi_mat =  self.pi_vec.reshape(n_classes, n_samples)
         
        # Get pairwise componentwise L1-distances to the input training set
        dx = manhattan_distances(X, Y=self.X, sum_over_features=False)
        
        # test mean
        test_means_ = np.zeros((n_eval, n_classes))
        
        K_test = []
        for c in range(n_classes):
            # Test correlation model
            r = self.corr(self.theta_[c,:], dx).reshape(n_eval, n_samples)
            K_test.append(r)
            
            # Compute test mean
            tmp = (y_mat[c,:] - pi_mat[c,:]).reshape(n_samples, 1)
            test_means_[:, c] = np.dot(r, tmp).reshape(n_eval)
        
        # Compute test covariance
        test_covars_ = np.zeros((n_eval, n_classes, n_classes))
        D = np.zeros((1, n_features))
        for c in range(n_classes):
            # For each test point x_0, for each class c compute the vector of k_c(x_0, x_0)
            #  should be 1 for most interesting correlations
            r = self.corr(self.theta_[c,:], D)
            
            b = np.dot(self.E[c], K_test[c].T)
            # corrected based on errata 
            #  note: first term could be saved during training
            v = np.dot(np.dot(self.E[c], self.F_inv), b)
            
            for idx in range(n_eval): 
                for c_p in range(n_classes):
                    sigma = np.dot(K_test[c_p][idx, :], v[:,idx])
                    test_covars_[idx,c,c_p] = sigma
                test_covars_[idx][c,c] = test_covars_[idx][c,c] + r - np.dot(b[:, idx].T, K_test[c][idx, :])                       
        
        # decision
        p = np.zeros((n_eval, n_classes))
        # if montecarlo (MC) is specified then sample from the posterior
        if self.mc_iter is None:
            for idx in range(n_eval):
                sv = linalg.svd(test_covars_[idx], compute_uv=False)
                # only need to care about the component with maximum variance
                max_std_idx = np.argmax(sv)
                max_std = sv[max_std_idx]
                
                # Compute output standard deviation (perhaps should do this in a feedback loop)
                f_vec = np.random.multivariate_normal(test_means_[idx], test_covars_[idx], 100).T
                f_vec = f_vec.reshape(f_vec.size)
                pi_vec = self.soft_max(f_vec)
                pi_vec = pi_vec.reshape(n_classes, 100)
                
                out_std = np.std(pi_vec[max_std_idx,:])
                confidence = 1. / stats.norm.cdf((1 + .95)/2)
                N = int((max_std * confidence / out_std) ** 2)
                if self.verbose:
                    print "Computed MC samples for 95%% confidense is %s" % N
                
                if N > 100:
                    f_vec = np.random.multivariate_normal(test_means_[idx], test_covars_[idx], N - 100).T
                    f_vec = f_vec.reshape(f_vec.size)
                    pi_vec2 = self.soft_max(f_vec)
                    pi_vec2 = pi_vec2.reshape(n_classes, N - 100)
                    # Looka at them all
                    pi_vec = np.hstack((pi_vec, pi_vec2))
                    
                p[idx] = pi_vec.mean(axis = 1)
                    
        elif self.mc_iter > 0:
            # If a specific number is given, use it
            for idx in range(n_eval):
                f_vec = np.random.multivariate_normal(test_means_[idx], test_covars_[idx], self.mc_iter).T
                f_vec = f_vec.reshape(f_vec.size)
                pi_vec = self.soft_max(f_vec)
                pi_vec = pi_vec.reshape(n_classes, self.mc_iter)
                p[idx] = pi_vec.mean(axis = 1)
        else:
            # TODO: use a integrals calculation
            raise ValueError("Only monte carlo method implemented at this stage.")
                
        return p

    def soft_max(self, f_vec):
        """
        Soft-max decision function
        
        Parameters
        ----------
        f_vec : array_like shape=(n_classes * n_samples)
            latent function values
            
        Returns
        -------
        pi_vec : array_like
            An array with shape (n_classes * n_samples) with 
            squashed probabilities for decision function values
        """
        
        # loggit softmax
        loggit_soft_max = (lambda c,i,f: np.exp(f[c,i]) / np.exp(logsumexp(f[:,i])))
        n_classes = self.n_classes
        n_samples = f_vec.size / n_classes
        
        # Reshape is cheap
        f_vec = f_vec.reshape(n_classes, n_samples)
        # decision function values: - of size Cn
        pi_vec = np.zeros((n_classes, n_samples))
        for i in range(n_samples):
            for c in range(n_classes):
                pi_vec[c, i] = loggit_soft_max(c, i, f_vec)
        pi_vec = pi_vec.reshape((n_classes * n_samples, 1))
        return pi_vec

    def _remove_cache(self):
        """
        Remove cached values
        """
        self._cache = {'R': None, 'D': None, 'ij': None}
        
    def _get_R(self):
        """
        Create and cache R, return from cache if available
        Returns
        -------
            R - Sparse matrix of stacked identity matrices
        """
        
        n_samples = self.X.shape[0]
        n_classes = self.n_classes
        R_sparse = self._cache['R']
        if R_sparse is None:    
            # Matrix of stacked identity matrices
            R_sparse = []
            for c in range(n_classes):
                R = sparse.csr_matrix(np.eye(n_samples))
                R_sparse.append(R)
            R_sparse = sparse.vstack(R_sparse).tocsr()
            self._cache['R'] = R_sparse
        
        return R_sparse
    
    def _update_func(self, f_vec, K, K_sparse, R_sparse):
        """
        Update function
        
        Parameters
        ----------
        f_vec : array_like shape=(n_classes * n_samples, 1)
            latent function values
        
        K : Kernel matrix in dense form, shape=(n_classes * n_samples, n_samples)
        K_sparse : Kernel matrix in sparse form, shape=(n_classes * n_samples, n_samples)
        
        R_sparse: Stacked identity matrix, in sparse form
        """
        
        n_classes = self.n_classes
        n_samples = f_vec.size / n_classes
        
        pi_vec = self.soft_max(f_vec)
        pi_mat = pi_vec.reshape(n_classes, n_samples)
        
        # Stacked diagonal matrices: p of size Cn x n
        P = np.zeros((n_classes * n_samples, n_samples)) 
        
        # Stack vertically to build p 
        for c in range(n_classes):
            P[c * n_samples : (c + 1) * n_samples,:] = np.diagflat(pi_mat[c,:])
        
        # W = D - p x p^T
        W = np.diagflat(pi_vec) - np.dot(P, P.T)
        y_pi_diff = self.y - pi_vec
        b = np.dot(W, f_vec) + y_pi_diff

        # E is block diagonal
        E = []
        
        F = np.zeros((n_samples, n_samples))
        # B = I + W^1/2 x K x W^1/2
        # Compute B^-1 and det(B)
        B_log_det = 0
        for c in range(n_classes):
            D_c_sqrt = np.diagflat(pi_mat[c,:] ** .5)
            B_c = np.eye(n_samples) + np.dot(np.dot(D_c_sqrt, K[c]), D_c_sqrt)
            # B_c is stable
            L = linalg.cholesky(B_c, lower=True)
            L_inv = inv_triangular(L, lower = True)
            # det(L x L^T) = det(L)^2
            B_log_det = B_log_det + 2 * np.sum(np.log(np.diagonal(L)))
            # B_c = L x L^T => B_c^-1 = (L^T)^-1 x L^-1 = (L^-1)^T x L^-1 
            B_c_inv = np.dot(L_inv.T, L_inv)
            E_c = np.dot(np.dot(D_c_sqrt, B_c_inv), D_c_sqrt)
            E.append(E_c)
            # F = R^T x E x R = sum{E_c}
            F = F + E_c
        
        L = linalg.cholesky(F, lower=True)
        L_inv = inv_triangular(L, lower = True)
        F_inv = np.dot(L_inv.T, L_inv)

        E_sparse = sparse.block_diag(E, format='csr')
        c = E_sparse.dot(K_sparse).dot(b)
        ERM = E_sparse.dot(R_sparse).dot(F_inv)
        RTC = R_sparse.T.dot(c)
        # a = K^-1 x f
        a = b - c + np.dot(ERM, RTC)
        f_vec = K_sparse.dot(a)
        # gradiant = -K^-1 x f + y - - 
        gradiant = -a + y_pi_diff
        # At maximum must have f = K x (y - -) and gradiant = 0
        
        return f_vec, E, F_inv, B_log_det, a, pi_vec, gradiant
        
    def reduced_likelihood_function(self, theta=None, f_vec=None):
        """
        This function determines the laplace approximated marginal
        likelihood function for the given autocorrelation parameters theta.

        Maximizing this function wrt the latent function values is
        equivalent to maximizing the likelihood of the assumed joint Gaussian
        distribution of the observations y evaluated onto the design of
        experiments X.

        Parameters
        ----------
        theta : array_like, optional
            An array containing the autocorrelation parameters at which the
            Gaussian Process model parameters should be determined.
            Default uses the built-in autocorrelation parameters
            (ie ``theta = self.theta_``).
        
        f_vec: array_like, optional
            Starting values for latent function values.
            Default uses vector of zeros.

        Returns
        -------
        reduced_likelihood_function_value : double
            The value of the reduced likelihood function associated to the
            given latent function values f_vec.
            
        par : dict
            A dictionary containing the requested Gaussian Process model
            parameters:

                f_vec
                        Latent function values
                E
                        Computed block-diagonal matrix to find hessian inverse.
                F_inv
                        Internal square matrix to find hessian inverse.
                pi_vec
                        Squashed decision function values.
        """

        n_classes = self.n_classes
        n_samples = self.X.shape[0]

        if theta is None:
            # Use built-in autocorrelation parameters
            theta = self.theta_
        
        if f_vec is None:
            # Column vector of latent function values at all training 
            #     points for all classes, shape = (Cn,1)
            f_vec = np.zeros((n_classes * n_samples, 1))

        D, ij = self._cache['D'], self._cache['ij']
        
        # Numerical considerations:
        # should avoid computing K^-1, formulating based on fmin_ncg needs this term  
        #       see if we can change the cost function, with new hessian and prime
        #       in a way that K^-1 is avoided (hopefully replaced by B^-1)
        #       then we could use fmin_ncg
        
        # Kernel is block diagonal with equivalant dense shape of Cn x Cn
        K = []
        for c in range(n_classes):
            # Assumption is that C latent processes are uncorrelated (between classes)
            # Calculate matrix of distances D between samples for each class
            
            # Correlation model
            r = self.corr(theta[c,:], D)
            # Set up kernel matrix
            R = np.eye(n_samples) * (1. + self.nugget)
            R[ij[:, 0], ij[:, 1]] = r
            R[ij[:, 1], ij[:, 0]] = r
            
            K.append(R)

        K_sparse = sparse.block_diag(K, format='csr')

        # Get R from cache
        R_sparse = self._get_R()

        # Initialize output
        reduced_likelihood_function_value = - np.inf

        k = 0
        xtol = len(f_vec) * self.tol
        update = [2 * xtol]
        
        f_vec_old = f_vec
        objective = 0
        while k < self.max_iter and np.sum(np.abs(update)) > xtol:
            f_vec, E, F_inv, B_log_det, a, pi_vec, gradiant = self._update_func(f_vec, K, K_sparse, R_sparse)
            # TODO: use gradiant to update step-size and feed it back to _update_func
            
            # update to x is inherent to the algorithm above 
            update = f_vec - f_vec_old
            f_vec_old = f_vec
            
            objective = -0.5 * np.dot(a.T, f_vec) + np.dot(self.y.T, f_vec)
            objective = np.asarray(objective).squeeze()
            
            f_mat = f_vec.reshape(n_classes, n_samples)
            logsquash_term = np.sum(logsumexp(f_mat, axis = 0))
            objective = objective - logsquash_term
            
            k = k + 1
            
        if self.verbose:
            print "Maximum objective reached %s" % objective
            print "Number of iterations %s (of total %s)" % (k, self.max_iter)
        
        # Reduced likelyhood value 
        reduced_likelihood_function_value = objective - 0.5 * B_log_det    
        
        # Save those that are used for prediction
        #  soome can be computed from f_vec for a lightweight edition
        par = {}
        par['f_vec'] = f_vec
        par['E'] = E
        par['F_inv'] = F_inv
        par['pi_vec'] = pi_vec
        
        return reduced_likelihood_function_value, par


    def _check_params(self, n_classes = None, n_samples = None):
        '''
        Check to make sure parameters are valid
        '''
        
        if self.method != 'laplace':
            raise ValueError("Method can only be 'laplace' currently.")
        
        if (n_classes is not None
                and n_classes < 2):
            raise ValueError("at least two classes are needed for training.")
        
        theta0 = array2d(self.theta0)
        if n_classes is not None:
            # if for single class given, extend it to all classes
            if theta0.shape[0] == 1: 
                theta0 = np.tile(theta0[0,:], (n_classes,1))
            elif theta0.shape[0] != n_classes:
                raise ValueError("first dimension of theta0 (if non-single) must be number of classs (= %s)." % n_classes)
        self.theta0 = theta0
        
        # Check nugget value
        self.nugget = np.asarray(self.nugget)
        if np.any(self.nugget) < 0.:
            raise ValueError("nugget must be positive or zero.")
        
        if (n_samples is not None
                and self.nugget.shape not in [(), (n_samples,)]):
            raise ValueError("nugget must be either a scalar "
                             "or array of length n_samples.")
            
        if (n_samples is not None
                and n_classes is not None):
                if self.max_iter < 1:
                    self.max_iter = 10 * n_samples * n_classes
        
        # Check correlation model
        if not callable(self.corr):
            if self.corr in self._correlation_types:
                self.corr = self._correlation_types[self.corr]
            else:
                raise ValueError("corr should be one of %s or callable, "
                                 "%s was given."
                                 % (self._correlation_types.keys(), self.corr))

