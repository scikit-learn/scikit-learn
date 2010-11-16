#!/usr/bin/python
# -*- coding: utf-8 -*-

################
# Dependencies #
################

import numpy as np
from scipy import linalg, optimize, rand
from ..base import BaseEstimator
from .regression import regpoly0
from .correlation import correxp2, corriid
machine_epsilon = np.finfo(np.double).eps
if hasattr(linalg, 'solve_triangular'):
    # only in scipy since 0.9
    solve_triangular = linalg.solve_triangular
else:
    # slower, but works
    solve_triangular = linalg.solve


##############################
# The Gaussian Process class #
##############################


class GaussianProcess(BaseEstimator):
    """
    A class that implements scalar Gaussian Process based prediction (also
    known as Kriging).

    Example
    -------
    import numpy as np
    from scikits.learn.gaussian_process import GaussianProcess
    import pylab as pl

    f = lambda x: x*np.sin(x)
    X = np.array([1., 3., 5., 6., 7., 8.])
    Y = f(X)
    gp = GaussianProcess(regr=regpoly0, corr=correxp2, theta0=1e-1, \
                         thetaL=1e-3, thetaU=1e0, random_start=100)
    gp.fit(X, Y)

    pl.figure(1)
    x = np.linspace(0,10,1000)
    y, MSE = gp.predict(x, eval_MSE=True)
    sigma = np.sqrt(MSE)
    pl.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
    pl.plot(X, Y, 'r.', markersize=10, label=u'Observations')
    pl.plot(x, y, 'k-', label=u'$\widehat{f}(x) = {\\rm BLUP}(x)$')
    pl.fill(np.concatenate([x, x[::-1]]), \
            np.concatenate([y - 1.9600 * sigma, (y + 1.9600 * sigma)[::-1]]), \
            alpha=.5, fc='b', ec='None', label=u'95\% confidence interval')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.legend(loc='upper left')

    pl.figure(2)
    theta_values = np.logspace(np.log10(gp.thetaL), np.log10(gp.thetaU),100)
    psi_values = []
    for t in theta_values:
            psi_values.append(gp.reduced_likelihood_function(theta=t)[0])
    pl.plot(theta_values, psi_values)
    pl.xlabel(u'$\\theta$')
    pl.ylabel(u'Score')
    pl.xscale('log')

    pl.show()

    Methods
    -------
    fit(X, y) : self
        Fit the model.

    predict(X) : array
        Predict using the model.

    Todo
    ----
    o Add the 'sparse' storage mode for which the correlation matrix is stored
      in its sparse eigen decomposition format instead of its full Cholesky
      decomposition.

    Implementation details
    ----------------------
    The presentation implementation is based on a translation of the DACE
    Matlab toolbox.

    See references:
    [1] H.B. Nielsen, S.N. Lophaven, H. B. Nielsen and J. Sondergaard (2002).
        DACE - A MATLAB Kriging Toolbox.
        http://www2.imm.dtu.dk/~hbn/dace/dace.pdf
    """

    def __init__(self, regr=regpoly0, corr=correxp2, beta0=None, \
                 storage_mode='full', verbose=True, theta0=1e-1, \
                 thetaL=None, thetaU=None, optimizer='fmin_cobyla', \
                 random_start=1, normalize=True, \
                 nugget=10. * machine_epsilon):
        """
        The Gaussian Process model constructor.

        Parameters
        ----------
        regr : lambda function, optional
            A regression function returning an array of outputs of the linear
            regression functional basis. The number of observations n_samples
            should be greater than the size p of this basis.
            Default assumes a simple constant regression trend (see regpoly0).

        corr : lambda function, optional
            A stationary autocorrelation function returning the autocorrelation
            between two points x and x'.
            Default assumes a squared-exponential autocorrelation model (see
            correxp2).

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
            Default is verbose = True.

        theta0 : double array_like, optional
            An array with shape (n_features, ) or (1, ).
            The parameters in the autocorrelation model.
            If thetaL and thetaU are also specified, theta0 is considered as
            the starting point for the maximum likelihood rstimation of the
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
            Design sites X and observations y are centered and reduced wrt
            means and standard deviations estimated from the n_samples
            observations provided.
            Default is normalize = True so that data is normalized to ease
            maximum likelihood estimation.

        nugget : double, optional
            Introduce a nugget effect to allow smooth predictions from noisy
            data.
            Default assumes a nugget close to machine precision for the sake of
            robustness (nugget = 10.*machine_epsilon).

        optimizer : string, optional
            A string specifying the optimization algorithm to be used
            ('fmin_cobyla' is the sole algorithm implemented yet).
            Default uses 'fmin_cobyla' algorithm from scipy.optimize.

        random_start : int, optional
            The number of times the Maximum Likelihood Estimation should be
            performed from a random starting point.
            The first MLE always uses the specified starting point (theta0),
            the next starting points are picked at random according to an
            exponential distribution (log-uniform on [thetaL, thetaU]).
            Default does not use random starting point (random_start = 1).

        Returns
        -------
        gp : self
            A Gaussian Process model object awaiting data to be fitted to.
        """

        self.regr = regr
        self.corr = corr
        self.beta0 = beta0

        # Check storage mode
        if storage_mode != 'full' and storage_mode != 'light':
            if storage_mode == 'sparse':
                raise NotImplementedError("The 'sparse' storage mode is not " \
                      + "supported yet. Please contribute!")
            else:
                raise ValueError("Storage mode should either be 'full' or " \
                      + "'light'. Unknown storage mode: " + str(storage_mode))
        else:
            self.storage_mode = storage_mode

        self.verbose = verbose

        # Check correlation parameters
        self.theta0 = np.atleast_2d(np.asanyarray(theta0, dtype=np.float))
        self.thetaL, self.thetaU = thetaL, thetaU
        lth = self.theta0.size

        if self.thetaL is not None and self.thetaU is not None:
            # Parameters optimization case
            self.thetaL = np.atleast_2d(np.asanyarray(thetaL, dtype=np.float))
            self.thetaU = np.atleast_2d(np.asanyarray(thetaU, dtype=np.float))

            if self.thetaL.size != lth or self.thetaU.size != lth:
                raise ValueError("theta0, thetaL and thetaU must have the " \
                      + "same length")
            if np.any(self.thetaL <= 0) or np.any(self.thetaU < self.thetaL):
                raise ValueError("The bounds must satisfy O < thetaL <= " \
                      + "thetaU")

        elif self.thetaL is None and self.thetaU is None:
            # Given parameters case
            if np.any(self.theta0 <= 0):
                raise ValueError("theta0 must be strictly positive")

        elif self.thetaL is None or self.thetaU is None:
            # Error
            raise Exception("thetaL and thetaU should either be both or " \
                  + "neither specified")

        # Store other parameters
        self.normalize = normalize
        self.nugget = nugget
        self.optimizer = optimizer
        self.random_start = int(random_start)

    def fit(self, X, y):
        """
        The Gaussian Process model fitting method.

        Parameters
        ----------
        X : double array_like
            An array with shape (n_samples, n_features) with the design sites
            at which observations were made.

        y : double array_like
            An array with shape (n_features, ) with the observations of the
            scalar output to be predicted.

        Returns
        -------
        gp : self
            A fitted Gaussian Process model object awaiting data to perform
            predictions.
        """

        # Force data to numpy.array type (from coding guidelines)
        X = np.asanyarray(X, dtype=np.float)
        y = np.asanyarray(y, dtype=np.float)

        # Check design sites & observations
        n_samples_X = X.shape[0]
        if X.ndim > 1:
            n_features = X.shape[1]
        else:
            n_features = 1
        X = X.reshape(n_samples_X, n_features)

        n_samples_y = y.shape[0]
        if y.ndim > 1:
            raise NotImplementedError("y has more than one dimension. This " \
                  + "is not supported yet (scalar output prediction only). " \
                  + "Please contribute!")
        y = y.reshape(n_samples_y, 1)

        if n_samples_X != n_samples_y:
            raise Exception("X and y must have the same number of rows!")
        else:
            n_samples = n_samples_X

        # Normalize data or don't
        if self.normalize:
            mean_X = np.mean(X, axis=0)
            std_X = np.std(X, axis=0)
            mean_y = np.mean(y, axis=0)
            std_y = np.std(y, axis=0)
            std_X[std_X == 0.] = 1.
            std_y[std_y == 0.] = 1.
        else:
            mean_X = np.array([0.])
            std_X = np.array([1.])
            mean_y = np.array([0.])
            std_y = np.array([1.])

        X = (X - mean_X) / std_X
        y = (y - mean_y) / std_y

        # Calculate matrix of distances D between samples
        mzmax = n_samples * (n_samples - 1) / 2
        ij = np.zeros([mzmax, 2])
        D = np.zeros([mzmax, n_features])
        ll = np.array([-1])
        for k in range(n_samples-1):
            ll = ll[-1] + 1 + range(n_samples - k - 1)
            ij[ll] = np.concatenate([[np.repeat(k, n_samples - k - 1, 0)], \
                                     [np.arange(k + 1, n_samples).T]]).T
            D[ll] = X[k] - X[(k + 1):n_samples]
        if np.min(np.sum(np.abs(D), 1)) == 0. and self.corr != corriid:
            raise Exception("Multiple X are not allowed")

        # Regression matrix and parameters
        F = self.regr(X)
        n_samples_F = F.shape[0]
        if F.ndim > 1:
            p = F.shape[1]
        else:
            p = 1
        if n_samples_F != n_samples:
            raise Exception("Number of rows in F and X do not match. Most " \
                          + "likely something is going wrong with the " \
                          + "regression model.")
        if p > n_samples_F:
            raise Exception("Ordinary least squares problem is undetermined " \
                           + "n_samples=%d must be greater than the " \
                           + "regression model size p=%d!" % (n_samples, p))
        if self.beta0 is not None \
           and (self.beta0.shape[0] != p or self.beta0.ndim > 1):
            raise Exception("Shapes of beta0 and F do not match.")

        # Set attributes
        self.X = X
        self.y = y
        self.D = D
        self.ij = ij
        self.F = F
        self.X_sc = np.concatenate([[mean_X], [std_X]])
        self.y_sc = np.concatenate([[mean_y], [std_y]])

        # Determine Gaussian Process model parameters
        if self.thetaL is not None and self.thetaU is not None:
            # Maximum Likelihood Estimation of the parameters
            if self.verbose:
                print "Performing Maximum Likelihood Estimation of the " \
                    + "autocorrelation parameters..."
            self.theta, self.reduced_likelihood_function_value, self.par = \
                self.arg_max_reduced_likelihood_function()
            if np.isinf(self.reduced_likelihood_function_value):
                raise Exception("Bad parameter region. " \
                              + "Try increasing upper bound")

        else:
            # Given parameters
            if self.verbose:
                print "Given autocorrelation parameters. " \
                    + "Computing Gaussian Process model parameters..."
            self.theta = self.theta0
            self.reduced_likelihood_function_value, self.par = \
                self.reduced_likelihood_function()
            if np.isinf(self.reduced_likelihood_function_value):
                raise Exception("Bad point. Try increasing theta0")

        if self.storage_mode == 'light':
            # Delete heavy data (it will be computed again if required)
            # (it is required only when MSE is wanted in self.predict)
            if self.verbose:
                print "Light storage mode specified. " \
                    + "Flushing autocorrelation matrix..."
            self.D = None
            self.ij = None
            self.F = None
            self.par['C'] = None
            self.par['Ft'] = None
            self.par['G'] = None

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

        verbose : boolean, optional
            A boolean specifying the verbose level.
            Default is verbose = True.

        batch_size : integer, optional
            An integer giving the maximum number of points that can be
            evaluated simulatneously (depending on the available memory).
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

        # Check
        if np.any(np.isnan(self.par['beta'])):
            raise Exception("Not a valid GaussianProcess. " \
                          + "Try fitting it again with different parameters " \
                          + "theta.")

        # Check design & trial sites
        X = np.asanyarray(X, dtype=np.float)
        n_samples = self.X.shape[0]
        if self.X.ndim > 1:
            n_features = self.X.shape[1]
        else:
            n_features = 1
        n_eval = X.shape[0]
        if X.ndim > 1:
            n_features_X = X.shape[1]
        else:
            n_features_X = 1
        X = X.reshape(n_eval, n_features_X)

        if n_features_X != n_features:
            raise ValueError("The number of features in X (X.shape[1] = %d) " \
                           + "should match the sample size used for fit() " \
                           + "which is %d." % (n_features_X, n_features))

        if batch_size is None:
            # No memory management
            # (evaluates all given points in a single batch run)

            # Normalize trial sites
            X = (X - self.X_sc[0][:]) / self.X_sc[1][:]

            # Initialize output
            y = np.zeros(n_eval)
            if eval_MSE:
                MSE = np.zeros(n_eval)

            # Get distances to design sites
            dx = np.zeros([n_eval * n_samples, n_features])
            kk = np.arange(n_samples).astype(int)
            for k in range(n_eval):
                dx[kk] = X[k] - self.X
                kk = kk + n_samples

            # Get regression function and correlation
            f = self.regr(X)
            r = self.corr(self.theta, dx).reshape(n_eval, n_samples)

            # Scaled predictor
            y_ = np.dot(f, self.par['beta']) \
               + np.dot(r, self.par['gamma'])

            # Predictor
            y = (self.y_sc[0] + self.y_sc[1] * y_).ravel()

            # Mean Squared Error
            if eval_MSE:
                par = self.par
                if par['C'] is None:
                    # Light storage mode (need to recompute C, F, Ft and G)
                    if self.verbose:
                        print "This GaussianProcess used light storage mode " \
                            + "at instanciation. Need to recompute " \
                            + "autocorrelation matrix..."
                    reduced_likelihood_function_value, par = \
                        self.reduced_likelihood_function()

                rt = solve_triangular(par['C'], r.T)
                if self.beta0 is None:
                    # Universal Kriging
                    u = solve_triangular(self.par['G'].T, \
                                         np.dot(self.par['Ft'].T, rt) - f.T)
                else:
                    # Ordinary Kriging
                    u = np.zeros(y.shape)

                MSE = self.par['sigma2'] * (1. - (rt ** 2.).sum(axis=0) \
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

                y, MSE = np.array([]), np.array([])
                for k in range(n_eval / batch_size):
                    batch_from = k * batch_size
                    batch_to = min([(k + 1) * batch_size + 1, n_eval + 1])
                    X_batch = X[batch_from:batch_to][:]
                    y_batch, MSE_batch = \
                        self.predict(X_batch, eval_MSE=eval_MSE, \
                                     batch_size=None)
                    y.append(y_batch)
                    MSE.append(MSE_batch)

                return y, MSE

            else:

                y = np.array([])
                for k in range(n_eval / batch_size):
                    batch_from = k * batch_size
                    batch_to = min([(k + 1) * batch_size + 1, n_eval + 1])
                    X_batch = X[batch_from:batch_to][:]
                    y_batch = \
                        self.predict(X_batch, eval_MSE=eval_MSE, \
                                     batch_size=None)

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
            (ie theta = self.theta).

        Returns
        -------
        reduced_likelihood_function_value : double
            The value of the reduced likelihood function associated to the
            given autocorrelation parameters theta.

        par : dict
            A dictionary containing the requested Gaussian Process model
            parameters:

            par['sigma2'] : Gaussian Process variance.
            par['beta'] : Generalized least-squares regression weights for
                          Universal Kriging or given beta0 for Ordinary
                          Kriging.
            par['gamma'] : Gaussian Process weights.
            par['C'] : Cholesky decomposition of the correlation matrix [R].
            par['Ft'] : Solution of the linear equation system : [R] x Ft = F
            par['G'] : QR decomposition of the matrix Ft.
            par['detR'] : Determinant of the correlation matrix raised at power
                          1/n_samples.
        """

        if theta is None:
            # Use built-in autocorrelation parameters
            theta = self.theta

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
            if self.X.ndim > 1:
                n_features = self.X.shape[1]
            else:
                n_features = 1
            mzmax = n_samples * (n_samples - 1) / 2
            ij = np.zeros([mzmax, n_features])
            D = np.zeros([mzmax, n_features])
            ll = np.array([-1])
            for k in range(n_samples-1):
                ll = ll[-1] + 1 + range(n_samples - k - 1)
                ij[ll] = \
                    np.concatenate([[np.repeat(k, n_samples - k - 1, 0)], \
                                    [np.arange(k + 1, n_samples).T]]).T
                D[ll] = self.X[k] - self.X[(k + 1):n_samples]
            if min(sum(abs(D), 1)) == 0. and self.corr != corriid:
                raise Exception("Multiple X are not allowed")
            F = self.regr(self.X)

        # Set up R
        r = self.corr(theta, D)
        R = np.eye(n_samples) * (1. + self.nugget)
        R[ij.astype(int)[:, 0], ij.astype(int)[:, 1]] = r
        R[ij.astype(int)[:, 1], ij.astype(int)[:, 0]] = r

        # Cholesky decomposition of R
        try:
            C = linalg.cholesky(R, lower=True)
        except linalg.LinAlgError:
            return reduced_likelihood_function_value, par

        # Get generalized least squares solution
        Ft = solve_triangular(C, F)
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
                raise Exception("F is too ill conditioned. Poor combination " \
                              + "of regression model and observations.")
            else:
                # Ft is too ill conditioned, get out (try different theta)
                return reduced_likelihood_function_value, par

        Yt = solve_triangular(C, self.y)
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
        par['sigma2'] = sigma2 * self.y_sc[1] ** 2.
        par['beta'] = beta
        par['gamma'] = solve_triangular(C.T, rho)
        par['C'] = C
        par['Ft'] = Ft
        par['G'] = G

        return reduced_likelihood_function_value, par

    def arg_max_reduced_likelihood_function(self):
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
            print "The chosen optimizer is: " + str(self.optimizer)
            if self.random_start > 1:
                print str(self.random_start) + " random starts are required."

        percent_completed = 0.

        if self.optimizer == 'fmin_cobyla':

            minus_reduced_likelihood_function = lambda log10t: \
                - self.reduced_likelihood_function(theta=10. ** log10t)[0]

            constraints = []
            for i in range(self.theta0.size):
                constraints.append(lambda log10t: \
                            log10t[i] - np.log10(self.thetaL[0, i]))
                constraints.append(lambda log10t: \
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
                log10_optimal_theta = \
                    optimize.fmin_cobyla(minus_reduced_likelihood_function, \
                                    np.log10(theta0), constraints, iprint=0)

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
                        print str(5 * percent_completed) + "% completed"

        else:

            raise NotImplementedError("This optimizer ('%s') is not " \
                                    + "implemented yet. Please contribute!" \
                                    % self.optimizer)

        return best_optimal_theta, best_optimal_rlf_value, best_optimal_par

    def score(self, X_test, y_test):
        """
        This score function returns the deviations of the Gaussian Process
        model evaluated onto a test dataset.

        Parameters
        ----------
        X_test : array_like
            The feature test dataset with shape (n_tests, n_features).

        y_test : array_like
            The target test dataset (n_tests, ).

        Returns
        -------
        score_values : array_like
            The deviations between the prediction and the targets:
            y_pred - y_test.
        """

        return np.ravel(self.predict(X_test, eval_MSE=False)) - y_test
