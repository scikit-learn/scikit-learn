"""L1 trend filtering implementation"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import logging
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from sklearn.base import BaseEstimator, TransformerMixin


class L1TrendFiltering(BaseEstimator, TransformerMixin):
    """
    Represents a class implementing L1 trend filtering using iterative optimization.

    This class aims to fit and transform input data by applying L1 trend filtering,
    designed to reconstruct signals by reducing noise while preserving important
    trends. The optimization process involves solving convex problems iteratively
    using primal-dual algorithms and backtracking line searches.

    Attributes
    ----------
    _lambda : float
        Regularization parameter controlling the smoothness of the filtered signal.
    _alpha : float
        Initial step size for optimization algorithms like gradient descent.
    _beta : float
        Line search backtracking parameter defining step size reductions.
    _mu : float
        Scaling parameter used in the optimization process.
    _max_iter : int
        Maximum iterations allowed for the optimization procedure.
    _max_ls_iter : int
        Maximum iterations allowed for backtracking line search.
    _tol : float
        Convergence tolerance to determine stopping criteria.
    """
    _logger = logging.getLogger(__name__)

    def __init__(self, lambda_=1.0, alpha=0.01, beta=0.5, mu=2.0, max_iter=40, max_ls_iter=20, tol=1e-4):
        """
        Represents an optimization configuration class for setting various hyperparameters
        used in optimization routines. This class allows the user to control parameters
        such as regularization strength, step size, and tolerance, among others, for
        fine-tuning optimization behavior.

        Attributes
        ----------
        lambda_ : float
            Regularization parameter that controls the strength of regularization.
        alpha : float
            Initial step size for optimization algorithms such as gradient descent.
        beta : float
            Backtracking line search parameter that governs the step size reduction.
        mu : float
            Parameter that controls the update rule scaling or other problem-dependent factors.
        max_iter : int
            Maximum number of iterations allowed for the optimization process.
        max_ls_iter : int
            Maximum number of iterations allowed for line search in an optimization routine.
        tol : float
            Tolerance criteria for convergence. Optimization stops when this value is achieved.
        """
        self._lambda = lambda_
        self._alpha = alpha
        self._beta = beta
        self._mu = mu
        self._max_iter = max_iter
        self._max_ls_iter = max_ls_iter
        self._tol = tol

    @classmethod
    def _second_difference_matrix(cls, n):
        """
        Constructs a second-order finite difference matrix.

        This method generates a sparse second-order difference matrix, which is used
        to compute the second numeric derivatives in finite differences. The resulting
        matrix is created using diagonals with specific values and offsets, providing
        an efficient representation of the second derivative operator in a
        discretized form.

        Parameters
        ----------
        n : int
            The size of the finite grid or the number of discrete points.

        Returns
        -------
        numpy.ndarray
            A 2D array where each row represents the second derivative finite
            difference operator applied to the grid. The array's shape is (n - 2, n).
        """
        diagonals = [np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)]
        offsets = [0, 1, 2]
        return diags(diagonals, offsets, shape=(n - 2, n)).toarray()

    def _objective(self, D, y, z, mu1, mu2, t):
        """
        Computes the primary objective, dual objective, and duality gap for an optimization
        problem using input data and specified parameters. The function solves the objective
        based on a linear system and uses various mathematical operations to determine the
        results.

        Parameters
        ----------
        D : numpy.ndarray
            The input matrix for the optimization problem, used in matrix calculations.
        y : numpy.ndarray
            The input vector, generally representing the target or observed data.
        z : numpy.ndarray
            The dual variable vector used in duality and optimization computations.
        mu1 : numpy.ndarray
            The first Lagrange multiplier vector for regularization terms.
        mu2 : numpy.ndarray
            The second Lagrange multiplier vector for regularization terms.
        t : float
            The scalar parameter controlling precision or regularization.

        Returns
        -------
        pobj : float
            The primal objective value, representing the minimization objective.
        dobj : float
            The dual objective value, representing the alternate maximization objective.
        gap : float
            The duality gap, representing the difference between primal and dual objectives.
        """
        DTz = D.T @ z
        DDTz = D @ DTz
        Dy = D @ y

        w = Dy - (mu1 - mu2)
        pobj1 = 0.5 * np.dot(w, np.linalg.solve(D @ D.T, w)) + self._lambda * np.sum(mu1 + mu2)
        pobj2 = 0.5 * np.dot(DTz, DTz) + self._lambda * np.sum(np.abs(Dy - DDTz))
        pobj = min(pobj1, pobj2)
        dobj = -0.5 * np.dot(DTz, DTz) + np.dot(Dy, z)
        gap = pobj - dobj
        return pobj, dobj, gap

    def fit(self, X, y=None):
        """
        Fits the model according to the given training data.

        This method trains the model using the input features `X` and, optionally,
        the target labels `y`. The model learns patterns from the data to make
        accurate predictions or decisions. The input data should have
        appropriate dimensions and align with the training requirements of the model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input training data where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like of shape (n_samples,), default=None
            The target values corresponding to the training data. If None, the method
            assumes an unsupervised learning scenario.

        Returns
        -------
        self : object
            Returns the instance of the model after fitting. This allows
            method chaining and further operations on the fitted model.
        """
        return self

    def transform(self, X):
        """
        Perform the L1 trend filtering optimization.
        :param X: Input signal (n_samples x n_features).
        :return: Filtered signal with the same shape as X.
        """
        if len(X.shape) == 1:
            X = X.reshape(-1, 1)

        results = []
        for col in X.T:
            results.append(self._l1tf_1d(col))

        return np.array(results).T

    def _l1tf_1d(self, y):
        """
        Computes the solution to a one-dimensional total variation denoising problem using
        an iterative primal-dual algorithm with line search. The method minimizes the
        objective function for L1 trend filtering, utilizing second difference matrices
        and optimization via Newton's method.

        Parameters
        ----------
        y : numpy.ndarray
            Input data for which the L1 trend filtering is applied. It is assumed to be
            a one-dimensional array.

        Returns
        -------
        numpy.ndarray
            The smoothed or denoised version of input `y`, obtained as the primal
            optimal point for the L1 trend filtering problem.

        Raises
        ------
        ValueError
            If parameters or interim computations result in invalid states or conditions
            inappropriate for solving the optimization problem.
        """
        n = len(y)
        m = n - 2
        D = self._second_difference_matrix(n)
        Dy = D @ y
        DDT = D @ D.T

        z = np.zeros(m)
        mu1 = np.ones(m)
        mu2 = np.ones(m)
        t = 1e-10

        for iters in range(self._max_iter):
            # Compute dual and residuals
            DTz = D.T @ z
            DDTz = D @ DTz
            w = Dy - (mu1 - mu2)

            pobj, dobj, gap = self._objective(D, y, z, mu1, mu2, t)
            if gap <= self._tol:
                return y - D.T @ z  # Primal optimal point

            if gap >= 0.2:
                t = max(2 * m * self._mu / gap, 1.2 * t)

            # Compute Newton step
            f1 = z - self._lambda
            f2 = -z - self._lambda
            S = DDT - diags(mu1 / f1 + mu2 / f2)
            r = -DDTz + Dy + (1 / t) / f1 - (1 / t) / f2
            dz = spsolve(S, r)
            dmu1 = -(mu1 + ((1 / t) + dz * mu1) / f1)
            dmu2 = -(mu2 + ((1 / t) - dz * mu2) / f2)

            # line search
            step = 1
            for _ in range(self._max_ls_iter):
                newz = z + step * dz
                newmu1 = mu1 + step * dmu1
                newmu2 = mu2 + step * dmu2
                newf1 = newz - self._lambda
                newf2 = -newz - self._lambda

                if np.all(newf1 < 0) and np.all(newf2 < 0):
                    z, mu1, mu2 = newz, newmu1, newmu2
                    break
                step *= self._beta

        if iters >= self._max_iter:
            message = f'maxiter {self._max_iter} exceeded'
            self._logger.warning(message)

        # last computed primal variable
        return y - D.T @ z
