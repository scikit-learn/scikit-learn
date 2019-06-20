import numpy

from ..base import BaseEstimator, RegressorMixin
from ..utils import check_X_y

# Author: Ben Hammel <bdhammel@gmail.com>
#         Nick Sullivan-Molina


class LsqEllipse(BaseEstimator, RegressorMixin):
    """Lest Squares fitting algorithm for Elliptical data

    References
    ----------
    (*) Halir, R., Flusser, J.: 'Numerically Stable Direct Least Squares Fitting of
    Ellipses'
    (**) http://mathworld.wolfram.com/Ellipse.html

    Attributes
    ----------
    coef_ : array, shape (n_features, ) or (n_targets, n_features)
        Estimated coefficients for the linear regression problem.
        If multiple targets are passed during the fit (y 2D), this
        is a 2D array of shape (n_targets, n_features), while if only
        one target is passed, this is a 1D array of length n_features.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LsqEllipse
    >>> X = np.array([])
    >>> y = np.array([])
    >>> reg = LsqEllipse().fit(X, y)
    # >>> reg.score(X, y)
    # 1.0
    >>> reg.get_parameters()
    """

    def __init__(self):
        self._coef = None

    def fit(self, X, y):
        """Fit the data

        Parameters
        ----------
        X : numpy.ndarray
            X data values for the x-y data pairs to fit
        y : numpu.ndarry
            y data values for the x-y data pairs to fit
        """
        X, y = check_X_y(
            X, y, y_numeric=True, ensure_2d=False, ensure_min_samples=5
        )

        # Quadratic part of design matrix [eqn. 15] from (*)
        D1 = numpy.mat(numpy.vstack([X**2, X * y, y**2])).T
        # Linear part of design matrix [eqn. 16] from (*)
        D2 = numpy.mat(numpy.vstack([X, y, numpy.ones(len(X))])).T

        # forming scatter matrix [eqn. 17] from (*)
        S1 = D1.T * D1
        S2 = D1.T * D2
        S3 = D2.T * D2

        # Constraint matrix [eqn. 18]
        C1 = numpy.mat('0. 0. 2.; 0. -1. 0.; 2. 0. 0.')

        # Reduced scatter matrix [eqn. 29]
        M = C1.I * (S1 - S2 * S3.I * S2.T)

        # M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors from this
        # equation [eqn. 28]
        eigval, eigvec = numpy.linalg.eig(M)

        # eigenvector must meet constraint 4ac - b^2 to be valid.
        cond = 4 * numpy.multiply(eigvec[0, :], eigvec[2, :]) - numpy.power(eigvec[1, :], 2)
        a1 = eigvec[:, numpy.nonzero(cond.A > 0)[1]]

        # |d f g> = -S3^(-1) * S2^(T)*|a b c> [eqn. 24]
        a2 = -S3.I * S2.T * a1

        # eigenvectors |a b c d f g>
        # list of the coefficients describing an ellipse [a,b,c,d,f,g]
        # corresponding to ax**2+2bxy+cy**2+2dx+2fy+g
        self._coef = numpy.vstack([a1, a2])

    def get_parameters(self):
        """finds the important parameters of the fitted ellipse

        Theory taken form http://mathworld.wolfram

        Args
        -----
        coef (list): list of the coefficients describing an ellipse
           [a,b,c,d,f,g] corresponding to ax**2+2bxy+cy**2+2dx+2fy+g

        Returns
        _______
        center (List): of the form [x0, y0]
        width (float): major axis
        height (float): minor axis
        phi (float): rotation of major axis form the x-axis in radians
        """

        # eigenvectors are the coefficients of an ellipse in general form
        # a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0
        # [eqn. 15) from (**) or (***)
        a = self._coef[0, 0]
        b = self._coef[1, 0] / 2.
        c = self._coef[2, 0]
        d = self._coef[3, 0] / 2.
        f = self._coef[4, 0] / 2.
        g = self._coef[5, 0]

        # finding center of ellipse [eqn.19 and 20] from (**)
        x0 = (c * d - b * f) / (b**2. - a * c)
        y0 = (a * f - b * d) / (b**2. - a * c)

        # Find the semi-axes lengths [eqn. 21 and 22] from (**)
        numerator = 2 * (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g)
        denominator1 = (b * b - a * c) * ((c - a) * numpy.sqrt(1 + 4 * b * b / ((a - c) * (a - c))) - (c + a))
        denominator2 = (b * b - a * c) * ((a - c) * numpy.sqrt(1 + 4 * b * b / ((a - c) * (a - c))) - (c + a))
        width = numpy.sqrt(numerator / denominator1)
        height = numpy.sqrt(numerator / denominator2)

        # angle of counterclockwise rotation of major-axis of ellipse to x-axis
        # [eqn. 23] from (**) or [eqn. 26] from (***).
        phi = .5 * numpy.arctan((2. * b) / (a - c))

        self._center = [x0, y0]
        self._width = width
        self._height = height
        # angle of counterclockwise rotation of major-axis of ellipse to
        # x-axis [eqn. 23] from (**)
        self._phi = phi

        return self._center, self._width, self._height, self._phi
