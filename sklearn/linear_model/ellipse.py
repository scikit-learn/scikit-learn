import numpy as np

from ..base import BaseEstimator
from ..utils import check_X_y

# Author: Ben Hammel <bdhammel@gmail.com>
#         Nick Sullivan-Molina


class LsqEllipse(BaseEstimator):
    """
    Lest Squares fitting of Elliptical data

    References
    ----------
    (*) Halir R., Flusser J. 'Numerically Stable Direct Least Squares
    Fitting of Ellipses'
    (**) Weisstein, Eric W. "Ellipse." From MathWorld--A Wolfram Web Resource.
    http://mathworld.wolfram.com/Ellipse.html

    Attributes
    ----------
    _coef : array
        Estimated coefficients for the Least squares fit to the elliptical data
        containing the values [a,b,c,d,f,g].T corresponding to
        ax**2 + 2bxy + cy**2 + 2dx + 2fy + g

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LsqEllipse
    >>> X = np.array([ 1.,  0., -1., -0.,  1.])
    >>> y = np.array([ 0. ,  0.5,  0. , -0.5, -0. ])
    >>> reg = LsqEllipse().fit(X, y)
    >>> reg.as_parameters()
    """

    def __init__(self):
        self._coef = None

    def fit(self, X, y):
        """
        Fit the data

        Parameters
        ----------
        X : array, shape (n_points, )
            X data values for the x-y data pairs to fit
        y : array, shape (n_point, )
            y data values for the x-y data pairs to fit

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_X_y(
            X, y, y_numeric=True, ensure_2d=False, ensure_min_samples=5
        )

        # Quadratic part of design matrix [eqn. 15] from (*)
        D1 = np.mat(np.vstack([X**2, X * y, y**2])).T
        # Linear part of design matrix [eqn. 16] from (*)
        D2 = np.mat(np.vstack([X, y, np.ones_like(X)])).T

        # Forming scatter matrix [eqn. 17] from (*)
        S1 = D1.T * D1
        S2 = D1.T * D2
        S3 = D2.T * D2

        # Constraint matrix [eqn. 18]
        C1 = np.mat('0. 0. 2.; 0. -1. 0.; 2. 0. 0.')

        # Reduced scatter matrix [eqn. 29]
        M = C1.I * (S1 - S2 * S3.I * S2.T)

        # M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors from this
        # equation [eqn. 28]
        eigval, eigvec = np.linalg.eig(M)

        # Eigenvector must meet constraint 4ac - b^2 to be valid.
        cond = (
            4*np.multiply(eigvec[0, :], eigvec[2, :])
            - np.power(eigvec[1, :], 2)
        )
        a1 = eigvec[:, np.nonzero(cond.A > 0)[1]]

        # |d f g> = -S3^(-1) * S2^(T)*|a b c> [eqn. 24]
        a2 = -S3.I * S2.T * a1

        # Eigenvectors |a b c d f g>
        # list of the coefficients describing an ellipse [a,b,c,d,f,g]
        # corresponding to ax**2 + 2bxy + cy**2 + 2dx + 2fy + g
        self._coef = np.vstack([a1, a2])

        return self

    @property
    def coefficients(self):
        """
        List of the coefficients describing the fitted ellipse

        Returns
        -------
        [a,b,c,d,f,g] corresponding to ax**2 + 2bxy + cy**2 + 2dx + 2fy + g
        """
        return np.asarray(self._coef).ravel()

    def as_parameters(self):
        """Returns the definition of the fitted ellipse as localized parameters

        Returns
        _______
        center : list
            [x0, y0]
        width : float
            Semimajor axis
        height : float
            Semiminor axis
        phi : float
            The counterclockwise angle of rotation from the x-axis to the major
            axis of the ellipse
        """

        # Eigenvectors are the coefficients of an ellipse in general form
        # a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0
        # [eqn. 15) from (**) or (***)
        a = self.coefficients[0]
        b = self.coefficients[1] / 2.
        c = self.coefficients[2]
        d = self.coefficients[3] / 2.
        f = self.coefficients[4] / 2.
        g = self.coefficients[5]

        # Finding center of ellipse [eqn.19 and 20] from (**)
        x0 = (c*d - b*f) / (b**2. - a*c)
        y0 = (a*f - b*d) / (b**2. - a*c)
        center = [x0, y0]

        # Find the semi-axes lengths [eqn. 21 and 22] from (**)
        numerator = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
        denominator1 = (b * b - a * c) * (
            (c - a) * np.sqrt(1 + 4*b*b / ((a - c)*(a - c))) - (c + a)
        )
        denominator2 = (b*b - a*c) * (
            (a - c) * np.sqrt(1 + 4*b*b / ((a - c) * (a - c))) - (c + a)
        )
        width = np.sqrt(numerator / denominator1)
        height = np.sqrt(numerator / denominator2)

        # Angle of counterclockwise rotation of major-axis of ellipse to x-axis
        # [eqn. 23] from (**) or [eqn. 26] from (***).
        phi = .5 * np.arctan((2.*b) / (a - c))

        return center, width, height, phi

    def predict(self, n_points=None, t=None):
        """
        Return the X, Y values of the predicted ellipse
        Points are returned along the parametric curve of the ellipse as evenly
        spaced points starting at t=0 to t=2pi

        Parameters
        ---------
        n_points : int
            Number of points to return
        t : array
            Parametric points used to generate x-y pairs, If provided,
            `n_points` will be ignored

        Returns
        -------
        X : array, shape (n_points, )
            X data values for the x-y data pairs
        y : array, shape (n_point, )
            y data values for the x-y data pairs
        """
        if n_points is None and t is None:
            raise AttributeError("A value for `n_points` or `t` must be",
                                 "provided")

        if t is None:
            t = np.linspace(0, 2 * np.pi, n_points)

        center, width, height, phi = self.as_parameters()

        X = (center[0]
             + width * np.cos(t) * np.cos(phi)
             - height * np.sin(t) * np.sin(phi))
        y = (center[1]
             + width * np.cos(t) * np.sin(phi)
             + height * np.sin(t) * np.cos(phi))

        return X, y
