import math
import numpy as np
from scipy import optimize
from .._shared.utils import check_random_state


def _check_data_dim(data, dim):
    if data.ndim != 2 or data.shape[1] != dim:
        raise ValueError('Input data must have shape (N, %d).' % dim)


def _check_data_atleast_2D(data):
    if data.ndim < 2 or data.shape[1] < 2:
        raise ValueError('Input data must be at least 2D.')


def _norm_along_axis(x, axis):
    """NumPy < 1.8 does not support the `axis` argument for `np.linalg.norm`."""
    return np.sqrt(np.einsum('ij,ij->i', x, x))


class BaseModel(object):

    def __init__(self):
        self.params = None


class LineModelND(BaseModel):
    """Total least squares estimator for N-dimensional lines.

    In contrast to ordinary least squares line estimation, this estimator
    minimizes the orthogonal distances of points to the estimated line.

    Lines are defined by a point (origin) and a unit vector (direction)
    according to the following vector equation::

        X = origin + lambda * direction

    Attributes
    ----------
    params : tuple
        Line model parameters in the following order `origin`, `direction`.

    Examples
    --------
    >>> x = np.linspace(1, 2, 25)
    >>> y = 1.5 * x + 3
    >>> lm = LineModelND()
    >>> lm.estimate(np.array([x, y]).T)
    True
    >>> tuple(np.round(lm.params, 5))
    (array([ 1.5 ,  5.25]), array([ 0.5547 ,  0.83205]))
    >>> res = lm.residuals(np.array([x, y]).T)
    >>> np.abs(np.round(res, 9))
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
    >>> np.round(lm.predict_y(x[:5]), 3)
    array([ 4.5  ,  4.562,  4.625,  4.688,  4.75 ])
    >>> np.round(lm.predict_x(y[:5]), 3)
    array([ 1.   ,  1.042,  1.083,  1.125,  1.167])

    """

    def estimate(self, data):
        """Estimate line model from data.

        This minimizes the sum of shortest (orthogonal) distances
        from the given data points to the estimated line.

        Parameters
        ----------
        data : (N, dim) array
            N points in a space of dimensionality dim >= 2.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.
        """
        _check_data_atleast_2D(data)

        origin = data.mean(axis=0)
        data = data - origin

        if data.shape[0] == 2:  # well determined
            direction = data[1] - data[0]
            norm = np.linalg.norm(direction)
            if norm != 0:  # this should not happen to be norm 0
                direction /= norm
        elif data.shape[0] > 2:  # over-determined
            # Note: with full_matrices=1 Python dies with joblib parallel_for.
            _, _, v = np.linalg.svd(data, full_matrices=False)
            direction = v[0]
        else:  # under-determined
            raise ValueError('At least 2 input points needed.')

        self.params = (origin, direction)

        return True

    def residuals(self, data, params=None):
        """Determine residuals of data to model.

        For each point, the shortest (orthogonal) distance to the line is
        returned. It is obtained by projecting the data onto the line.

        Parameters
        ----------
        data : (N, dim) array
            N points in a space of dimension dim.
        params : (2, ) array, optional
            Optional custom parameter set in the form (`origin`, `direction`).

        Returns
        -------
        residuals : (N, ) array
            Residual for each data point.
        """
        _check_data_atleast_2D(data)
        if params is None:
            params = self.params
        assert params is not None
        if len(params) != 2:
            raise ValueError('Parameters are defined by 2 sets.')

        origin, direction = params
        res = (data - origin) - \
              np.dot(data - origin, direction)[..., np.newaxis] * direction
        return _norm_along_axis(res, axis=1)

    def predict(self, x, axis=0, params=None):
        """Predict intersection of the estimated line model with a hyperplane
        orthogonal to a given axis.

        Parameters
        ----------
        x : (n, 1) array
            Coordinates along an axis.
        axis : int
            Axis orthogonal to the hyperplane intersecting the line.
        params : (2, ) array, optional
            Optional custom parameter set in the form (`origin`, `direction`).

        Returns
        -------
        data : (n, m) array
            Predicted coordinates.

        Raises
        ------
        ValueError
            If the line is parallel to the given axis.
        """
        if params is None:
            params = self.params
        assert params is not None
        if len(params) != 2:
            raise ValueError('Parameters are defined by 2 sets.')

        origin, direction = params

        if direction[axis] == 0:
            # line parallel to axis
            raise ValueError('Line parallel to axis %s' % axis)

        l = (x - origin[axis]) / direction[axis]
        data = origin + l[..., np.newaxis] * direction
        return data

    def predict_x(self, y, params=None):
        """Predict x-coordinates for 2D lines using the estimated model.

        Alias for::

            predict(y, axis=1)[:, 0]

        Parameters
        ----------
        y : array
            y-coordinates.
        params : (2, ) array, optional
            Optional custom parameter set in the form (`origin`, `direction`).

        Returns
        -------
        x : array
            Predicted x-coordinates.

        """
        x = self.predict(y, axis=1, params=params)[:, 0]
        return x

    def predict_y(self, x, params=None):
        """Predict y-coordinates for 2D lines using the estimated model.

        Alias for::

            predict(x, axis=0)[:, 1]

        Parameters
        ----------
        x : array
            x-coordinates.
        params : (2, ) array, optional
            Optional custom parameter set in the form (`origin`, `direction`).

        Returns
        -------
        y : array
            Predicted y-coordinates.

        """
        y = self.predict(x, axis=0, params=params)[:, 1]
        return y


class CircleModel(BaseModel):

    """Total least squares estimator for 2D circles.

    The functional model of the circle is::

        r**2 = (x - xc)**2 + (y - yc)**2

    This estimator minimizes the squared distances from all points to the
    circle::

        min{ sum((r - sqrt((x_i - xc)**2 + (y_i - yc)**2))**2) }

    A minimum number of 3 points is required to solve for the parameters.

    Attributes
    ----------
    params : tuple
        Circle model parameters in the following order `xc`, `yc`, `r`.

    Examples
    --------
    >>> t = np.linspace(0, 2 * np.pi, 25)
    >>> xy = CircleModel().predict_xy(t, params=(2, 3, 4))
    >>> model = CircleModel()
    >>> model.estimate(xy)
    True
    >>> tuple(np.round(model.params, 5))
    (2.0, 3.0, 4.0)
    >>> res = model.residuals(xy)
    >>> np.abs(np.round(res, 9))
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

    """

    def estimate(self, data):
        """Estimate circle model from data using total least squares.

        Parameters
        ----------
        data : (N, 2) array
            N points with ``(x, y)`` coordinates, respectively.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        _check_data_dim(data, dim=2)

        x = data[:, 0]
        y = data[:, 1]

        # http://www.had2know.com/academics/best-fit-circle-least-squares.html
        x2y2 = (x ** 2 + y ** 2)
        sum_x = np.sum(x)
        sum_y = np.sum(y)
        sum_xy = np.sum(x * y)
        m1 = np.array([[np.sum(x ** 2), sum_xy, sum_x],
                       [sum_xy, np.sum(y ** 2), sum_y],
                       [sum_x, sum_y, float(len(x))]])
        m2 = np.array([[np.sum(x * x2y2),
                        np.sum(y * x2y2),
                        np.sum(x2y2)]]).T
        a, b, c = np.linalg.pinv(m1).dot(m2)
        a, b, c = a[0], b[0], c[0]
        xc = a / 2
        yc = b / 2
        r = np.sqrt(4 * c + a ** 2 + b ** 2) / 2

        self.params = (xc, yc, r)

        return True

    def residuals(self, data):
        """Determine residuals of data to model.

        For each point the shortest distance to the circle is returned.

        Parameters
        ----------
        data : (N, 2) array
            N points with ``(x, y)`` coordinates, respectively.

        Returns
        -------
        residuals : (N, ) array
            Residual for each data point.

        """

        _check_data_dim(data, dim=2)

        xc, yc, r = self.params

        x = data[:, 0]
        y = data[:, 1]

        return r - np.sqrt((x - xc)**2 + (y - yc)**2)

    def predict_xy(self, t, params=None):
        """Predict x- and y-coordinates using the estimated model.

        Parameters
        ----------
        t : array
            Angles in circle in radians. Angles start to count from positive
            x-axis to positive y-axis in a right-handed system.
        params : (3, ) array, optional
            Optional custom parameter set.

        Returns
        -------
        xy : (..., 2) array
            Predicted x- and y-coordinates.

        """
        if params is None:
            params = self.params
        xc, yc, r = params

        x = xc + r * np.cos(t)
        y = yc + r * np.sin(t)

        return np.concatenate((x[..., None], y[..., None]), axis=t.ndim)


class EllipseModel(BaseModel):
    """Total least squares estimator for 2D ellipses.

    The functional model of the ellipse is::

        xt = xc + a*cos(theta)*cos(t) - b*sin(theta)*sin(t)
        yt = yc + a*sin(theta)*cos(t) + b*cos(theta)*sin(t)
        d = sqrt((x - xt)**2 + (y - yt)**2)

    where ``(xt, yt)`` is the closest point on the ellipse to ``(x, y)``. Thus
    d is the shortest distance from the point to the ellipse.

    The estimator is based on a least squares minimization. The optimal
    solution is computed directly, no iterations are required. This leads
    to a simple, stable and robust fitting method.

    The ``params`` attribute contains the parameters in the following order::

        xc, yc, a, b, theta

    Attributes
    ----------
    params : tuple
        Ellipse model parameters in the following order `xc`, `yc`, `a`, `b`,
        `theta`.

    Examples
    --------

    >>> xy = EllipseModel().predict_xy(np.linspace(0, 2 * np.pi, 25),
    ...                                params=(10, 15, 4, 8, np.deg2rad(30)))
    >>> ellipse = EllipseModel()
    >>> ellipse.estimate(xy)
    True
    >>> np.round(ellipse.params, 2)
    array([ 10.  ,  15.  ,   4.  ,   8.  ,   0.52])
    >>> np.round(abs(ellipse.residuals(xy)), 5)
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
    """

    def estimate(self, data):
        """Estimate circle model from data using total least squares.

        Parameters
        ----------
        data : (N, 2) array
            N points with ``(x, y)`` coordinates, respectively.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.


        References
        ----------
        .. [1] Halir, R.; Flusser, J. "Numerically stable direct least squares
               fitting of ellipses". In Proc. 6th International Conference in
               Central Europe on Computer Graphics and Visualization.
               WSCG (Vol. 98, pp. 125-132).

        """
        # Original Implementation: Ben Hammel, Nick Sullivan-Molina
        # another REFERENCE: [2] http://mathworld.wolfram.com/Ellipse.html
        _check_data_dim(data, dim=2)

        x = data[:, 0]
        y = data[:, 1]

        # Quadratic part of design matrix [eqn. 15] from [1]
        D1 = np.vstack([x ** 2, x * y, y ** 2]).T
        # Linear part of design matrix [eqn. 16] from [1]
        D2 = np.vstack([x, y, np.ones(len(x))]).T

        # forming scatter matrix [eqn. 17] from [1]
        S1 = np.dot(D1.T, D1)
        S2 = np.dot(D1.T, D2)
        S3 = np.dot(D2.T, D2)

        # Constraint matrix [eqn. 18]
        C1 = np.array([[0., 0., 2.], [0., -1., 0.], [2., 0., 0.]])

        try:
            # Reduced scatter matrix [eqn. 29]
            M = np.linalg.inv(C1).dot(
                S1 - np.dot(S2, np.linalg.inv(S3)).dot(S2.T))
        except np.linalg.LinAlgError:  # LinAlgError: Singular matrix
            return False

        # M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors
        # from this equation [eqn. 28]
        eig_vals, eig_vecs = np.linalg.eig(M)

        # eigenvector must meet constraint 4ac - b^2 to be valid.
        cond = 4 * np.multiply(eig_vecs[0, :], eig_vecs[2, :]) \
               - np.power(eig_vecs[1, :], 2)
        a1 = eig_vecs[:, (cond > 0)]
        # seeks for empty matrix
        if 0 in a1.shape or len(a1.ravel()) != 3:
            return False
        a, b, c = a1.ravel()

        # |d f g> = -S3^(-1)*S2^(T)*|a b c> [eqn. 24]
        a2 = np.dot(-np.linalg.inv(S3), S2.T).dot(a1)
        d, f, g = a2.ravel()

        # eigenvectors are the coefficients of an ellipse in general form
        # a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0 (eqn. 15) from [2]
        b /= 2.
        d /= 2.
        f /= 2.

        # finding center of ellipse [eqn.19 and 20] from [2]
        x0 = (c * d - b * f) / (b ** 2. - a * c)
        y0 = (a * f - b * d) / (b ** 2. - a * c)

        # Find the semi-axes lengths [eqn. 21 and 22] from [2]
        numerator = a * f ** 2 + c * d ** 2 + g * b ** 2 \
                    - 2 * b * d * f - a * c * g
        term = np.sqrt((a - c) ** 2 + 4 * b ** 2)
        denominator1 = (b ** 2 - a * c) * (term - (a + c))
        denominator2 = (b ** 2 - a * c) * (- term - (a + c))
        width = np.sqrt(2 * numerator / denominator1)
        height = np.sqrt(2 * numerator / denominator2)

        # angle of counterclockwise rotation of major-axis of ellipse
        # to x-axis [eqn. 23] from [2].
        phi = 0.5 * np.arctan((2. * b) / (a - c))
        if a > c:
            phi += 0.5 * np.pi

        self.params = np.nan_to_num([x0, y0, width, height, phi]).tolist()
        self.params = [float(np.real(x)) for x in self.params]
        return True

    def residuals(self, data):
        """Determine residuals of data to model.

        For each point the shortest distance to the ellipse is returned.

        Parameters
        ----------
        data : (N, 2) array
            N points with ``(x, y)`` coordinates, respectively.

        Returns
        -------
        residuals : (N, ) array
            Residual for each data point.

        """

        _check_data_dim(data, dim=2)

        xc, yc, a, b, theta = self.params

        ctheta = math.cos(theta)
        stheta = math.sin(theta)

        x = data[:, 0]
        y = data[:, 1]

        N = data.shape[0]

        def fun(t, xi, yi):
            ct = math.cos(t)
            st = math.sin(t)
            xt = xc + a * ctheta * ct - b * stheta * st
            yt = yc + a * stheta * ct + b * ctheta * st
            return (xi - xt) ** 2 + (yi - yt) ** 2

        # def Dfun(t, xi, yi):
        #     ct = math.cos(t)
        #     st = math.sin(t)
        #     xt = xc + a * ctheta * ct - b * stheta * st
        #     yt = yc + a * stheta * ct + b * ctheta * st
        #     dfx_t = - 2 * (xi - xt) * (- a * ctheta * st
        #                                - b * stheta * ct)
        #     dfy_t = - 2 * (yi - yt) * (- a * stheta * st
        #                                + b * ctheta * ct)
        #     return [dfx_t + dfy_t]

        residuals = np.empty((N, ), dtype=np.double)

        # initial guess for parameter t of closest point on ellipse
        t0 = np.arctan2(y - yc, x - xc) - theta

        # determine shortest distance to ellipse for each point
        for i in range(N):
            xi = x[i]
            yi = y[i]
            # faster without Dfun, because of the python overhead
            t, _ = optimize.leastsq(fun, t0[i], args=(xi, yi))
            residuals[i] = np.sqrt(fun(t, xi, yi))

        return residuals

    def predict_xy(self, t, params=None):
        """Predict x- and y-coordinates using the estimated model.

        Parameters
        ----------
        t : array
            Angles in circle in radians. Angles start to count from positive
            x-axis to positive y-axis in a right-handed system.
        params : (5, ) array, optional
            Optional custom parameter set.

        Returns
        -------
        xy : (..., 2) array
            Predicted x- and y-coordinates.

        """

        if params is None:
            params = self.params

        xc, yc, a, b, theta = params

        ct = np.cos(t)
        st = np.sin(t)
        ctheta = math.cos(theta)
        stheta = math.sin(theta)

        x = xc + a * ctheta * ct - b * stheta * st
        y = yc + a * stheta * ct + b * ctheta * st

        return np.concatenate((x[..., None], y[..., None]), axis=t.ndim)


def _dynamic_max_trials(n_inliers, n_samples, min_samples, probability):
    """Determine number trials such that at least one outlier-free subset is
    sampled for the given inlier/outlier ratio.
    Parameters
    ----------
    n_inliers : int
        Number of inliers in the data.
    n_samples : int
        Total number of samples in the data.
    min_samples : int
        Minimum number of samples chosen randomly from original data.
    probability : float
        Probability (confidence) that one outlier-free sample is generated.
    Returns
    -------
    trials : int
        Number of trials.
    """
    if n_inliers == 0:
        return np.inf

    nom = 1 - probability
    if nom == 0:
        return np.inf

    inlier_ratio = n_inliers / float(n_samples)
    denom = 1 - inlier_ratio ** min_samples
    if denom == 0:
        return 1
    elif denom == 1:
        return np.inf

    nom = np.log(nom)
    denom = np.log(denom)
    if denom == 0:
        return 0

    return int(np.ceil(nom / denom))


def ransac(data, model_class, min_samples, residual_threshold,
           is_data_valid=None, is_model_valid=None,
           max_trials=100, stop_sample_num=np.inf, stop_residuals_sum=0,
           stop_probability=1, random_state=None):
    """Fit a model to data with the RANSAC (random sample consensus) algorithm.

    RANSAC is an iterative algorithm for the robust estimation of parameters
    from a subset of inliers from the complete data set. Each iteration
    performs the following tasks:

    1. Select `min_samples` random samples from the original data and check
       whether the set of data is valid (see `is_data_valid`).
    2. Estimate a model to the random subset
       (`model_cls.estimate(*data[random_subset]`) and check whether the
       estimated model is valid (see `is_model_valid`).
    3. Classify all data as inliers or outliers by calculating the residuals
       to the estimated model (`model_cls.residuals(*data)`) - all data samples
       with residuals smaller than the `residual_threshold` are considered as
       inliers.
    4. Save estimated model as best model if number of inlier samples is
       maximal. In case the current estimated model has the same number of
       inliers, it is only considered as the best model if it has less sum of
       residuals.

    These steps are performed either a maximum number of times or until one of
    the special stop criteria are met. The final model is estimated using all
    inlier samples of the previously determined best model.

    Parameters
    ----------
    data : [list, tuple of] (N, D) array
        Data set to which the model is fitted, where N is the number of data
        points and D the dimensionality of the data.
        If the model class requires multiple input data arrays (e.g. source and
        destination coordinates of  ``skimage.transform.AffineTransform``),
        they can be optionally passed as tuple or list. Note, that in this case
        the functions ``estimate(*data)``, ``residuals(*data)``,
        ``is_model_valid(model, *random_data)`` and
        ``is_data_valid(*random_data)`` must all take each data array as
        separate arguments.
    model_class : object
        Object with the following object methods:

         * ``success = estimate(*data)``
         * ``residuals(*data)``

        where `success` indicates whether the model estimation succeeded
        (`True` or `None` for success, `False` for failure).
    min_samples : int
        The minimum number of data points to fit a model to.
    residual_threshold : float
        Maximum distance for a data point to be classified as an inlier.
    is_data_valid : function, optional
        This function is called with the randomly selected data before the
        model is fitted to it: `is_data_valid(*random_data)`.
    is_model_valid : function, optional
        This function is called with the estimated model and the randomly
        selected data: `is_model_valid(model, *random_data)`, .
    max_trials : int, optional
        Maximum number of iterations for random sample selection.
    stop_sample_num : int, optional
        Stop iteration if at least this number of inliers are found.
    stop_residuals_sum : float, optional
        Stop iteration if sum of residuals is less than or equal to this
        threshold.
    stop_probability : float in range [0, 1], optional
        RANSAC iteration stops if at least one outlier-free set of the
        training data is sampled with ``probability >= stop_probability``,
        depending on the current best model's inlier ratio and the number
        of trials. This requires to generate at least N samples (trials):

            N >= log(1 - probability) / log(1 - e**m)

        where the probability (confidence) is typically set to a high value
        such as 0.99, and e is the current fraction of inliers w.r.t. the
        total number of samples.
    random_state : int, RandomState instance or None, optional
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.


    Returns
    -------
    model : object
        Best model with largest consensus set.
    inliers : (N, ) array
        Boolean mask of inliers classified as ``True``.

    References
    ----------
    .. [1] "RANSAC", Wikipedia, http://en.wikipedia.org/wiki/RANSAC

    Examples
    --------

    Generate ellipse data without tilt and add noise:

    >>> t = np.linspace(0, 2 * np.pi, 50)
    >>> xc, yc = 20, 30
    >>> a, b = 5, 10
    >>> x = xc + a * np.cos(t)
    >>> y = yc + b * np.sin(t)
    >>> data = np.column_stack([x, y])
    >>> np.random.seed(seed=1234)
    >>> data += np.random.normal(size=data.shape)

    Add some faulty data:

    >>> data[0] = (100, 100)
    >>> data[1] = (110, 120)
    >>> data[2] = (120, 130)
    >>> data[3] = (140, 130)

    Estimate ellipse model using all available data:

    >>> model = EllipseModel()
    >>> model.estimate(data)
    True
    >>> np.round(model.params)  # doctest: +SKIP
    array([ 72.,  75.,  77.,  14.,   1.])

    Estimate ellipse model using RANSAC:

    >>> ransac_model, inliers = ransac(data, EllipseModel, 20, 3, max_trials=50)
    >>> abs(np.round(ransac_model.params))
    array([ 20.,  30.,   5.,  10.,   0.])
    >>> inliers # doctest: +SKIP
    array([False, False, False, False,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True], dtype=bool)
    >>> sum(inliers) > 40
    True

    Robustly estimate geometric transformation:

    >>> from skimage.transform import SimilarityTransform
    >>> np.random.seed(0)
    >>> src = 100 * np.random.rand(50, 2)
    >>> model0 = SimilarityTransform(scale=0.5, rotation=1,
    ...                              translation=(10, 20))
    >>> dst = model0(src)
    >>> dst[0] = (10000, 10000)
    >>> dst[1] = (-100, 100)
    >>> dst[2] = (50, 50)
    >>> model, inliers = ransac((src, dst), SimilarityTransform, 2, 10)
    >>> inliers
    array([False, False, False,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True], dtype=bool)

    """

    best_model = None
    best_inlier_num = 0
    best_inlier_residuals_sum = np.inf
    best_inliers = None

    random_state = check_random_state(random_state)

    if min_samples < 0:
        raise ValueError("`min_samples` must be greater than zero")

    if max_trials < 0:
        raise ValueError("`max_trials` must be greater than zero")

    if stop_probability < 0 or stop_probability > 1:
        raise ValueError("`stop_probability` must be in range [0, 1]")

    if not isinstance(data, list) and not isinstance(data, tuple):
        data = [data]

    # make sure data is list and not tuple, so it can be modified below
    data = list(data)
    # number of samples
    num_samples = data[0].shape[0]

    for num_trials in range(max_trials):

        # choose random sample set
        samples = []
        random_idxs = random_state.randint(0, num_samples, min_samples)
        for d in data:
            samples.append(d[random_idxs])

        # check if random sample set is valid
        if is_data_valid is not None and not is_data_valid(*samples):
            continue

        # estimate model for current random sample set
        sample_model = model_class()

        success = sample_model.estimate(*samples)

        if success is not None:  # backwards compatibility
            if not success:
                continue

        # check if estimated model is valid
        if is_model_valid is not None \
                and not is_model_valid(sample_model, *samples):
            continue

        sample_model_residuals = np.abs(sample_model.residuals(*data))
        # consensus set / inliers
        sample_model_inliers = sample_model_residuals < residual_threshold
        sample_model_residuals_sum = np.sum(sample_model_residuals**2)

        # choose as new best model if number of inliers is maximal
        sample_inlier_num = np.sum(sample_model_inliers)
        if (
            # more inliers
            sample_inlier_num > best_inlier_num
            # same number of inliers but less "error" in terms of residuals
            or (sample_inlier_num == best_inlier_num
                and sample_model_residuals_sum < best_inlier_residuals_sum)
        ):
            best_model = sample_model
            best_inlier_num = sample_inlier_num
            best_inlier_residuals_sum = sample_model_residuals_sum
            best_inliers = sample_model_inliers
            if (
                best_inlier_num >= stop_sample_num
                or best_inlier_residuals_sum <= stop_residuals_sum
                or num_trials
                    >= _dynamic_max_trials(best_inlier_num, num_samples,
                                           min_samples, stop_probability)
            ):
                break

    # estimate final model using all inliers
    if best_inliers is not None:
        # select inliers for each data array
        for i in range(len(data)):
            data[i] = data[i][best_inliers]
        best_model.estimate(*data)

    return best_model, best_inliers
