from scipy.optimize import leastsq
import numpy as np
from scipy.linalg import solve_triangular

from ._base import RegressorMixin
from ..preprocessing import PolynomialFeatures


def exp_kernel(z):
    """
    Implements the exponential kernel $$ e^{-|z|^2/2}$$
    as a vectorized function.
    :param z: a numpy array of any dimension
    :return: the exponential kernel evaluated on the array z assuming
             the last axis is the dimension of the elements of z
    """
    return np.exp(-np.linalg.norm(np.atleast_1d(z), axis=-1) ** 2 / 2)


class LocalRegression(RegressorMixin):
    """
    This class performs local polynomial regression of dimension $d$
    and degree $p$, i.e. given a training set $(X,y)$ of $N$ samples
    $x_i \in \R^d$ and $y_i \in \R$, it makes a new prediction at
    $x \in \R^d$ as $x=\beta^0$, where $\beta$ is the value minimizing
    the cost functional
    \begin{align*}
        J(\beta) & := \sum_{i=1}^{N}{(y_i - \beta ^0 + j_1(x_i))^2 w_i},
        j_1(x_i) & := \sum_{1 \leq |\alpha| \leq p}{(x-x_i)^\alpha}, \\
        w_i & := K_h(x-x_i),
    \end{align*}
    where $K_h$ is a kernel function scaled with bandwith $h$. The minimum
    of the cost function is computed via QR decomposition (or analytically).
    """

    def __init__(self, degree, kernel=exp_kernel, warm_start=True):
        """
        :param degree: an int specifying the degree of the polynomial
        :param kernel: a kernel function for the weight calculation
        :param warm_start: if True the fitted kernel is preserved enabling
                           partial_fit() for new y
        """
        self.degree = degree
        self.kernel = kernel
        self.warm_start = warm_start
        self.bandwidth = None
        self.fitted_kernel = None
        self.method = None
        self.X_train_ = None
        self.y_train = None

    def fit(self, X, y, bandwidth=None):
        """
        Fits the regressor to the data. As the concept of local regression is
        to fit the data to each prediction, this function only stores the data
        and either sets a fixed bandwidth or estimates an optimal one.
        :param X: a numpy array of shape (N, d)
        :param y: a numpy array of shape (N,)
        :param bandwidth: a scalar or None
        :return: self
        """
        self.X_train_ = X
        self.y_train = y
        self._set_bandwidth(bandwidth)
        return self

    def predict(self, X, method=None):
        """
        Performs the prediction for each value x in the prediction set X.
        If $d=1$ and $p=0,1$ the cost functional can be minimized analytically.
        :param X: a numpy array of dimension (M, d)
        :param method: can be 'analytic' or 'qr' or 'leastsq' or None. Method is
                       set automatically if None.
        :return: a numpy array y of dimension (M,) with the predictions
        """
        if method is None:
            method = self._determine_method()
        if self.fitted_kernel is not None:
            self.fitted_kernel = None
            self.method = None
        return self._predict_with_method(X, method)

    def fit_partial(self, y, bandwidth=None):
        """ Re-fits only the y values of the regression. Only works if
            warm_start==True and a previous (full) fit and predict
            has already been performed.
        """
        if self.warm_start and self.fitted_kernel is not None and self.method in ['analytic', 'qr']:
            self.y_train = y
            self._set_bandwidth(bandwidth)
        else:
            raise ValueError("The fit_partial method can only be invoked\
                if fit and predict have been invoked previously with method \
                    `analytic´ or `qr` and warm_start is set to True")

    def predict_partial(self):
        """ Predicts on the last value of y_train set by fit() or fit_partial()
            and the last X that has been used for prediction. """
        if self.fitted_kernel is not None and self.method is not None:
            return self._predict_with_method(X=None, method=self.method)
        else:
            raise ValueError("The method predict_partial requires a full\
                              prior run of predict with method `analytic´ \
                              or `qr´.")

    def _set_bandwidth(self, bandwidth):
        """ Sets the bandwidth in the fitting."""
        if bandwidth is None:
            self.bandwidth = self._silverman()
        else:
            self.bandwidth = bandwidth

    def _silverman(self):
        """
        This function implements Silverman's Rule of Thumb
        \begin{align*}
            h = \Big( \frac{4}{3n} \Big)^{\frac{1}{5}} \hat \sigma_Y
        \end{align*}
        to estimate the optimal bandwidth of the training data y.
        :return: bandwidth h
        """
        sigma_y = np.std(self.y_train)
        n = self.y_train.shape[0]
        return (4 / (3 * n)) ** (1 / 5) * sigma_y

    def _scaled_kernel(self):
        """
        Scales the kernel function self.kernel by the
        bandwidth self.bandwidth.
        :return: scaled kernel function
        """

        def kh(z):
            d = self.X_train_.shape[1]
            return self.kernel(z / self.bandwidth) / self.bandwidth ** d

        return kh

    def _predict_with_method(self, X, method):
        """ Performs the prediction based on the `method´ flag.
        :param method: a string chosen from 'analytic', 'leastsq', 'qr'
        """
        self.method = method
        num_dims = self.X_train_.shape[1]
        if method == 'analytic':
            if self.degree == 0:
                return self._predict_nadaraya_watson(X)
            elif num_dims == 1 and self.degree == 1:
                return self._predict_locally_linear(X)
            else:
                raise ValueError('Method `analytic´ is only available if \
                                  self.degree=0 or self.degree=1 and X.shape[1] == 1.')
        elif method == 'qr':
            return self._predict_qr(X)
        elif method == 'leastsq':
            return self._predict_leastsq(X)
        else:
            raise ValueError('Parameter `method´ has to be `analytic´ or\
                `qr´ or `leastsq´ or None, but is currently set to: %s' % method)

    def _determine_method(self):
        """ Automatically selects the method for prediction based on the dimension
            of the training data self_X_train_.
        """
        num_dims = self.X_train_.shape[1]
        if self.degree == 0 or (self.degree == 1 and num_dims == 1):
            return 'analytic'
        else:
            return 'qr'

    def _predict_nadaraya_watson(self, X):
        """
        Performs local regression of degree $p=0$ in dimension $d=1$.
        In this case, the cost functional can be minimized analytically
        and for any $x \in \R$, the estimate $y$ is given by
        \begin{align*}
            y &= \sum_{i=1}^{N}{W^i_h(x) y_i}, \\
            W^0_i(x) &= \frac{K_h(x - x_i)}{\sum_{j=1}^N{K_h(x - x_j)}}
        \end{align*}
        :param X: a numpy array of dimension (M, d) at which to predict
        :return: a numpy array of dimension (M,) with the M predicted y's
        """
        if (self.warm_start and self.fitted_kernel is None) or (not self.warm_start):
            k_h = self._scaled_kernel()
            nwk = k_h(np.array([x - self.X_train_ for x in X]))
            self.fitted_kernel = nwk / np.sum(nwk, axis=1)[:, np.newaxis]
        return np.sum(self.fitted_kernel * self.y_train.squeeze(), axis=1)

    def _predict_locally_linear(self, X):
        """
        Performs local regression of degree $p=1$ in dimension $d=1$.
        In this case, the cost functional can be minimized analytically
        and for any $x \in \R$, the estimate $y$ is given by
        \begin{align*}
            y &= \sum_{i=1}^{N}{W^i_h(x) y_i}, \\
            W^1_i(x) &:= \frac{K_h(x-x_i)}{N}\frac{s_2(x)-s_1(x)(x-x_i)}{x_2(x)s_0(x)} \\
            s_r(x) &:= \frac{1}{N} \sum_{i=1}^N{(x-x_i)^rK_h(x-x_i)}.
        \end{align*}
        :param X: a numpy array of dimension (M, d) at which to predict
        :return: a numpy array of dimension (M,) with the M predicted y's
        """
        if (self.warm_start and self.fitted_kernel is None) or (not self.warm_start):
            n = self.X_train_.shape[0]
            m = X.shape[0]
            self.fitted_kernel = np.zeros((m, n))
            k_h = self._scaled_kernel()
            X_ = self.X_train_.squeeze()
            for i in range(m):
                x = X[i]
                llk = k_h(x - self.X_train_).squeeze()
                s0 = np.mean(llk)
                s1 = np.mean((x - X_) * llk)
                s2 = np.mean((x - X_) ** 2 * llk)
                s = (s2 - s1 * (x - X_)) / (s2 * s0 - s1 ** 2) / n
                self.fitted_kernel[i, :] = s * llk
        return np.sum(self.fitted_kernel * self.y_train.squeeze(), axis=1)

    def _predict_qr(self, X):
        """
        Performs a prediction for each x in X by solving the associated
        normal equations via QR decomposition.

        :param X: a numpy array of dimension (M, d) at which to predict
        :return: a numpy array of dimension (M,) with the M predicted y's
        """
        if (self.warm_start and self.fitted_kernel is None) or (not self.warm_start):
            self.fitted_kernel = []
            poly = PolynomialFeatures(degree=self.degree)
            for i in range(X.shape[0]):
                x = X[i, :]
                phi = poly.fit_transform(self.X_train_ - x)
                m = phi.shape[1]
                kh = self._scaled_kernel()
                w = kh(x - self.X_train_)
                w_mat = np.diag(np.sqrt(w))
                a = w_mat @ phi
                q, r = np.linalg.qr(a, mode='complete')
                r = r[:m, :m]
                self.fitted_kernel.append((q.transpose() @ w_mat, r))
        n = len(self.fitted_kernel)
        y_pred = np.zeros(n)
        for i in range(n):
            qw, r = self.fitted_kernel[i]
            m = r.shape[1]
            c = (qw @ self.y_train)[:m]
            beta = solve_triangular(r, c)
            y_pred[i] = beta[0]
        return y_pred

    def _predict_leastsq(self, X):
        """
        Performs a prediction for each x in X by minimizing the cost
        functional $J$.

        :param X: a numpy array of dimension (M, d) at which to predict
        :return: a numpy array of dimension (M,) with the M predicted y's
        """
        poly = PolynomialFeatures(degree=self.degree)
        x0 = poly.fit_transform(np.zeros((1, X[0].shape[0])))[0]
        return np.array([leastsq(func=self._cost_functional(poly, x),
                                 x0=x0)[0][0] for x in X])

    def _cost_functional(self, poly, x):
        """
        Creates the cost functional $J$ for optimization using the parameters.
        :param poly: an instance of PolynomialFeatures
        :param x: a numpy array of shape (d,)
        :return: cost funtional $J$
        """

        def cost(beta):
            res = np.sum(poly.fit_transform(self.X_train_ - x) * beta, axis=1) - self.y_train
            kh = self._scaled_kernel()
            w = kh(x - self.X_train_)
            res *= np.sqrt(np.abs(w))
            return res

        return cost
