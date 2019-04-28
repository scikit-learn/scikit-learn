# Authors: Alex Li <7Alex7Li@gmail.com>
#          Siyuan Ma <Siyuan.ma9@gmail.com>

import numpy as np
import scipy as sp
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin
from sklearn.metrics.pairwise import pairwise_kernels, euclidean_distances
from sklearn.utils import check_random_state
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.validation import check_is_fitted, check_X_y


class FKR_EigenPro(BaseEstimator, RegressorMixin):
    """Fast kernel regression using EigenPro iteration.

    Train least squared kernel regression model with mini-batch EigenPro
    iteration.

    Parameters
    ----------
    batch_size : int, default = 'auto'
        Mini-batch size for gradient descent.

    n_epoch : int, default = 1
        The number of passes over the training data.

    n_components : int, default = 1000
        the maximum number of eigendirections used in modifying the kernel
        operator. Convergence rate speedup over normal gradient descent is
        approximately the largest eigenvalue over the n_componentth
        eigenvalue, however, it may take time to compute eigenvalues for
        large n_components

    subsample_size : int, default = 'auto'
        The number of subsamples used for estimating the largest
        n_component eigenvalues and eigenvectors. When it is set to 'auto',
        it will be 4000 if there are less than 100,000 samples
        (for training), and otherwise 10000.

    kernel : string or callable, default = "gaussian"
        Kernel mapping used internally. Strings can be anything supported
        by sklearn's library, however, it is recommended to use a radial
        kernel. There is special support for gaussian, laplace, and cauchy
        kernels. A callable should accept two arguments and return a
        floating point number.

    bandwidth : float, default=5
        Bandwidth to use with the gaussian, laplacian, and cauchy kernels.
        Ignored by other kernels.

    gamma : float, default=None
        Gamma parameter for the RBF, polynomial, exponential chi2 and
        sigmoid kernels. Interpretation of the default value is left to
        the kernel; see the documentation for sklearn.metrics.pairwise.
        Ignored by other kernels.

    degree : float, default=3
        Degree of the polynomial kernel. Ignored by other kernels.

    coef0 : float, default=1
        Zero coefficient for polynomial and sigmoid kernels.
        Ignored by other kernels.

    kernel_params : mapping of string to any
        Additional parameters (keyword arguments) for kernel function
        passed as callable object.

    random_state : int, RandomState instance or None, (default=None)
        The seed of the pseudo random number generator to use when
        shuffling the data.  If int, random_state is the seed used by the
        random number generator; If RandomState instance, random_state is
        the random number generator; If None, the random number generator
        is the RandomState instance used by `np.random`.

    References
    ----------
    * Siyuan Ma, Mikhail Belkin
      "Diving into the shallows: a computational perspective on
      large-scale machine learning", NIPS 2017.

    Examples
    --------
    >>> from sklearn.fast_kernel import FKR_EigenPro
    >>> import numpy as np
    >>> n_samples, n_features, n_targets = 4000, 20, 3
    >>> rng = np.random.RandomState(1)
    >>> x_train = rng.randn(n_samples, n_features)
    >>> y_train = rng.randn(n_samples, n_targets)
    >>> rgs = FKR_EigenPro(n_epoch=3, bandwidth=1, subsample_size=50)
    >>> rgs.fit(x_train, y_train)
    FKR_EigenPro(bandwidth=1, batch_size='auto', coef0=1, degree=3, gamma=None,
                 kernel='gaussian', kernel_params=None, n_components=1000,
                 n_epoch=3, random_state=None, subsample_size=50)
    >>> y_pred = rgs.predict(x_train)
    >>> loss = np.mean(np.square(y_train - y_pred))
    """

    def __init__(self, batch_size="auto", n_epoch=2, n_components=1000,
                 subsample_size="auto", kernel="gaussian",
                 bandwidth=5, gamma=None, degree=3, coef0=1,
                 kernel_params=None, random_state=None):
        self.batch_size = batch_size
        self.n_epoch = n_epoch
        self.n_components = n_components
        self.subsample_size = subsample_size
        self.kernel = kernel
        self.bandwidth = bandwidth
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params
        self.random_state = random_state

    def _kernel(self, X, Y, Y_squared=None):
        """Calculate the kernel matrix

        Parameters
        ---------
        X : {float, array}, shape = [n_samples, n_features]
            Input data.

        Y : {float, array}, shape = [n_centers, n_targets]
            Kernel centers.

        Y_squared : {float, array}, shape = [1, n_centers]
            Square of L2 norms of centers.

        Returns
        -------
        K : {float, array}, shape = [n_samples, n_centers]
            Kernel matrix.
        """
        if (self.kernel != "gaussian" and
                self.kernel != "laplace" and
                self.kernel != "cauchy"):
            if callable(self.kernel):
                params = self.kernel_params or {}
            else:
                params = {"gamma": self.gamma,
                          "degree": self.degree,
                          "coef0": self.coef0}
            return pairwise_kernels(X, Y, metric=self.kernel,
                                    filter_params=True, **params)
        distance = euclidean_distances(X, Y, squared=True,
                                       Y_norm_squared=Y_squared)
        bandwidth = np.float32(self.bandwidth)
        if self.kernel == "gaussian":
            K = np.exp(-distance / (2 * np.square(bandwidth)))
        elif self.kernel == "laplace":
            d = np.maximum(distance, 0)
            K = np.exp(-np.sqrt(d) / bandwidth)
        else:  # self.kernel == "cauchy":
            K = 1 / (1 + distance / np.square(bandwidth))
        return K

    def _nystrom_svd(self, X, n_components):
        """Compute the top eigensystem of a kernel
        operator using Nystrom method

        Parameters
        ----------
        X : {float, array}, shape = [n_subsamples, n_features]
            Subsample feature matrix.

        n_components : int
            Number of top eigencomponents to be restored.

        Returns
        -------
        S : {float, array}, shape = [k]
            Top eigenvalues.

        V : {float, array}, shape = [n_subsamples, k]
            Top eigenvectors of a subsample kernel matrix (which can be
            directly used to approximate the eigenfunctions of the kernel
            operator).
        """
        m, _ = X.shape
        K = self._kernel(X, X)
        W = K / m
        S, V = sp.linalg.eigh(W, eigvals=(m - n_components, m - 1))

        # Flip so eigenvalues are in descending order.
        S = np.maximum(np.float32(1e-7), np.flipud(S))
        V = np.fliplr(V)[:, :n_components] / np.sqrt(m, dtype='float32')

        return S, V

    def _setup(self, feat, max_components, mG, alpha):
        """Compute preconditioner and scale factors for EigenPro iteration

        Parameters
        ----------
        feat : {float, array}, shape = [n_samples, n_features]
            Feature matrix (normally from training data).

        max_components : float
            Maximum number of components to be used in EigenPro iteration.

        mG : int
            Maximum batch size to fit in memory.

        alpha : float
            Exponential factor (< 1) for eigenvalue ratio.

        Returns
        -------
        max_S : float
            Normalized largest eigenvalue.

        max_kxx : float
            Maximum of k(x,x) where k is the EigenPro kernel.
        """
        alpha = np.float32(alpha)

        # Estimate eigenvalues (S) and eigenvectors (V) of the kernel matrix
        # corresponding to the feature matrix.
        S, V = self._nystrom_svd(feat, max_components)
        n_subsamples = feat.shape[0]

        # Calculate the number of components to be used such that the
        # corresponding batch size is bounded by the subsample size and the
        # memory size.
        max_bs = min(max(n_subsamples / 5, mG), n_subsamples)
        n_components = np.sum(np.power(1 / S, alpha) < max_bs) - 1
        if n_components < 2:
            n_components = min(S.shape[0] - 1, 2)

        self.V_ = V[:, :n_components]
        scale = np.power(S[0] / S[n_components], alpha)

        # Compute part of the preconditioner for step 2 of gradient descent in
        # the eigenpro model
        self.Q_ = (1 - np.power(S[n_components] / S[:n_components],
                                alpha)) / S[:n_components]

        max_S = S[0].astype(np.float32)
        kxx = 1 - np.sum(self.V_ ** 2, axis=1) * n_subsamples
        return max_S / scale, np.max(kxx)

    def _initialize_params(self, X, Y):
        """Validate parameters passed to the model, choose parameters
        that have not been passed in, and run setup for EigenPro iteration.
        """
        self.random_state_ = check_random_state(self.random_state)
        n, d = X.shape
        n_label = 1 if len(Y.shape) == 1 else Y.shape[1]
        self.centers_ = X

        # Calculate the subsample size to be used.
        if self.subsample_size == "auto":
            if n < 100000:
                sample_size = min(n, 4000)
            else:
                sample_size = 12000
        else:
            sample_size = min(n, self.subsample_size)

        n_components = min(sample_size - 1, self.n_components)
        n_components = max(1, n_components)

        # Each batch will require about 1 gb memory
        mem_bytes = 1024 ** 3
        mem_usages = (d + n_label + 2 * np.arange(sample_size)) * n * 4
        mG = np.sum(mem_usages < mem_bytes)

        # Calculate largest eigenvalue and max{k(x,x)} using subsamples.
        self.pinx_ = self.random_state_.choice(n, sample_size,
                                               replace=False).astype('int32')
        max_S, beta = self._setup(X[self.pinx_], n_components, mG, alpha=.95)
        # Calculate best batch size.
        if self.batch_size == "auto":
            bs = min(np.int32(beta / max_S + 1), mG)
        else:
            bs = self.batch_size
        self.bs_ = min(bs, n)

        # Calculate best step size.
        if self.bs_ < beta / max_S + 1:
            eta = self.bs_ / beta
        elif self.bs_ < n:
            eta = 2. * self.bs_ / (beta + (self.bs_ - 1) * max_S)
        else:
            eta = 0.95 * 2 / max_S
        self.eta_ = np.float32(eta)
        # Remember the shape of Y for predict() and ensure it's shape is 2-D.
        self.was_1D_ = False
        if len(Y.shape) == 1:
            Y = np.reshape(Y, (Y.shape[0], 1))
            self.was_1D_ = True
        return Y

    def fit(self, X, Y):
        """Train fast kernel regression model

        Parameters
        ----------
        X : {float, array}, shape = [n_samples, n_features]
            Training data.

        Y : {float, array}, shape = [n_samples, n_targets]
            Training targets.

        Returns
        -------
        self : returns an instance of self.
        """
        X, Y = check_X_y(X, Y, dtype=np.float32, multi_output=True,
                         ensure_min_samples=3, y_numeric=True)
        Y = Y.astype(np.float32)  # check_X_y does not seem to do this
        """Parameter Initialization"""
        Y = self._initialize_params(X, Y)

        """Training loop"""
        n = self.centers_.shape[0]

        self.coef_ = np.zeros((n, Y.shape[1]), dtype=np.float32)
        self.centers_squared_ = \
            np.square(self.centers_).sum(axis=1, keepdims=True).T
        step = np.float32(self.eta_ / self.bs_)
        for epoch in range(0, self.n_epoch):
            epoch_inds = \
                self.random_state_.choice(n, n // self.bs_ * self.bs_,
                                          replace=False).astype('int32')

            for batch_inds in np.array_split(epoch_inds, n // self.bs_):
                batch_x = self.centers_[batch_inds]
                kfeat = self._kernel(batch_x, self.centers_,
                                     Y_squared=self.centers_squared_)

                batch_y = Y[batch_inds]

                # Update 1: Sampled Coordinate Block.
                gradient = np.dot(kfeat, self.coef_) - batch_y

                self.coef_[batch_inds] = \
                    self.coef_[batch_inds] - step * gradient

                # Update 2: Fixed Coordinate Block
                delta = np.dot(self.V_ * self.Q_,
                               np.dot(self.V_.T, np.dot(
                                   kfeat[:, self.pinx_].T, gradient)))
                self.coef_[self.pinx_] += step * delta
        return self

    def predict(self, X):
        """Predict using the kernel regression model

        Parameters
        ----------
        X : {float, array}, shape = [n_samples, n_features]
            Samples.

        Returns
        -------
        Y : {float, array}, shape = [n_samples, n_targets]
            Predicted targets.
        """
        check_is_fitted(self, ["bs_", "centers_", "centers_squared_", "coef_",
                               "eta_", "random_state_", "pinx_",
                               "Q_", "V_", "was_1D_"])
        X = np.asarray(X, dtype=np.float64)

        if len(X.shape) == 1:
            raise ValueError("Reshape your data. X should be a matrix of shape"
                             " (n_samples, n_features).")
        n = X.shape[0]

        Ys = []
        for batch_inds in np.array_split(range(n), max(1, n // self.bs_)):
            batch_x = X[batch_inds]
            kfeat = self._kernel(batch_x, self.centers_,
                                 Y_squared=self.centers_squared_)

            pred = np.dot(kfeat, self.coef_)
            Ys.append(pred)
        Y = np.vstack(Ys)
        if self.was_1D_:
            Y = np.reshape(Y, Y.shape[0])
        return Y

    def _get_tags(self):
        return {'multioutput': True}


class FKC_EigenPro(BaseEstimator, ClassifierMixin):
    """Fast kernel classification using EigenPro iteration.

    Train least squared kernel classification model with mini-batch EigenPro
    iteration.

    Parameters
    ----------
    batch_size : int, default = 'auto'
        Mini-batch size for gradient descent.

    n_epoch : int, default = 1
        The number of passes over the training data.

    n_components : int, default = 1000
        the maximum number of eigendirections used in modifying the
        kernel operator. Convergence rate speedup over normal gradient
        descent is approximately the largest eigenvalue over the
        n_componenth eigenvalue, however, it may take time to compute
        eigenvalues for large n_components

    subsample_size : int, default = 'auto'
        The size of subsamples used for estimating the largest
        n_component eigenvalues and eigenvectors. When it is set to
        'auto', it will be 4000 if there are less than 100,000 samples
        (for training), and otherwise 10000.

    kernel : string or callable, default = "gaussian"
        Kernel mapping used internally. Strings can be anything
        supported by sklearn's library, however, it is recommended to
        use a radial kernel. There is special support for gaussian,
        laplace, and cauchy kernels. A callable should accept two
        arguments and return a floating point number.

    bandwidth : float, default=5
        Bandwidth to use with the gaussian, laplacian, and cauchy
        kernels. Ignored by other kernels.

    gamma : float, default=None
        Gamma parameter for the RBF, polynomial, exponential chi2
        and sigmoid kernels. Interpretation of the default value is left
        to the kernel; see the documentation for
        sklearn.metrics.pairwise. Ignored by other kernels.

    degree : float, default=3
        Degree of the polynomial kernel. Ignored by other kernels.

    coef0 : float, default=1
        Zero coefficient for polynomial and sigmoid kernels. Ignored by
        other kernels.

    kernel_params : mapping of string to any
        Additional parameters (keyword arguments) for kernel function
        passed as callable object.

    random_state : int, RandomState instance or None (default=None)
        The seed of the pseudo random number generator to use when
        shuffling the data.  If int, random_state is the seed used by
        the random number generator; If RandomState instance,
        random_state is the random number generator;
        If None, the random number generator is the RandomState
        instance used by `np.random`.

    References
    ----------
    * Siyuan Ma, Mikhail Belkin
      "Diving into the shallows: a computational perspective on
      large-scale machine learning", NIPS 2017.

    Examples
    --------
    >>> from sklearn.fast_kernel import FKC_EigenPro
    >>> import numpy as np
    >>> n_samples, n_features, n_targets = 4000, 20, 3
    >>> rng = np.random.RandomState(1)
    >>> x_train = rng.randn(n_samples, n_features)
    >>> y_train = rng.randint(n_targets, size=n_samples)
    >>> rgs = FKC_EigenPro(n_epoch=3, bandwidth=1, subsample_size=50)
    >>> rgs.fit(x_train, y_train)
    FKC_EigenPro(bandwidth=1, batch_size='auto', coef0=1, degree=3, gamma=None,
                 kernel='gaussian', kernel_params=None, n_components=1000,
                 n_epoch=3, random_state=None, subsample_size=50)
    >>> y_pred = rgs.predict(x_train)
    >>> loss = np.mean(y_train != y_pred)
    """

    def __init__(self, batch_size="auto", n_epoch=2, n_components=1000,
                 subsample_size="auto", kernel="gaussian",
                 bandwidth=5, gamma=None, degree=3, coef0=1,
                 kernel_params=None, random_state=None):
        self.batch_size = batch_size
        self.n_epoch = n_epoch
        self.n_components = n_components
        self.subsample_size = subsample_size
        self.kernel = kernel
        self.bandwidth = bandwidth
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params
        self.random_state = random_state

    def fit(self, X, Y):
        """ Train fast kernel classification model

        Parameters
        ----------
        X : {float, array}, shape = [n_samples, n_raw_feature]
            The raw input feature matrix.

        Y : {float, array}, shape =[n_samples]
            The labels corresponding to the features of X.

        Returns
        -------
        self : returns an instance of self.
       """
        self.regressor_ = FKR_EigenPro(
            batch_size=self.batch_size, n_epoch=self.n_epoch,
            n_components=self.n_components,
            subsample_size=self.subsample_size, kernel=self.kernel,
            bandwidth=self.bandwidth, gamma=self.gamma,
            degree=self.degree, coef0=self.coef0,
            kernel_params=self.kernel_params, random_state=self.random_state)
        X, Y = check_X_y(X, Y, multi_output=False, ensure_min_samples=3)
        check_classification_targets(Y)
        self.classes_ = np.unique(Y)

        loc = {}
        for ind, label in enumerate(self.classes_):
            loc[label] = ind

        class_matrix = np.zeros((Y.shape[0], self.classes_.shape[0]))

        for ind, label in enumerate(Y):
            class_matrix[ind][loc[label]] = 1
        self.regressor_.fit(X, class_matrix)

        return self

    def predict(self, X):
        """Predict using the kernel classification model

        Parameters
        ----------
        X : {float, array}, shape = [n_samples, n_features]
            Samples.

        Returns
        -------
        y : {float, array}, shape = [n_samples]
            Predicted labels.
        """
        check_is_fitted(self, ["regressor_"])
        Y = self.regressor_.predict(X)
        return self.classes_[np.argmax(Y, axis=1)]

    def _get_tags(self):
        return {'multioutput': True}
