import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_X_y
from sklearn.utils.multiclass import check_classification_targets
from sklearn.fast_kernel_regression import FastKernelRegression


class FastKernelClassification(BaseEstimator, ClassifierMixin):
    """Fast kernel classification.

       Train least squared kernel classification model with mini-batch EigenPro
       iteration.

       Parameters
       ----------
           bs: int, default = 'auto'
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

           mem_gb : int, default = 12
               Physical device memory in GB.

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

           kernel_params : mapping of string to any, optional
               Additional parameters (keyword arguments) for kernel function
               passed as callable object.

           random_state : int, RandomState instance or None, optional (default=None)
               The seed of the pseudo random number generator to use when shuffling
               the data.  If int, random_state is the seed used by the random number
               generator; If RandomState instance, random_state is the random number
               generator; If None, the random number generator is the RandomState
               instance used by `np.random`.

       References
       ----------
       * Siyuan Ma, Mikhail Belkin
         "Diving into the shallows: a computational perspective on
         large-scale machine learning", NIPS 2017.

       Examples
       --------
           >>> from sklearn.fast_kernel_classification import FastKernelClassification
           >>> import numpy as np
           >>> n_samples, n_features, n_targets = 4000, 20, 3
           >>> rng = np.random.RandomState(1)
           >>> x_train = rng.randn(n_samples, n_features)
           >>> y_train = rng.randint(n_targets, size=n_samples)
           >>> rgs = FastKernelClassification(n_epoch=3, bandwidth=1,subsample_size=50)
           >>> rgs.fit(x_train, y_train)
           FastKernelClassification(bandwidth=1, bs='auto', coef0=1, degree=3,
                        gamma=None, kernel='gaussian', kernel_params=None, mem_gb=12,
                        n_components=1000, n_epoch=3, random_state=None,
                        subsample_size=50)

           >>> y_pred = rgs.predict(x_train)
           >>> loss = np.mean(np.square(y_train - y_pred))
       """

    def __init__(self, bs="auto", n_epoch=1, n_components=1000,
                 subsample_size="auto", mem_gb=12, kernel="gaussian",
                 bandwidth=5, gamma=None, degree=3, coef0=1,
                 kernel_params=None, random_state=None):
        self.bs = bs
        self.n_epoch = n_epoch
        self.n_components = n_components
        self.subsample_size = subsample_size
        self.mem_gb = mem_gb
        self.kernel = kernel
        self.bandwidth = bandwidth
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params
        self.random_state = random_state

    def fit(self, X, y):
        """ Train fast kernel classification model

        Parameters
        ----------
            X : {float, array}, shape = [n_samples, n_raw_feature]
                The raw input feature matrix.

            y : {float, array}, shape =[n_samples]
                The labels corresponding to the features of X.

        Returns
        -------
            self : returns an instance of self.
       """
        self.regressor_ = FastKernelRegression(
            bs=self.bs, n_epoch=self.n_epoch, n_components=self.n_components,
            subsample_size=self.subsample_size, mem_gb=self.mem_gb, kernel=self.kernel,
            bandwidth=self.bandwidth, gamma=self.gamma, degree=self.degree, coef0=self.coef0,
            kernel_params=self.kernel_params, random_state=self.random_state)
        X, y = check_X_y(X, y, multi_output=False, ensure_min_samples=3)
        check_classification_targets(y)
        self.classes_, ind = np.unique(y, return_inverse=True)

        loc = {}
        for ind, label in enumerate(self.classes_):
            loc[label] = ind

        class_matrix = np.zeros((y.shape[0], self.classes_.shape[0]))

        for ind, label in enumerate(y):
            class_matrix[ind][loc[label]] = 1

        self.regressor_.fit(X, class_matrix)

        return self

    def predict(self, X):
        """Predict using the kernel classification model

        Parameters
        ----------
            X: {float, array}, shape = [n_samples, n_features]
                Samples.

        Returns
        -------
            y : {float, array}, shape = [n_samples]
                Predicted labels.
        """
        check_is_fitted(self, ["regressor_"])
        Y = self.regressor_.predict(X)
        return self.classes_[np.argmax(Y, axis=1)]
