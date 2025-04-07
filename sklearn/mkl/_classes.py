"""Multiple Kernel Learning (MKL) classes."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Real

from ..base import (
    ClassifierMixin,
    OutlierMixin,
    RegressorMixin,
)
from ..utils._param_validation import Interval
from ._base import BaseMKL
from ._svm import SVC, SVR, OneClassSVM


class MKLC(ClassifierMixin, BaseMKL):
    """
    Multiple Kernel Learning for Classification.

    This classifier learns an optimal linear combination of multiple kernels and uses
    it within a support vector machine (SVM) for classification. After fitting, it can
    also be used as a transformer to apply the learned kernel combination to new data.

    Read more in the :ref:`User Guide <mkl>`.

    Parameters
    ----------
    kernels : list of {str, callable}, default=None
        List of kernel functions to combine. Each kernel can be either a string
        identifying a valid kernel name from
        `sklearn.metrics.pairwise.PAIRWISE_KERNEL_FUNCTIONS` or a callable.
        If a string is provided, the corresponding function from
        scikit-learn's kernel set is used. If a callable is provided, it must take
        two arrays `X` and `Y` as input (of shapes (n_samples_X, n_features) and
        (n_samples_Y, n_features)) and return a kernel matrix of shape
        (n_samples_X, n_samples_Y), as in scikit-learn's pairwise kernels.
        If None, it is assumed that `X` is a list of precomputed kernel matrices
        (see `fit` for more details).

    kernels_scopes : list of {"single", "all"}, default=None
        Scope of each kernel. Must have the same length as `kernels`.
        For each kernel, if the corresponding element is "single", the kernel is
        applied independently to each feature, resulting in `n_features` kernels being
        generated. If the element is "all", the kernel is applied to the full feature
        set. If None, all kernels use the "all" scope.

    kernels_param_grids : list of dict, default=None
        List of dictionaries specifying parameter grids for each kernel. Must have the
        same length as `kernels`. Keys must be strings matching parameter names, and
        values must be lists or arrays of values to test. If None, default parameter
        search is performed.

    C : float, default=1.0
        Regularization parameter of the internal SVM classifier.

    algo : {"average", "sum", "simple"}, default="simple"
        Algorithm used for learning kernel weights.

        - "sum"
            The algorithm computes a weighted sum of the kernel matrices, where each
            kernel has a fixed weight of 1.
        - "average"
            The algorithm computes a weighted sum of the kernel matrices, where each
            kernel has a fixed weight of 1 / M, with M being the number of kernels.
        - "simple"
            This algorithm uses the SimpleMKL approach, where the kernel weights are
            learned using a convex optimization procedure to optimize a linear
            combination of kernels. For more details, refer to the paper
            [Rakotomamonjy2008]_.

    svm_params : dict, default=None
        Additional parameters to pass to the internal `sklearn.svm.SVC`. Parameters
        `C`, `kernel`, and `random_state` are handled directly by MKLC and should not
        be specified here.

    precompute_kernels : bool, default=None
        Whether to precompute the kernel matrices before training. If set to None,
        kernel matrices are precomputed automatically if enough memory is available.

    tol : float, default=None
        Stopping tolerance for the MKL optimization. If None, a default is set depending
        on the algorithm and problem type :

        - For "average" and "sum" algorithms, it is unused.
        - For "simple" algorithm, it is set to 1e-2 for binary classification and
          1e-1 for multiclass classification.

    numeric_tol : float, default=1e-8
        Numerical threshold below which a kernel weight is considered zero.

    verbose : bool, default=False
        Controls verbosity of the optimization process.

    max_iter : int, default=200
        Maximum number of iterations allowed for the MKL optimization (only relevant
        for iterative algorithms like "simple").

    random_state : int or RandomState instance, default=None
        Controls the pseudo-random number generation for certain algorithms.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        The unique class labels.

    weights_ : ndarray of shape (n_kernels,)
        Weights assigned to each kernel after training.

    n_kernels_ : int
        Number of kernels used in the model.

    n_samples_in_ : int
        Number of samples seen during :term:`fit`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    n_iter_ : int
        Number of iterations run by the MKL algorithm (if applicable).

    See Also
    --------
    sklearn.svm.SVC : The SVM classifier used internally.
    sklearn.metrics.pairwise : Contains built-in kernel functions.

    References
    ----------
    .. [Rakotomamonjy2008] `Rakotomamonjy, A., Bach, F., Canu, S., & Grandvalet, Y.
        (2008). SimpleMKL. "Journal of Machine Learning Research", 9, 2491-2521.
        <http://www.jmlr.org/papers/volume9/rakotomamonjy08a/rakotomamonjy08a.pdf>`_

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.preprocessing import StandardScaler
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.mkl import MKLC
    >>> mkl = MKLC(
    ...     kernels=["linear", "rbf"],
    ...     kernels_scopes=["single", "all"],
    ...     kernels_param_grids=[{}, {"gamma": [0.1, 1.0]}],
    ...     algo="simple",
    ...     C=1.0,
    ... )
    >>> clf = make_pipeline(StandardScaler(), mkl)
    >>> clf.fit(X, y)
    Pipeline(steps=[('standardscaler', StandardScaler()),
                    ('mklc',
                     MKLC(kernels=['linear', 'rbf'],
                          kernels_param_grids=[{}, {'gamma': [0.1, 1.0]}],
                          kernels_scopes=['single', 'all']))])

    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    For more examples, see the :ref:`User Guide <mkl>`.
    """

    _parameter_constraints: dict = {
        **BaseMKL._parameter_constraints,
        "C": [Interval(Real, 0.0, None, closed="right")],
        "svm_params": [dict, None],
        # "algo": [StrOptions({"algo1", "algo2", ...})] # Algorithms supported
    }

    def __init__(
        self,
        *,
        kernels=None,  # None or list of functions or list of strings
        kernels_scopes=None,  # None or list of {"single", "all"}
        kernels_param_grids=None,  # None or list of (str, dict)
        C=1.0,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,  # If none, it tries to compute the kernels
        tol=None,  # DOC: auto depending on algo
        numeric_tol=1e-8,
        verbose=False,
        max_iter=200,  # Maybe put None if new iterative algorithms are implemented
        random_state=None,
    ):
        super().__init__(
            algo=algo,
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self.C = C
        self._warn_svm_params(svm_params, ["C", "random_state"])
        self.svm_params = svm_params

    def decision_function(self, X):
        """Evaluate the decision function for the samples in X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) \
            or (n_kernels, n_samples, n_samples)
            The input data.

        Returns
        -------
        X : ndarray of shape (n_samples, n_classes * (n_classes-1) / 2)
            Returns the decision function of the sample for each class
            in the model.
            If decision_function_shape='ovr', the shape is (n_samples,
            n_classes).

        Notes
        -----
        This method calls the `decision_function` of the internal fitted SVM on the
        combined kernel matrix.

        If decision_function_shape='ovo', the function values are proportional
        to the distance of the samples X to the separating hyperplane. If the
        exact distances are required, divide the function values by the norm of
        the weight vector (``coef_``). See also `this question
        <https://stats.stackexchange.com/questions/14876/
        interpreting-distance-from-hyperplane-in-svm>`_ for further details.
        If decision_function_shape='ovr', the decision function is a monotonic
        transformation of ovo decision function.
        """
        kernel = self.transform(X)
        return self._svm.decision_function(kernel)

    def _set_svm(self):
        self._svm = SVC(
            kernel="precomputed",
            C=self.C,
            random_state=self.random_state,
            **({} if self.svm_params is None else self.svm_params),
        )

    def _post_learning_processing(self):
        super()._post_learning_processing()
        self.classes_ = self._svm.classes_


class MKLR(RegressorMixin, BaseMKL):
    """
    Multiple Kernel Learning for Regression.

    This regressor learns an optimal linear combination of multiple kernels and uses
    it within a support vector regressor (SVR) for regression. After fitting, it can
    also be used as a transformer to apply the learned kernel combination to new data.

    Read more in the :ref:`User Guide <mkl>`.

    Parameters
    ----------
    kernels : list of {str, callable}, default=None
        List of kernel functions to combine. Each kernel can be either a string
        identifying a valid kernel name from
        `sklearn.metrics.pairwise.PAIRWISE_KERNEL_FUNCTIONS` or a callable.
        If a string is provided, the corresponding function from
        scikit-learn's kernel set is used. If a callable is provided, it must take
        two arrays `X` and `Y` as input (of shapes (n_samples_X, n_features) and
        (n_samples_Y, n_features)) and return a kernel matrix of shape
        (n_samples_X, n_samples_Y), as in scikit-learn's pairwise kernels.
        If None, it is assumed that `X` is a list of precomputed kernel matrices
        (see `fit` for more details).

    kernels_scopes : list of {"single", "all"}, default=None
        Scope of each kernel. Must have the same length as `kernels`.
        For each kernel, if the corresponding element is "single", the kernel is
        applied independently to each feature, resulting in `n_features` kernels being
        generated. If the element is "all", the kernel is applied to the full feature
        set. If None, all kernels use the "all" scope.

    kernels_param_grids : list of dict, default=None
        List of dictionaries specifying parameter grids for each kernel. Must have the
        same length as `kernels`. Keys must be strings matching parameter names, and
        values must be lists or arrays of values to test. If None, default parameter
        search is performed.

    C : float, default=1.0
        Regularization parameter of the internal SVM regressor.

    epsilon : float, default=0.1
        Epsilon parameter of the internal SVM regressor. Specifies the epsilon-tube
        within which no penalty is associated in the training loss function.

    algo : {"average", "sum", "simple"}, default="simple"
        Algorithm used for learning kernel weights.

        - "sum"
            The algorithm computes a weighted sum of the kernel matrices, where each
            kernel has a fixed weight of 1.
        - "average"
            The algorithm computes a weighted sum of the kernel matrices, where each
            kernel has a fixed weight of 1 / M, with M being the number of kernels.
        - "simple"
            This algorithm uses the SimpleMKL approach, where the kernel weights are
            learned using a convex optimization procedure to optimize a linear
            combination of kernels. For more details, refer to the paper
            [Rakotomamonjy2008]_.

    svm_params : dict, default=None
        Additional parameters to pass to the internal `sklearn.svm.SVR`. Parameters
        `C`, `epsilon`, `kernel`, and `random_state` are handled directly by MKLR and
        should not be specified here.

    precompute_kernels : bool, default=None
        Whether to precompute the kernel matrices before training. If set to None,
        kernel matrices are precomputed automatically if enough memory is available.

    tol : float, default=None
        Stopping tolerance for the MKL optimization. If None, a default is set depending
        on the algorithm.

        - For "average" and "sum" algorithms, it is unused.
        - For "simple" algorithm, it is set to 1e-2.

    numeric_tol : float, default=1e-8
        Numerical threshold below which a kernel weight is considered zero.

    verbose : bool, default=False
        Controls verbosity of the optimization process.

    max_iter : int, default=200
        Maximum number of iterations allowed for the MKL optimization (only relevant
        for iterative algorithms like "simple").

    random_state : int or RandomState instance, default=None
        Controls the pseudo-random number generation for certain algorithms.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    weights_ : ndarray of shape (n_kernels,)
        Weights assigned to each kernel after training.

    n_kernels_ : int
        Number of kernels used in the model.

    n_samples_in_ : int
        Number of samples seen during :term:`fit`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    n_iter_ : int
        Number of iterations run by the MKL algorithm (if applicable).

    See Also
    --------
    sklearn.svm.SVR : The SVR regressor used internally.
    sklearn.metrics.pairwise : Contains built-in kernel functions.

    References
    ----------
    .. [Rakotomamonjy2008] `Rakotomamonjy, A., Bach, F., Canu, S., & Grandvalet, Y.
        (2008). SimpleMKL. "Journal of Machine Learning Research", 9, 2491-2521.
        <http://www.jmlr.org/papers/volume9/rakotomamonjy08a/rakotomamonjy08a.pdf>`_

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.preprocessing import StandardScaler
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([-2.0, -3.0, 3.0, 4.0])
    >>> from sklearn.mkl import MKLR
    >>> mkl = MKLR(
    ...     kernels=["linear", "rbf"],
    ...     kernels_scopes=["single", "all"],
    ...     kernels_param_grids=[{}, {"gamma": [0.1, 1.0]}],
    ...     algo="simple",
    ...     C=1.0,
    ...     epsilon=0.1,
    ... )
    >>> reg = make_pipeline(StandardScaler(), mkl)
    >>> reg.fit(X, y)
    Pipeline(steps=[('standardscaler', StandardScaler()),
                ('mklr',
                 MKLR(kernels=['linear', 'rbf'],
                      kernels_param_grids=[{}, {'gamma': [0.1, 1.0]}],
                      kernels_scopes=['single', 'all']))])

    >>> print(reg.predict([[0, 0]]))
    [0.5]

    For more examples, see the :ref:`User Guide <mkl>`.
    """

    _parameter_constraints: dict = {
        **BaseMKL._parameter_constraints,
        "C": [Interval(Real, 0.0, None, closed="right")],
        "epsilon": [Interval(Real, 0.0, None, closed="left")],
        "svm_params": [dict, None],
        # "algo": [StrOptions({"algo1", "algo2", ...})] # Algorithms supported
    }

    def __init__(
        self,
        *,
        kernels=None,
        kernels_scopes=None,
        kernels_param_grids=None,
        C=1.0,
        epsilon=0.1,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,
        tol=None,
        numeric_tol=1e-8,
        verbose=False,
        max_iter=200,  # Maybe put None if new iterative algorithms are implemented
        random_state=None,
    ):
        super().__init__(
            algo=algo,
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self.C = C
        self.epsilon = epsilon
        self._warn_svm_params(svm_params, ["C", "epsilon"])
        self.svm_params = svm_params

    def _set_svm(self):
        self._svm = SVR(
            kernel="precomputed",
            C=self.C,
            epsilon=self.epsilon,
            **({} if self.svm_params is None else self.svm_params),
        )


class OneClassMKL(OutlierMixin, BaseMKL):
    """
    Multiple Kernel Learning for Outlier Detection.

    This model learns an optimal linear combination of multiple kernels and uses
    it within a One-Class SVM for outlier or novelty detection. After fitting, it can
    also be used as a transformer to apply the learned kernel combination to new data.

    Read more in the :ref:`User Guide <mkl>`.

    Parameters
    ----------
    kernels : list of {str, callable}, default=None
        List of kernel functions to combine. Each kernel can be either a string
        identifying a valid kernel name from
        `sklearn.metrics.pairwise.PAIRWISE_KERNEL_FUNCTIONS` or a callable.
        If a string is provided, the corresponding function from
        scikit-learn's kernel set is used. If a callable is provided, it must take
        two arrays `X` and `Y` as input (of shapes (n_samples_X, n_features) and
        (n_samples_Y, n_features)) and return a kernel matrix of shape
        (n_samples_X, n_samples_Y), as in scikit-learn's pairwise kernels.
        If None, it is assumed that `X` is a list of precomputed kernel matrices
        (see `fit` for more details).

    kernels_scopes : list of {"single", "all"}, default=None
        Scope of each kernel. Must have the same length as `kernels`.
        For each kernel, if the corresponding element is "single", the kernel is
        applied independently to each feature, resulting in `n_features` kernels being
        generated. If the element is "all", the kernel is applied to the full feature
        set. If None, all kernels use the "all" scope.

    kernels_param_grids : list of dict, default=None
        List of dictionaries specifying parameter grids for each kernel. Must have the
        same length as `kernels`. Keys must be strings matching parameter names, and
        values must be lists or arrays of values to test. If None, default parameter
        search is performed.

    nu : float, default=0.5
        Anomaly proportion parameter of the internal One-Class SVM.
        Specifies an upper bound on the fraction of training errors and a lower bound
        on the fraction of support vectors. Must be in the interval (0, 1].

    algo : {"average", "sum", "simple"}, default="simple"
        Algorithm used for learning kernel weights.

        - "sum"
            The algorithm computes a weighted sum of the kernel matrices, where each
            kernel has a fixed weight of 1.
        - "average"
            The algorithm computes a weighted sum of the kernel matrices, where each
            kernel has a fixed weight of 1 / M, with M being the number of kernels.
        - "simple"
            This algorithm uses the SimpleMKL approach, where the kernel weights are
            learned using a convex optimization procedure to optimize a linear
            combination of kernels. For more details, refer to the paper
            [Rakotomamonjy2008]_.

    svm_params : dict, default=None
        Additional parameters to pass to the internal `sklearn.svm.OneClassSVM`.
        Parameters `kernel`, `nu`, and `random_state` are handled directly by
        OneClassMKL and should not be specified here.

    precompute_kernels : bool, default=None
        Whether to precompute the kernel matrices before training. If set to None,
        kernel matrices are precomputed automatically if enough memory is available.

    tol : float, default=None
        Stopping tolerance for the MKL optimization. If None, a default is set depending
        on the algorithm.

        - For "average" and "sum" algorithms, it is unused.
        - For "simple" algorithm, it is set to 1e-2.

    numeric_tol : float, default=1e-8
        Numerical threshold below which a kernel weight is considered zero.

    verbose : bool, default=False
        Controls verbosity of the optimization process.

    max_iter : int, default=200
        Maximum number of iterations allowed for the MKL optimization (only relevant
        for iterative algorithms like "simple").

    random_state : int or RandomState instance, default=None
        Controls the pseudo-random number generation for certain algorithms.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    weights_ : ndarray of shape (n_kernels,)
        Weights assigned to each kernel after training.

    n_kernels_ : int
        Number of kernels used in the model.

    n_samples_in_ : int
        Number of samples seen during :term:`fit`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    n_iter_ : int
        Number of iterations run by the MKL algorithm (if applicable).

    offset_ : float
        Offset used to define the decision function from the raw scores.
        We have the relation: decision_function = score_samples - `offset_`.

    See Also
    --------
    sklearn.svm.OneClassSVM : The One-Class SVM used internally.
    sklearn.metrics.pairwise : Contains built-in kernel functions.

    References
    ----------
    .. [Rakotomamonjy2008] `Rakotomamonjy, A., Bach, F., Canu, S., & Grandvalet, Y.
        (2008). SimpleMKL. "Journal of Machine Learning Research", 9, 2491-2521.
        <http://www.jmlr.org/papers/volume9/rakotomamonjy08a/rakotomamonjy08a.pdf>`_

    Examples
    --------
    >>> from sklearn.mkl import OneClassMKL
    >>> X = [[0], [0.44], [0.45], [0.46], [1]]
    >>> mkl = OneClassMKL(
    ...     kernels=["linear", "rbf"],
    ...     kernels_scopes=["single", "all"],
    ...     kernels_param_grids=[{}, {"gamma": [1.0, 10.0]}],
    ...     algo="simple",
    ...     nu=0.5,
    ... )
    >>> mkl.fit(X)
    >>> mkl.predict(X)
    array([-1,  1,  1,  1, -1])
    >>> mkl.score_samples(X)
    array([1.77987316, 2.05479873, 2.05560497, 2.05615569, 1.73328509])

    For more examples, see the :ref:`User Guide <mkl>`.
    """

    _parameter_constraints: dict = {
        **BaseMKL._parameter_constraints,
        "nu": [Interval(Real, 0.0, 1.0, closed="right")],
        "svm_params": [dict, None],
        # "algo": [StrOptions({"algo1", "algo2", ...})] # Algorithms supported
    }

    def __init__(
        self,
        *,
        kernels=None,
        kernels_scopes=None,
        kernels_param_grids=None,
        nu=0.5,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,
        tol=None,
        numeric_tol=1e-8,
        verbose=False,
        max_iter=200,  # Maybe put None if new iterative algorithms are implemented
        random_state=None,
    ):
        super().__init__(
            algo=algo,
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self.nu = nu
        self._warn_svm_params(svm_params, ["nu"])
        self.svm_params = svm_params

    def decision_function(self, X):
        """Signed distance to the separating hyperplane.

        Signed distance is positive for an inlier and negative for an outlier.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) \
            or (n_kernels, n_samples, n_samples)
            The input data.

        Returns
        -------
        dec : ndarray of shape (n_samples,)
            Returns the decision function of the samples.

        Notes
        -----
        This method calls the `decision_function` of the internal fitted SVM on the
        combined kernel matrix.
        """
        kernel = self.transform(X)
        return self._svm.decision_function(kernel)

    def score_samples(self, X):
        """Raw scoring function of the samples.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) \
            or (n_kernels, n_samples, n_samples)
            The input data.

        Returns
        -------
        score_samples : ndarray of shape (n_samples,)
            Returns the (unshifted) scoring function of the samples.

        Notes
        -----
        This method calls the `score_samples` of the internal fitted SVM on the
        combined kernel matrix.
        """
        kernel = self.transform(X)
        return self._svm.score_samples(kernel)

    def _set_svm(self):
        svm_params = {} if self.svm_params is None else self.svm_params
        self._svm = OneClassSVM(
            kernel="precomputed",
            nu=self.nu,
            **svm_params,
        )

    def _post_learning_processing(self):
        super()._post_learning_processing()
        self.offset_ = self._svm.offset_
