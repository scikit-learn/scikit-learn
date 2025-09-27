# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Polynomial Chaos Expansion"""

from warnings import warn

import numpy as np
from scipy.stats import uniform

from sklearn.base import BaseEstimator, RegressorMixin, _fit_context, clone
from sklearn.exceptions import DataConversionWarning
from sklearn.linear_model._base import LinearModel, LinearRegression
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.polynomial_chaos._adaptive import BasisIncrementStrategy
from sklearn.preprocessing._orthogonal import OrthogonalPolynomialFeatures
from sklearn.utils._orthogonal_polynomial import Polynomial
from sklearn.utils._param_validation import (
    HasMethods,
    Integral,
    Interval,
    Iterable,
    StrOptions,
)
from sklearn.utils.validation import (
    _get_feature_names,
    check_is_fitted,
    check_X_y,
    column_or_1d,
    validate_data,
)


class PolynomialChaosExpansion(RegressorMixin, BaseEstimator):
    """Polynomial Chaos Expansion.

    In addition to the standard scikit-learn estimator API, this estimator:

        * allows the computation of statistics (mean and variance) of the
          output, and
        * enables global sensitivity analysis through the calculation of main-
          and total-effect Sobol sensitivity indices.

    Read more in the :ref:`User Guide <polynomial_chaos>`.

    Parameters
    ----------
    distribution : scpiy.stats distribution or tuple (distribution_1, \
        distribution_2, ...), default=None
        Distribution of the input parameter(s) of the Polynomial Chaos
        expansion. If a single distribution is given, it specifies the
        distributions for all input features. If a tuple (`distribution_1`,
        `distribution_2`, ...) is passed, then `distribution_1` is the
        distribution of the first feature, `distribution_2` is the distribution
        of the second feature, and so on. Note that when a tuple is passed, its
        length must be consistent with the the number of input features and the
        number of columns in `multiindices` (when the latter is provided). If
        `distribution = None`, we will assume a uniform distribution where the
        lower and upper bounds are extracted from the input features.

    degree : int, default=2
        The maximum degree of the Polynomial Chaos expansion.

    truncation : {'full_tensor', 'total_degree', 'hyperbolic_cross', \
        'Zaremba_cross'}, default='total_degree'
        The truncation rule that should be used to determine the shape of the
        basis of the Polynomial Chaos expansion. The default is
        `'total_degree'`.

    weights : array_like, default=None
        Optional weights that can be used to select certain basis terms where
        higher-degree polynomials should be used in the basis. A larger value
        for the weight of a certain feature indicates that a higher-degree
        polynomial is used. The weights must be all positive. When `weights =
        None`, an unweighted multiindex set will be used. The default is
        `None`.

    estimator : LinearModel, default=LinearRegression(fit_intercept=False)
        The :class:`~sklearn.linear_model.LinearModel` used to solve for the
        coefficients of the Polynomial Chaos expansion. This should be another
        `LinearModel` that has a :term:`fit` method. Make sure to set
        `fit_intercept = False`.

    feature_selector : SelectorMixin, default=None
        The meta-transformer used to select basis terms after fitting. This
        is useful when using sparsity promoting linear estimators, such as
        :class:`~sklearn.linear_model.LassoCV`. The feature selector can be
        used to prune the basis terms with near-zero coefficients to speed-up
        and reduce memory usage at prediction time. If
        `feature_selector = None`, all features will be retained.

    multiindices : ndarray of shape \
        (n_output_features_, n_features_in_), default=None
        The combination of `degree`, `truncation` and `weights` provides a
        flexible way to define the Polynomial Chaos basis. To allow for even
        more fine-grained control, this optional argument allows to specify an
        arbitrary set of basis terms that will be used instead. When this
        argument is provided, it supersedes the values in `degree`,
        `truncation` and `weights`. If `multiindices = None`, then the
        multiindex set shape given in `truncation` will be used. The default is
        `None`.

    scale_outputs : bool, default=True
        When true, the output features are scaled to have zero mean and unit
        variance before the Polynomial Chaos expansion is fitted. In certain
        cases, this option can improve the numerics.

    Attributes
    ----------
    coef_ : array-like of length (n_terms,) if n_features_in_ == 1, \
        array-like of shape (n_output_features_, n_terms_) otherwise
        The coefficients of this Polynomial Chaos expansion.

    distributions_ : array-like of length (n_terms,)
        The input distributions seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (n_features_in_,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    multiindices_ : array-like of shape (n_terms, n_features_in_)
        An array with the combinations of input features and polynomial degrees
        that constitute the Polynomial Chaos basis. Every row in this array
        contains a single multiindex.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    n_outputs_ : int
        Number of outputs seen during :term:`fit`.

    norms_ : array-like of shape (n_terms,)
        The norm of each polynomial basis terms.

    output_mean_ : array-like of shape (n_outputs_,)
        The mean of the output (when `scale_outputs = True`).

    output_std_ : array-like of shape (n_outputs_,)
        The standard deviation of the output (when `scale_outputs = True`).

    pipeline_ : sklearn.pipeline.Pipeline
        The pipeline constructed during :term:`fit`. The pipeline has the
        following steps:
        - `"basis_transformer"`: an \
            :class:`~sklearn.preprocessing.OrthogonalPolynomialFeatures` \
            transformer
        - `"feature_selector"`: a \
            :ref:`User Guide <univariate_feature_selection>` \
            meta-transformer such as \
            :class:`~sklearn.feature_selection.SelectFromModel` (only when \
            `feature_selector` is not `None`)
        - `"estimator"`: a :class:`~sklearn.linear_model.LinearModel`

    polynomials_ : array-like of shape (n_features_in_,)
        The orthogonal polynomial associated with each input feature.

    strategy_ : skleanr.polynomial_chaos.BasisIncrementStrategy
        The adaptive basis growth strategy used during :term:`fit`.

    See Also
    --------
    :class:`~sklearn.preprocessing.OrthogonalPolynomialFeatures`: Transformer
        that maps features into orthogonal polynomial features.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import uniform
    >>> from sklearn.polynomial_chaos import PolynomialChaosExpansion
    >>> X = np.linspace(0, 1).reshape(-1, 2)
    >>> y = np.prod((3 * X**2 + 1) / 2, axis=1)
    >>> pce = PolynomialChaosExpansion(uniform(), degree=4)
    >>> _ = pce.fit(X, y)
    >>> sens = pce.total_sens()
    """

    _parameter_constraints: dict = {
        "distribution": [HasMethods("dist"), "array-like", None],
        "degree": [Interval(Integral, 0, None, closed="left")],
        "truncation": [
            StrOptions(
                {
                    "full_tensor",
                    "total_degree",
                    "hyperbolic_cross",
                    "Zaremba_cross",
                }
            )
        ],
        "weights": ["array-like", None],
        "estimator": [HasMethods("fit"), None],
        "feature_selector": [HasMethods("fit"), None],
        "multiindices": ["array-like", None],
        "scale_outputs": [bool],
    }

    def __init__(
        self,
        distribution=None,
        degree=2,
        truncation="total_degree",
        weights=None,
        estimator=None,
        feature_selector=None,
        multiindices=None,
        scale_outputs=True,
    ):
        self.distribution = distribution
        self.degree = degree
        self.truncation = truncation
        self.weights = weights
        self.estimator = estimator
        self.feature_selector = feature_selector
        self.multiindices = multiindices
        self.scale_outputs = scale_outputs

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, strategy="gerstner_griebel", n_iter=1):
        """Fit Polynomial Chaos Expansion model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Feature vectors or other representations of training data.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values.

        strategy : {'gerstner_griebel'}, default='gerstner_griebel'
            Basis increment strategy. After the first :term:`fit`, we can use a
            basis increment strategy to add more polynomials to the basis where
            further refinement is needed. Only the adaptive basis growth
            procedure from [Gerstner2003]_ is implemented at the moment. In
            this method, we use the variance contribution of each polynomial
            term to select the multiindices that are suitable for further
            refinement.

        n_iter : int
            Maximum number of iterations with the adaptive algorithm. The
            default is `1`, meaning that the basis will not be incremented.

        Returns
        -------
        self : object
            PolynomialChaosExpansion class instance.
        """
        # Check input features
        self.feature_names_in_ = _get_feature_names(X)
        X, y = check_X_y(X, y, multi_output=True, y_numeric=True)
        self.n_features_in_ = X.shape[1]

        # Cannot perform regression with less than 2 data points
        if len(y) < 2:
            raise ValueError(f"expected more than 1 sample, got {len(y)}")

        # Check distributions
        if self.distribution is None:
            data_min = np.amin(X, axis=0)
            data_max = np.amax(X, axis=0)
            self.distributions_ = [
                uniform(loc=a, scale=b - a) for a, b in zip(data_min, data_max)
            ]
        elif hasattr(self.distribution, "dist"):
            self.distributions_ = [self.distribution] * self.n_features_in_
        elif isinstance(self.distribution, Iterable):
            self.distributions_ = list(self.distribution)
            if not len(self.distributions_) == self.n_features_in_:
                raise ValueError(
                    "the number of distributions does not match the number "
                    f"of input features, got {len(self.distributions_)} but "
                    f"expected {self.n_features_in_}"
                )
            for j, distribution in enumerate(self.distributions_):
                if not hasattr(distribution, "dist"):
                    raise ValueError(
                        "distributions must be all of type 'scipy.stats' "
                        "frozen distribution, but the distribution at index "
                        f"{j} has type '{type(distribution)}'"
                    )
        else:
            raise ValueError(
                "distribution must be a 'scipy.stats' frozen distribution "
                f"or a tuple/list, got '{type(self.distribution)}'"
            )

        # Check estimator
        if self.estimator is None:
            estimator = LinearRegression(fit_intercept=False)
        else:
            if not isinstance(self.estimator, BaseEstimator):
                raise ValueError(
                    "estimator must be an instance of 'sklearn.BaseEstimator', got"
                    f" '{type(self.estimator)}'"
                )
            if isinstance(self.estimator, LinearModel) or isinstance(
                self.estimator, MultiOutputRegressor
            ):
                if (
                    isinstance(self.estimator, LinearModel)
                    and self.estimator.fit_intercept
                ) or (
                    isinstance(self.estimator, MultiOutputRegressor)
                    and self.estimator.estimator.fit_intercept
                ):
                    raise ValueError(
                        "make sure to set 'fit_intercept=False' in estimator"
                    )

            else:
                warn(
                    (
                        "cannot explicitly check the value of 'fit_intercept' in"
                        " estimator, please ensure 'fit_intercept=False'"
                    ),
                    UserWarning,
                )
            estimator = clone(self.estimator)

        # Check feature_selector
        if self.feature_selector is not None:
            if not (
                hasattr(self.feature_selector, "fit")
                and hasattr(self.feature_selector, "transform")
            ):
                raise ValueError(
                    "feature_selector must be a transformer with fit and transform "
                    f"methods, got '{type(self.feature_selector)}'"
                )
            feature_selector = clone(self.feature_selector)
        else:
            feature_selector = None

        # Check outputs
        y = np.atleast_1d(y)
        if y.ndim == 2 and y.shape[1] == 1:
            warn(
                (
                    "A column-vector y was passed when a 1d array was"
                    " expected. Please change the shape of y to "
                    "(n_samples, ), for example using ravel()."
                ),
                DataConversionWarning,
                stacklevel=2,
            )
        self.n_outputs_ = 1 if y.ndim == 1 else y.shape[1]

        # Get orthogonal polynomials for each distribution
        self.polynomials_ = list()
        for distribution in self.distributions_:
            self.polynomials_.append(Polynomial.from_distribution(distribution))

        # Scale input features
        X_scaled = np.zeros_like(X)
        for j, (distribution, polynomial) in enumerate(
            zip(self.distributions_, self.polynomials_)
        ):
            X_scaled[:, j] = polynomial.scale_features_from_distribution(
                X[:, j], distribution
            )

        # Scale outputs
        self.output_mean_ = (
            np.mean(y, axis=0) if self.scale_outputs else np.full(self.n_outputs_, 0)
        )
        self.output_std_ = (
            np.std(y, axis=0) if self.scale_outputs else np.full(self.n_outputs_, 1)
        )
        if self.n_outputs_ > 1:
            self.output_std_[self.output_std_ == 0] = 1  # avoid / 0 error
        y = y - self.output_mean_
        y = y / self.output_std_

        self.multiindices_ = self.multiindices
        self.strategy_ = BasisIncrementStrategy.from_string(strategy)

        # Adaptive basis growth
        for iter in range(n_iter):
            # Create orthogonal polynomial basis transformer
            basis = OrthogonalPolynomialFeatures(
                self.degree,
                [str(polynomial) for polynomial in self.polynomials_],
                normalize=True,
                truncation=self.truncation,
                weights=self.weights,
                multiindices=self.multiindices_,
            )

            # Build pipeline
            steps = [("basis_transformer", basis)]
            if feature_selector is not None:
                steps.append(("feature_selector", feature_selector))
            steps.append(("estimator", estimator))

            self.pipeline_ = Pipeline(steps)

            # Fit pipeline
            self.pipeline_.fit(X_scaled, y)

            # Convenient access to multiindices and norms
            self.multiindices_ = self.pipeline_["basis_transformer"].multiindices_

            if feature_selector is not None:
                # After fitting, feature_selector has support mask
                support = self.pipeline_["feature_selector"].get_support()
                self.multiindices_ = self.multiindices_[support]

            # Convenient access to coefficients
            if hasattr(self.pipeline_["estimator"], "coef_"):
                self.coef_ = self.pipeline_["estimator"].coef_
            elif hasattr(self.pipeline_["estimator"], "estimators_"):
                self.coef_ = np.vstack(
                    [
                        estimator.coef_
                        for estimator in self.pipeline_["estimator"].estimators_
                    ]
                )
            else:
                raise ValueError(
                    "unable to access coefficients from estimator, please make sure"
                    " the estimator is a 'LinearModel' with a 'coef_' attribute, or a"
                    " 'MultiOutputRegressor' with an estimator that has a 'coef_'"
                    " attribute"
                )

            # Do not update multiindices in last iteration
            if iter < n_iter - 1:
                self.multiindices_ = self.strategy_.propose(self)

        # By convention, fit returns itself
        return self

    def predict(self, X):
        """Predict using this Polynomial Chaos model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Points where the Polynomial Chaos expansion should be evaluated.

        Returns
        -------
        y : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            Polynomial Chaos predictions in query points.
        """
        # Check if this estimator is fitted
        check_is_fitted(self)

        # Check input features
        X = validate_data(self, X, reset=False, ensure_2d=True)

        # Scale features
        X_scaled = np.zeros_like(X)
        for j, (distribution, polynomial) in enumerate(
            zip(self.distributions_, self.polynomials_)
        ):
            X_scaled[:, j] = polynomial.scale_features_from_distribution(
                X[:, j], distribution
            )

        y_fit = self.pipeline_.predict(X_scaled)
        return np.squeeze(y_fit * self.output_std_ + self.output_mean_)

    def joint_sens(self, *features):
        r"""Return the joint sensivitivity indices.

        Given a Polynomial Chaos expansion

        .. math::
            y = \sum_u c_u \Psi_u(\xi)

        where :math:`c_u` are the coefficients and :math:`\Psi_u` are the
        orthogonal polynomials in the expansion, the Polynomial Chaos-based
        Sobol indices can be computed as

        .. math::
            S_{\{i_0, i_1, \ldots, i_s\}} = \sum_{u \in \mathcal{I}} c_u^2
            \mathbb{E}[\Psi_u^2] / \mathbb{V}[y]

        where :math:`\mathcal{I}` is the set of multiindices for which all
        components at indices :math:`i_0, i_1, \ldots, i_s` are positive, and
        all other components are :math:`0`.

        When only one input feature is provided, this corresponds to the main
        Sobol sensitivity index. When two or more input features are provided,
        this corresponds to the joint Sobol sensitivity index.

        When the features are named, this method accepts any combination of
        input feature names. If the features are unnamed, this method accepts
        the index of the input features. Duplicate features or feature names
        are not accepted.

        Parameters
        ----------
        *features : int or str
            Input features to compute the joint sensitivity index of. If a
            single `int` or `str` argument is given, we return the main
            sensitivity index of the input feature with the given index or
            name. If a `tuple` of `int`\ s or `str`\ s is given, we return the
            joint sensitivity index of the input features with the given
            indices or names. In the latter case, the input features must be
            all `int`\ s or all `str`\ s. Feature indices or names must be
            unique.

        Returns
        -------
        sensitivity_index : (n_outputs,)
            The joint sensitivity index for the given set of input features.
        """
        check_is_fitted(self)

        # Input checking
        dtype = type(features[0])
        for arg in features:
            if not isinstance(arg, dtype):
                raise ValueError("inputs must be all string or all int")
        if dtype == str:
            if self.feature_names_in_ is None:
                raise ValueError("feature names have not been set")
            idcs = np.zeros(len(features), dtype=int)
            for j, arg in enumerate(features):
                if arg not in self.feature_names_in_:
                    raise ValueError(f"feature '{arg}' not found")
                idcs[j] = np.where(self.feature_names_in_ == arg)[0][0]
        else:
            idcs = column_or_1d(features)
        if len(np.unique(idcs)) != len(idcs):
            raise ValueError("features must be unique")
        if len(features) > self.n_features_in_ or np.any(idcs > self.n_features_in_):
            raise ValueError(f"this model has only {self.n_features_in_} features")

        # Actually compute the joint sensitivity index
        joint_sens = np.zeros(self.n_outputs_)
        for j, index in enumerate(self.multiindices_):
            if np.sum(index != 0) == len(idcs) and all(i > 0 for i in index[idcs]):
                joint_sens += self.output_std_**2 * np.atleast_2d(self.coef_)[:, j] ** 2

        return np.divide(joint_sens.T, self.var()).T

    def main_sens(self):
        r"""Return the main sensitivity indices.

        Given a Polynomial Chaos expansion

        .. math::
            y = \sum_u c_u \Psi_u(\xi)

        where :math:`c_u` are the coefficients and :math:`\Psi_u` are the
        orthogonal polynomials in the expansion, the main Polynomial
        Chaos-based Sobol indices can be computed as

        .. math::
            S_{j} = \sum_{u \in \mathcal{I}_j} c_u^2
            \mathbb{E}[\Psi_u^2] / \mathbb{V}[y]

        where :math:`\mathcal{I}_j` is the set of multiindices for which the
        :math:`j`\ th coordinate (and only that coordinate) is larger than 0.
        This sensitivity index expresses the effect on the output variance of
        varying only the :math:`j`\ th parameter .

        Returns
        -------
        sensitivity_indices : array-like (n_outputs, n_features)
            The main-effectt Sobol sensitivity indices for all input features.
        """
        check_is_fitted(self)
        return np.vstack([self.joint_sens(i) for i in range(self.n_features_in_)]).T

    def total_sens(self):
        r"""Return the total sensitivity indices.

        Given a Polynomial Chaos expansion

        .. math::
            y = \sum_u c_u \Psi_u(\xi)

        where :math:`c_u` are the coefficients and :math:`\Psi_u` are the
        orthogonal polynomials in the expansion, the main Polynomial
        Chaos-based Sobol indices can be computed as

        .. math::
            S_{j}^T = \sum_{u \in \mathcal{I}_j^T} c_u^2
            \mathbb{E}[\Psi_u^2] / \mathbb{V}[y]

        where :math:`\mathcal{I}_j^T` is the set of multiindices for which the
        :math:`j`\ th coordinate, amongst others, is larger than 0. This
        sensitivity index expresses the effect on the output variance of
        varying the :math:`j`\ th parameter, including all interaction terms.

        Returns
        -------
        sensitivity_indices : array-like (n_outputs, n_features)
            The total-effectt Sobol sensitivity indices for all input features.
        """
        check_is_fitted(self)
        total_sens = np.zeros((self.n_outputs_, self.n_features_in_))
        for i in range(self.n_features_in_):
            for j, index in enumerate(self.multiindices_):
                if index[i] == 0:
                    continue
                total_sens[:, i] += (
                    self.output_std_**2 * np.atleast_2d(self.coef_)[:, j] ** 2
                )
        return np.divide(total_sens.T, self.var()).T

    def mean(self):
        r"""Return the approximation for the mean of the model output.

        Given a Polynomial Chaos expansion

        .. math::
            y = \sum_u c_u \Psi_u(\xi)

        where :math:`c_u` are the coefficients and :math:`\Psi_u` are the
        orthogonal polynomials in the expansion, the mean of the response
        :math:`y` is

        .. math::
            \mathbb{E}[y] = c_0

        Returns
        -------
        mean : array-like (n_outputs,)
            The approximation for the mean of the model output.
        """
        check_is_fitted(self)
        for j, index in enumerate(self.multiindices_):
            if np.all(index == 0):
                return (
                    np.atleast_2d(self.coef_)[:, j] * self.output_std_
                    + self.output_mean_
                )
        return 0

    def var(self):
        r"""Return the approximation for the variance of the model output.

        Given a Polynomial Chaos approximation

        .. math::
            y = \sum_u c_u \Psi_u(\xi)

        where :math:`c_u` are the coefficients and :math:`\Psi_u` are the
        orthogonal polynomials in the expansion, the variance of the response
        :math:`y` is

        .. math::
            \mathbb{V}[y] = \sum_{u > 0} c_u^2 \mathbb{E}[\Psi_u^2]

        Returns
        -------
        var : array-like (n_outputs,)
            The approximation for the variance of the model output.
        """
        check_is_fitted(self)
        var = 0
        for j, index in enumerate(self.multiindices_):
            if np.any(index > 0):
                var += self.output_std_**2 * np.atleast_2d(self.coef_)[:, j] ** 2
        return var
