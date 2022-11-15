# flake8: noqa
"""
=======================================
Release Highlights for scikit-learn 1.0
=======================================

.. currentmodule:: sklearn

We are very pleased to announce the release of scikit-learn 1.0! The library
has been stable for quite some time, releasing version 1.0 is recognizing that
and signalling it to our users. This release does not include any breaking
changes apart from the usual two-release deprecation cycle. For the future, we
do our best to keep this pattern.

This release includes some new key features as well as many improvements and
bug fixes. We detail below a few of the major features of this release. **For
an exhaustive list of all the changes**, please refer to the :ref:`release
notes <changes_1_0>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

##############################################################################
# Keyword and positional arguments
# ---------------------------------------------------------
# The scikit-learn API exposes many functions and methods which have many input
# parameters. For example, before this release, one could instantiate a
# :class:`~ensemble.HistGradientBoostingRegressor` as::
#
#         HistGradientBoostingRegressor("squared_error", 0.1, 100, 31, None,
#             20, 0.0, 255, None, None, False, "auto", "loss", 0.1, 10, 1e-7,
#             0, None)
#
# Understanding the above code requires the reader to go to the API
# documentation and to check each and every parameter for its position and
# its meaning. To improve the readability of code written based on scikit-learn,
# now users have to provide most parameters with their names, as keyword
# arguments, instead of positional arguments. For example, the above code would
# be::
#
#     HistGradientBoostingRegressor(
#         loss="squared_error",
#         learning_rate=0.1,
#         max_iter=100,
#         max_leaf_nodes=31,
#         max_depth=None,
#         min_samples_leaf=20,
#         l2_regularization=0.0,
#         max_bins=255,
#         categorical_features=None,
#         monotonic_cst=None,
#         warm_start=False,
#         early_stopping="auto",
#         scoring="loss",
#         validation_fraction=0.1,
#         n_iter_no_change=10,
#         tol=1e-7,
#         verbose=0,
#         random_state=None,
#     )
#
# which is much more readable. Positional arguments have been deprecated since
# version 0.23 and will now raise a ``TypeError``. A limited number of
# positional arguments are still allowed in some cases, for example in
# :class:`~decomposition.PCA`, where ``PCA(10)`` is still allowed, but ``PCA(10,
# False)`` is not allowed.

##############################################################################
# Spline Transformers
# ---------------------------------------------------------
# One way to add nonlinear terms to a dataset's feature set is to generate
# spline basis functions for continuous/numerical features with the new
# :class:`~preprocessing.SplineTransformer`. Splines are piecewise polynomials,
# parametrized by their polynomial degree and the positions of the knots. The
# :class:`~preprocessing.SplineTransformer` implements a B-spline basis.
#
# .. figure:: ../linear_model/images/sphx_glr_plot_polynomial_interpolation_001.png
#   :target: ../linear_model/plot_polynomial_interpolation.html
#   :align: center
#
# The following code shows splines in action, for more information, please
# refer to the :ref:`User Guide <spline_transformer>`.

import numpy as np
from sklearn.preprocessing import SplineTransformer

X = np.arange(5).reshape(5, 1)
spline = SplineTransformer(degree=2, n_knots=3)
spline.fit_transform(X)


##############################################################################
# Quantile Regressor
# --------------------------------------------------------------------------
# Quantile regression estimates the median or other quantiles of :math:`y`
# conditional on :math:`X`, while ordinary least squares (OLS) estimates the
# conditional mean.
#
# As a linear model, the new :class:`~linear_model.QuantileRegressor` gives
# linear predictions :math:`\hat{y}(w, X) = Xw` for the :math:`q`-th quantile,
# :math:`q \in (0, 1)`. The weights or coefficients :math:`w` are then found by
# the following minimization problem:
#
# .. math::
#     \min_{w} {\frac{1}{n_{\text{samples}}}
#     \sum_i PB_q(y_i - X_i w) + \alpha ||w||_1}.
#
# This consists of the pinball loss (also known as linear loss),
# see also :class:`~sklearn.metrics.mean_pinball_loss`,
#
# .. math::
#     PB_q(t) = q \max(t, 0) + (1 - q) \max(-t, 0) =
#     \begin{cases}
#         q t, & t > 0, \\
#         0,    & t = 0, \\
#         (1-q) t, & t < 0
#     \end{cases}
#
# and the L1 penalty controlled by parameter ``alpha``, similar to
# :class:`linear_model.Lasso`.
#
# Please check the following example to see how it works, and the :ref:`User
# Guide <quantile_regression>` for more details.
#
# .. figure:: ../linear_model/images/sphx_glr_plot_quantile_regression_002.png
#    :target: ../linear_model/plot_quantile_regression.html
#    :align: center
#    :scale: 50%

##############################################################################
# Feature Names Support
# --------------------------------------------------------------------------
# When an estimator is passed a `pandas' dataframe
# <https://pandas.pydata.org/docs/user_guide/dsintro.html#dataframe>`_ during
# :term:`fit`, the estimator will set a `feature_names_in_` attribute
# containing the feature names. Note that feature names support is only enabled
# when the column names in the dataframe are all strings. `feature_names_in_`
# is used to check that the column names of the dataframe passed in
# non-:term:`fit`, such as :term:`predict`, are consistent with features in
# :term:`fit`:
from sklearn.preprocessing import StandardScaler
import pandas as pd

X = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["a", "b", "c"])
scalar = StandardScaler().fit(X)
scalar.feature_names_in_

# %%
# The support of :term:`get_feature_names_out` is available for transformers
# that already had `get_feature_names` and transformers with a one-to-one
# correspondence between input and output such as
# :class:`~preprocessing.StandardScaler`. :term:`get_feature_names_out` support
# will be added to all other transformers in future releases. Additionally,
# :meth:`compose.ColumnTransformer.get_feature_names_out` is available to
# combine feature names of its transformers:
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
import pandas as pd

X = pd.DataFrame({"pet": ["dog", "cat", "fish"], "age": [3, 7, 1]})
preprocessor = ColumnTransformer(
    [
        ("numerical", StandardScaler(), ["age"]),
        ("categorical", OneHotEncoder(), ["pet"]),
    ],
    verbose_feature_names_out=False,
).fit(X)

preprocessor.get_feature_names_out()

# %%
# When this ``preprocessor`` is used with a pipeline, the feature names used
# by the classifier are obtained by slicing and calling
# :term:`get_feature_names_out`:
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline

y = [1, 0, 1]
pipe = make_pipeline(preprocessor, LogisticRegression())
pipe.fit(X, y)
pipe[:-1].get_feature_names_out()


##############################################################################
# A more flexible plotting API
# --------------------------------------------------------------------------
# :class:`metrics.ConfusionMatrixDisplay`,
# :class:`metrics.PrecisionRecallDisplay`, :class:`metrics.DetCurveDisplay`,
# and :class:`inspection.PartialDependenceDisplay` now expose two class
# methods: `from_estimator` and `from_predictions` which allow users to create
# a plot given the predictions or an estimator. This means the corresponding
# `plot_*` functions are deprecated. Please check :ref:`example one
# <sphx_glr_auto_examples_model_selection_plot_confusion_matrix.py>` and
# :ref:`example two
# <sphx_glr_auto_examples_classification_plot_digits_classification.py>` for
# how to use the new plotting functionalities.

##############################################################################
# Online One-Class SVM
# --------------------------------------------------------------------------
# The new class :class:`~linear_model.SGDOneClassSVM` implements an online
# linear version of the One-Class SVM using a stochastic gradient descent.
# Combined with kernel approximation techniques,
# :class:`~linear_model.SGDOneClassSVM` can be used to approximate the solution
# of a kernelized One-Class SVM, implemented in :class:`~svm.OneClassSVM`, with
# a fit time complexity linear in the number of samples. Note that the
# complexity of a kernelized One-Class SVM is at best quadratic in the number
# of samples. :class:`~linear_model.SGDOneClassSVM` is thus well suited for
# datasets with a large number of training samples (> 10,000) for which the SGD
# variant can be several orders of magnitude faster. Please check this
# :ref:`example
# <sphx_glr_auto_examples_miscellaneous_plot_anomaly_comparison.py>` to see how
# it's used, and the :ref:`User Guide <sgd_online_one_class_svm>` for more
# details.
#
# .. figure:: ../miscellaneous/images/sphx_glr_plot_anomaly_comparison_001.png
#    :target: ../miscellaneous/plot_anomaly_comparison.html
#    :align: center

##############################################################################
# Histogram-based Gradient Boosting Models are now stable
# --------------------------------------------------------------------------
# :class:`~sklearn.ensemble.HistGradientBoostingRegressor` and
# :class:`~ensemble.HistGradientBoostingClassifier` are no longer experimental
# and can simply be imported and used as::
#
#     from sklearn.ensemble import HistGradientBoostingClassifier

##############################################################################
# New documentation improvements
# ------------------------------
# This release includes many documentation improvements. Out of over 2100
# merged pull requests, about 800 of them are improvements to our
# documentation.
