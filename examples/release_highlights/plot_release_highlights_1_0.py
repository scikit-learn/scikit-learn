# flake8: noqa
"""
=======================================
Release Highlights for scikit-learn 1.0
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.0! The library has
been stable for quite some time, releasing version 1.0 is recognizing that and
signalling it to our users. The release does not include any breaking changes,
and as usual, and we do our best to follow a two release deprecation cycles for
any breaking changes.

This release includes many bug fixes and improvements, as well as some new key
features. We detail below a few of the major features of this release. **For an
exhaustive list of all the changes**, please refer to the :ref:`release notes
<changes_1_0>`.

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
# :class:`ensemble.HistGradientBoostingRegressor` as::
#
#         HistGradientBoostingRegressor("squared_error", 0.1, 100, 31, None,
#             20, 0.0, 255, None, None, False, "auto", "loss", 0.1, 10, 1e-7,
#             0, None)
#
# Understanding the above code requires the reader to go to the API
# documentation to check each parameter, and what they mean, with the position
# they have. To improve the readability of code written based on scikit-learn,
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
# which is much more readable.

##############################################################################
# Spline Transformers
# ---------------------------------------------------------
# One way to add nonlinear terms to the dataset's feature set is to generate
# spline basis functions for each feature with the new
# :class:`preprocessing.SplineTransformer`. Splines are piecewise polynomials,
# parametrized by their polynomial degree and the positions of the knots. The
# :class:`preprocessing.SplineTransformer` implements a B-spline basis

# .. figure:: ../linear_model/images/sphx_glr_plot_polynomial_interpolation_001.png
#   :target: ../linear_model/plot_successive_halving_iterations.html
#   :align: center

# The following code shows splines in action, for more information, please
# refer to :ref:`User Guide <spline_transformer>`

import numpy as np
from sklearn.preprocessing import SplineTransformer

X = np.arange(5).reshape(5, 1)
X
spline = SplineTransformer(degree=2, n_knots=3)
spline.fit_transform(X)


##############################################################################
# Quantile Regressor
# --------------------------------------------------------------------------
# Quantile regression estimates the median or other quantiles of :math:`y`
# conditional on :math:`X`, while ordinary least squares (OLS) estimates the
# conditional mean.
#
# As a linear model, the new :class:`QuantileRegressor` gives linear predictions
# :math:`\hat{y}(w, X) = Xw` for the :math:`q`-th quantile, :math:`q \in (0,
# 1)`. The weights or coefficients :math:`w` are then found by the following
# minimization problem:

# .. math::
#     \min_{w} {\frac{1}{n_{\text{samples}}}
#     \sum_i PB_q(y_i - X_i w) + \alpha ||w||_1}.

# This consists of the pinball loss (also known as linear loss),
# see also :class:`~sklearn.metrics.mean_pinball_loss`,

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

# .. figure:: /auto_examples/linear_model/images/sphx_glr_plot_quantile_regression_002.png
#    :target: ../auto_examples/linear_model/plot_quantile_regression.html
#    :align: center
#    :scale: 50%

##############################################################################
# `from_estimator` and `from_predictions` on `*Display` classes
# --------------------------------------------------------------------------
# :class:`metrics.ConfusionMatrixDisplay`,
# :class:`metrics.PrecisionRecallDisplay`, and :class:`metrics.DetCurveDisplay`
# now expose two class methods: `from_estimator` and `from_predictions` which
# allow users to create a plot given the predictions or an estimator. This
# means the corresponsing `plot_*` functions are deprecated. Please check
# :ref:`example one
# <sphx_glr_auto_examples_model_selection_plot_confusion_matrix.py>` and
# :ref:`example two
# <sphx_glr_auto_examples_classification_plot_digits_classification.py>` for
# examples of how to use the new plotting functionalities.

##############################################################################
# Online One-Class SVM
# --------------------------------------------------------------------------
# The class :class:`sklearn.linear_model.SGDOneClassSVM` implements an online
# linear version of the One-Class SVM using a stochastic gradient descent.
# Combined with kernel approximation techniques,
# :class:`sklearn.linear_model.SGDOneClassSVM` can be used to approximate the
# solution of a kernelized One-Class SVM, implemented in
# :class:`sklearn.svm.OneClassSVM`, with a linear complexity in the number of
# samples. Note that the complexity of a kernelized One-Class SVM is at best
# quadratic in the number of samples.
# :class:`sklearn.linear_model.SGDOneClassSVM` is thus well suited for datasets
# with a large number of training samples (> 10,000) for which the SGD variant
# can be several orders of magnitude faster. Please check this :ref:`example
# <sphx_glr_auto_examples_miscellaneous_plot_anomaly_comparison.py` to see how
# it's used, and :ref:`User Guide <sgd_online_one_class_svm>` for more details.

# .. figure:: /auto_examples/miscellaneous/images/sphx_glr_plot_anomaly_comparison_001.png
#    :target: ../auto_examples/miscellaneous/plot_anomaly_comparison.html
#    :align: center

##############################################################################
# Histogram Based Gradient Boosting Models
# --------------------------------------------------------------------------
# `ensemble.HistGradientBoostingRegressor` and
# `ensemble.HistGradientBoostingClassifier` are no longer experimental and can
# simply be imported and used as::
#
#     from sklearn.ensemble import HistGradientBoostingClassifier

##############################################################################
# New documentation improvements
# ------------------------------
# TODO
