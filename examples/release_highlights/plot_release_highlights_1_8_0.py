# ruff: noqa: CPY001
"""
=======================================
Release Highlights for scikit-learn 1.8
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.8! Many bug fixes
and improvements were added, as well as some key new features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_8>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Array API support
# -----------------
# The progressive adoption of the Python array API standard in
# scikit-learn means that PyTorch and CuPy input arrays
# are used directly. This means that in scikit-learn estimators
# and functions non-CPU devices, such as GPUs, can be used
# to perform the computation. As a result performance is improved
# and integration with these libraries is easier.
#
# In scikit-learn 1.8, several estimators and functions have been updated to
# support array API compatible inputs, for example PyTorch tensors and CuPy
# arrays.
#
# Array API support was added to the following estimators:
# :class:`preprocessing.StandardScaler`,
# :class:`preprocessing.PolynomialFeatures`, :class:`linear_model.RidgeCV`,
# :class:`linear_model.RidgeClassifierCV`, :class:`mixture.GaussianMixture` and
# :class:`calibration.CalibratedClassifierCV`.
#
# Array API support was also added to several metrics in :mod:`sklearn.metrics`
# module, see :ref:`array_api_supported` for more details.
#
# Please refer to the :ref:`array API support<array_api>` page for instructions
# to use scikit-learn with array API compatible libraries such as PyTorch or CuPy.
# Note: Array API support is experimental and must be explicitly enabled both
# in SciPy and scikit-learn.
#
# Here an excerpt of using :class:`calibration.CalibratedClassifierCV` and
# :class:`linear_model.RidgeCV` together on a GPU with the help of PyTorch:
#
# .. code-block:: python
#
#     ridge_pipeline_gpu = make_pipeline(
#         TableVectorizer(
#             numeric=make_pipeline(
#                 QuantileTransformer(),
#                 SplineTransformer(n_knots=10),
#             ),
#             high_cardinality=TargetEncoder(cv=5),
#         ),
#         FunctionTransformer(
#             lambda x: torch.tensor(x.to_numpy().astype(np.float32), device="cuda"))
#         ,
#         CalibratedClassifierCV(
#             RidgeClassifierCV(alphas=alphas), method="temperature"
#         ),
#     )
#     with sklearn.config_context(array_api_dispatch=True):
#         cv_results = cross_validate(ridge_pipeline_gpu, features, target)
#
#
# A [full notebook of this example on Google
# Colab](https://colab.research.google.com/drive/1ztH8gUPv31hSjEeR_8pw20qShTwViGRx?usp=sharing).
# For this example, using the colab GPU vs using a single CPU core leads to a
# 10x speedup which is quite typical for such workloads.

# %%
# Free-threaded CPython 3.14 support
# ----------------------------------
#
# scikit-learn has support for free-threaded CPython, in particular
# free-threaded wheels are available for all of our supported platforms on Python
# 3.14.
#
# Free-threaded (also known as nogil) CPython is a version of CPython that aims to
# enable efficient multi-threaded use cases by removing the Global Interpreter Lock
# (GIL).
#
# Please try your use cases with free-threaded CPython and `report issues
# <https://github.com/scikit-learn/scikit-learn/issues/new/choose>`_!
#
# For more details about free-threaded CPython see `py-free-threading doc
# <https://py-free-threading.github.io>`_, in particular `how to install a
# free-threaded CPython <https://py-free-threading.github.io/installing_cpython/>`_
# and `Ecosystem compatibility tracking <https://py-free-threading.github.io/tracking/>`_.
#
# The long term goal of free-threaded Python is to more efficiently leverage
# multi-core CPUs by using thread workers instead of subprocess workers for parallel
# computation when passing `n_jobs>1` in functions or estimators.
# Efficiency gains are expected by removing the need for inter-process communication.
# Note however that process-based parallelism is still the default joblib backend at
# the time of writing. You can call `joblib.parallel_config(backend="threading")` to
# change the default backend to "threading". Be aware that properly testing that
# everything is running smoothly when doing so is still an ongoing effort and that
# there are open issues to fix before considering making this the default.

# %%
# Temperature scaling in `CalibratedClassifierCV`
# -----------------------------------------------
# Probability calibration of classifiers with temperature scaling is available in
# :class:`calibration.CalibratedClassifierCV` by setting `method = "temperature"`.
# This method is particularly well suited for multiclass problems because it provides
# (better-) calibrated probabilities with just one free parameter, in contrast
# to using a "One-vs-Rest" scheme that adds more parameters for each class.

from sklearn.calibration import CalibratedClassifierCV
from sklearn.datasets import make_classification
from sklearn.naive_bayes import GaussianNB

X, y = make_classification(n_classes=3, n_informative=8, random_state=42)
clf = GaussianNB()
sig = CalibratedClassifierCV(clf, method="sigmoid", ensemble=False).fit(X, y)  # old
ts = CalibratedClassifierCV(clf, method="temperature", ensemble=False).fit(X, y)  # new

# %%
# The following example shows that temperature scaling can produce better-calibrated
# probabilities than sigmoid calibration in multi-class classification problem
# (3 classes).

import matplotlib.pyplot as plt

from sklearn.calibration import CalibrationDisplay

fig, axes = plt.subplots(
    figsize=(8, 4.5),
    ncols=3,
    sharey=True,
)
for i, c in enumerate(ts.classes_):
    CalibrationDisplay.from_predictions(
        y == c,
        ts.predict_proba(X)[:, i],
        name="Temperature scaling",
        ax=axes[i],
    )
    CalibrationDisplay.from_predictions(
        y == c,
        sig.predict_proba(X)[:, i],
        name="Sigmoid",
        ax=axes[i],
    )
    axes[i].set_title(f"Class {c}")
    axes[i].set_xlabel(None)
    axes[i].set_ylabel(None)
    axes[i].get_legend().remove()
fig.suptitle("Reliability Diagrams per Class")
fig.supxlabel("Mean Predicted Probability")
fig.supylabel("Fraction of Class")
fig.legend(*axes[0].get_legend_handles_labels(), loc=(0.72, 0.5))
plt.subplots_adjust(right=0.7)
_ = fig.show()

# %%
# Linear models improvements
# --------------------------
# There are two main developments going on for linear models: efficiency and API
# changes.
#
# Efficiency of squared error based models with L1 penalty
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# The first one is a massive improvement of efficiency, in particular a reduced fit
# time, for squared error based estimators with L1 penalty: `ElasticNet`, `Lasso`,
# `MultiTaskElasticNet`, `MultiTaskLasso` and their CV variants. The fit time
# improvement is mainly achieved by **gap safe screening rules**. They enable the
# coordinate descent solver to set feature coefficients to 0 early and not look at them
# again. The stronger the L1 penalty the earlier features can be excluded from further
# updates.

from time import time

from sklearn.datasets import make_regression
from sklearn.linear_model import ElasticNetCV

X, y = make_regression(n_features=5000)
model = ElasticNetCV()
tic = time()
model.fit(X, y)
toc = time()
print(f"Fitting ElasticNetCV took {toc - tic:.3} seconds.")

# %%
# API changes in logistic regression
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# After the deprecation (version 1.5) and removal (now in 1.8) of the `multi_class`
# parameter in `LogisticRegression` and `LogisticRegressionCV`, the bookkeeping goes
# on. The goal for `LogisticRegressionCV` is to fix types and shapes of a few fitted
# attributes: `C_`, `l1_ratio_`, `coefs_paths_`, `scores_`, `n_iter_`. To ease the
# transition, you can choose between new and old types and shapes with the new and
# temporary parameter `use_legacy_attributes`.

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegressionCV

X, y = make_classification(n_classes=3, n_informative=4)
model = LogisticRegressionCV(
    use_legacy_attributes=True, l1_ratios=(0,), solver="newton-cholesky"
).fit(X, y)
model.C_  # ndarray of shape (3,), 3 times the same value
model = LogisticRegressionCV(
    use_legacy_attributes=False, l1_ratios=(0,), solver="newton-cholesky"
).fit(X, y)
model.C_  # single float

# %%
# A further deprecation is related to the `penalty` parameter in both
# `LogisticRegression` and `LogisticRegressionCV`. It is redundant because `C` together
# with `l1_ratio` for `LogisticRegression` and `l1_ratios` for `LogisticRegressionCV`
# contains the same information. Removing `penalty` can ease specifying grid search
# parameter spaces and provides a tidier API.

# %%
# HTML representation of estimators
# ---------------------------------
# Hyperparameters in the dropdown table of the HTML representation now include
# links to the online documentation. Docstring descriptions are also shown as
# tooltips on hover.

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

clf = make_pipeline(StandardScaler(), LogisticRegression(random_state=0, C=10))

# %%
# Expand the estimator diagram below by clicking on "LogisticRegression" and then on
# "Parameters" below.

clf


# %%
# DecisionTreeRegressor with `criterion="absolute_error"`
# ------------------------------------------------------
# :class:`tree.DecisionTreeRegressor` with `criterion="absolute_error"`
# now runs much faster. It has now `O(n * log(n))` complexity compared to
# `O(n**2)` previously, which allows to scale to millions of data points.
#
# As an illustration, on a dataset with 100_000 samples and 1 feature, doing a
# single split takes of the order of 100 ms, compared to ~20 seconds before.

import time

from sklearn.datasets import make_regression
from sklearn.tree import DecisionTreeRegressor

X, y = make_regression(n_samples=100_000, n_features=1)
tree = DecisionTreeRegressor(criterion="absolute_error", max_depth=1)

tic = time.time()
tree.fit(X, y)
elapsed = time.time() - tic
print(f"Fit took {elapsed:.2f} seconds")

# %%
# ClassicalMDS
# ------------
# Classical MDS, also known as "Principal Coordinates Analysis (PCoA)"
# or "Torgerson's scaling" is now available within the `sklearn.manifold`
# module. Classical MDS is close to PCA and instead of of approximating
# distances, it approximates pairwise scalar products, which has an exact
# analytic solution in terms of eigendecomposition.
#
# Let's illustrate this new addition by using it on a S-curve dataset to
# get a low-dimensional representation of the data.

import matplotlib.pyplot as plt
from matplotlib import ticker

from sklearn import datasets, manifold

n_samples = 1500
S_points, S_color = datasets.make_s_curve(n_samples, random_state=0)
md_classical = manifold.ClassicalMDS(n_components=2)
S_scaling = md_classical.fit_transform(S_points)

fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
x, y, z = S_points.T
ax1.scatter(x, y, z, c=S_color, s=50, alpha=0.8)
ax1.set_title("Original S-curve samples", size=16)
ax1.view_init(azim=-60, elev=9)
for axis in (ax1.xaxis, ax1.yaxis, ax1.zaxis):
    axis.set_major_locator(ticker.MultipleLocator(1))

ax2 = fig.add_subplot(1, 2, 2)
x2, y2 = S_scaling.T
ax2.scatter(x2, y2, c=S_color, s=50, alpha=0.8)
ax2.set_title("Classical MDS", size=16)
for axis in (ax2.xaxis, ax2.yaxis):
    axis.set_major_formatter(ticker.NullFormatter())

plt.show()
