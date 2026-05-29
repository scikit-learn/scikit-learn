"""
=================================================
Array API support in CalibratedClassifierCV
=================================================

This example illustrates how :class:`~sklearn.calibration.CalibratedClassifierCV`
can be used with Array API-compatible inputs (such as PyTorch tensors,
CuPy arrays, or other libraries that implement the
`Array API standard <https://data-apis.org/array-api/latest/>`_).

When the input is an Array API namespace, the underlying estimator is fitted
in that namespace whenever it supports the Array API. The calibration step
itself uses NumPy internally (because the optimization routines and the
isotonic regression are NumPy/Cython-only), but inputs are converted to
NumPy and the resulting probabilities are converted back to the input
namespace so that the public API of :class:`CalibratedClassifierCV` stays
consistent with the input namespace.

In this example we use the ``array_api_strict`` namespace, a strict
reference implementation of the Array API that is useful for testing
namespace-agnostic code without requiring a GPU or another heavy
dependency. The same code works unchanged with NumPy, PyTorch, CuPy, or
any other compliant namespace.

.. note::

   To enable Array API dispatch, scikit-learn requires SciPy's own Array API
   support to be enabled. Set the ``SCIPY_ARRAY_API=1`` environment variable
   **before** importing scikit-learn or SciPy. For example, run this script
   with::

       SCIPY_ARRAY_API=1 python plot_calibration_array_api.py

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Generate a synthetic multiclass classification dataset
# ------------------------------------------------------
#
# We use a small synthetic dataset so that the example runs quickly.
import numpy as np

from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split

X, y = make_classification(
    n_samples=1000,
    n_features=10,
    n_informative=8,
    n_redundant=0,
    n_classes=3,
    n_clusters_per_class=1,
    class_sep=2.0,
    random_state=42,
)
X_train, X_calib, y_train, y_calib = train_test_split(X, y, random_state=42)

# %%
# Reference: NumPy path
# ---------------------
#
# We first fit a calibrated classifier in the usual NumPy way. This is the
# behavior every existing scikit-learn user is familiar with.
from sklearn.calibration import CalibratedClassifierCV
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.frozen import FrozenEstimator

base_estimator = LinearDiscriminantAnalysis().fit(X_train, y_train)

cal_numpy_sigmoid = CalibratedClassifierCV(
    FrozenEstimator(base_estimator), method="sigmoid", cv=3
).fit(X_calib, y_calib)
proba_numpy_sigmoid = cal_numpy_sigmoid.predict_proba(X_train)

cal_numpy_isotonic = CalibratedClassifierCV(
    FrozenEstimator(base_estimator), method="isotonic", cv=3
).fit(X_calib, y_calib)
proba_numpy_isotonic = cal_numpy_isotonic.predict_proba(X_train)

print("NumPy sigmoid probabilities (first row):", proba_numpy_sigmoid[0])
print("NumPy isotonic probabilities (first row):", proba_numpy_isotonic[0])

# %%
# Array API path
# --------------
#
# We now repeat the same workflow using arrays from a non-NumPy namespace.
# We pick ``array_api_strict`` for portability, but you can replace it with
# ``torch``, ``cupy``, or any other Array API-compatible namespace and the
# code below will work without modification (provided the namespace is
# installed and the device is available).
#
# Enabling Array API dispatch requires setting
# ``config_context(array_api_dispatch=True)``. Inside that context,
# scikit-learn estimators that support the Array API will preserve the
# input namespace and device in their outputs.
import array_api_strict as xp

from sklearn import config_context

X_train_xp = xp.asarray(X_train)
X_calib_xp = xp.asarray(X_calib)
y_train_xp = xp.asarray(y_train)
y_calib_xp = xp.asarray(y_calib)

with config_context(array_api_dispatch=True):
    base_estimator_xp = LinearDiscriminantAnalysis().fit(X_train_xp, y_train_xp)

    cal_xp_sigmoid = CalibratedClassifierCV(
        FrozenEstimator(base_estimator_xp), method="sigmoid", cv=3
    ).fit(X_calib_xp, y_calib_xp)
    proba_xp_sigmoid = cal_xp_sigmoid.predict_proba(X_train_xp)

    cal_xp_isotonic = CalibratedClassifierCV(
        FrozenEstimator(base_estimator_xp), method="isotonic", cv=3
    ).fit(X_calib_xp, y_calib_xp)
    proba_xp_isotonic = cal_xp_isotonic.predict_proba(X_train_xp)

# The output preserves the namespace of the input.
print("Array API sigmoid output type:", type(proba_xp_sigmoid).__module__)
print("Array API isotonic output type:", type(proba_xp_isotonic).__module__)

# %%
# Numerical equivalence between paths
# -----------------------------------
#
# Although the two paths use different array libraries, the calibrated
# probabilities are numerically equivalent (up to floating-point precision).
proba_xp_sigmoid_np = np.asarray(proba_xp_sigmoid)
proba_xp_isotonic_np = np.asarray(proba_xp_isotonic)

print(
    "Max sigmoid abs difference:",
    np.max(np.abs(proba_xp_sigmoid_np - proba_numpy_sigmoid)),
)
print(
    "Max isotonic abs difference:",
    np.max(np.abs(proba_xp_isotonic_np - proba_numpy_isotonic)),
)

# %%
# Summary
# -------
#
# :class:`~sklearn.calibration.CalibratedClassifierCV` returns numerically
# equivalent calibrated probabilities whether it receives NumPy arrays or
# arrays from another namespace that follows the Python array API standard
# [1]_, while preserving the input namespace and device in its output. The
# ``sigmoid`` and ``isotonic`` methods rely on NumPy and Cython routines
# internally, so their inputs are converted to NumPy and the results are
# converted back to the input namespace. Because the calibration step is
# usually cheap compared to fitting the base estimator, this conversion
# overhead is negligible.
#
# See the :ref:`Array API section of the calibration User Guide
# <calibration_array_api>` for further details and limitations.
#
# References
# ----------
#
# .. [1] `Python array API standard
#        <https://data-apis.org/array-api/latest/>`_
