"""
========================================================================
Accelerating a Mixed-Type Pipeline on GPU with the Python Array API
========================================================================

This example demonstrates how to use scikit-learn's experimental
:ref:`array_api` support to accelerate cross-validation of a non-trivial
regression pipeline on hardware accelerators (NVIDIA GPU via CUDA, Apple
Silicon via MPS, or Intel GPU via XPU) using PyTorch.

Real-world datasets often contain heterogeneous data types (both categorical
and numerical columns). At present, the Python Array API standard is designed
exclusively for numeric arrays and tensors; non-numeric data (such as string
categories) is outside its scope. In addition, some preprocessing transformers
(such as :class:`~sklearn.preprocessing.SplineTransformer`) currently only run
on CPU.

To work around these limitations in a realistic, production-friendly workflow:

1. A :class:`~sklearn.compose.ColumnTransformer` performs initial preprocessing
   on the **CPU** (one-hot encoding categorical features with
   :class:`~sklearn.preprocessing.OneHotEncoder` and generating non-linear
   splines with :class:`~sklearn.preprocessing.SplineTransformer`).
2. A custom :class:`~sklearn.preprocessing.FunctionTransformer` automatically
   detects the host's fastest hardware accelerator and transfers the dense
   numeric feature matrix to a PyTorch tensor on that device.
3. The downstream modeling steps—a non-linear kernel approximation via
   :class:`~sklearn.kernel_approximation.Nystroem` and regularized regression via
   :class:`~sklearn.linear_model.RidgeCV`—run with full hardware acceleration
   on the detected device via Array API dispatch.

Because the device detection happens dynamically at runtime, the pipeline
remains fully serializable (picklable) and can be saved and loaded across
machines with different hardware configurations without any code modifications.

.. note::
   Running this example on a GPU:

   Scikit-learn's documentation is rendered in a CPU-only environment on
   Continuous Integration (CI). While this example runs cleanly on CPU,
   speedups will only be visible on a machine equipped with a supported
   accelerator (CUDA, MPS, or XPU).

   You can try this example interactively with a free GPU on Google Colab:

   .. image:: https://colab.research.google.com/assets/colab-badge.svg
      :target: https://colab.research.google.com/github/scikit-learn/scikit-learn/blob/main/doc/auto_examples/applications/plot_array_api_gpu_pipeline.ipynb
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Setting Up Environment for Array API
# ------------------------------------
#
# Scikit-learn's Array API support currently requires setting the
# ``SCIPY_ARRAY_API=1`` environment variable before importing SciPy or
# scikit-learn.
#
# If running inside Google Colab, uncomment and run the following lines to
# install PyTorch and the latest development build of scikit-learn:
#
# .. code-block:: python
#
#     # !pip install torch
#     # !pip install --pre --extra-index \\
#     #     https://pypi.anaconda.org/scientific-python-nightly-wheels/simple \\
#     #     scikit-learn

import os

os.environ["SCIPY_ARRAY_API"] = "1"

import pickle
import time

import numpy as np
import pandas as pd
import torch

from sklearn import config_context
from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_openml
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import (
    FunctionTransformer,
    OneHotEncoder,
    SplineTransformer,
)

# %%
# Automatic Hardware Device Detection
# -----------------------------------
#
# We define helper functions to detect the fastest available compute device
# supported by PyTorch (CUDA > MPS > XPU > CPU) and move data accordingly.


def get_best_torch_device():
    """Detect the best available PyTorch compute device."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    elif hasattr(torch, "xpu") and torch.xpu.is_available():
        return torch.device("xpu")
    return torch.device("cpu")


def to_best_device(X):
    """Convert input array to a PyTorch tensor on the best detected device."""
    device = get_best_torch_device()
    if isinstance(X, torch.Tensor):
        return X.to(device)
    if hasattr(X, "toarray"):
        X = X.toarray()
    return torch.asarray(X, dtype=torch.float32, device=device)


# %%
# Loading a Mixed-Type Regression Dataset
# ---------------------------------------
#
# We load the Ames Housing dataset, which features a realistic combination of
# categorical variables (neighborhood, building type, zoning) and numerical
# measurements (living area, year built, lot area).

try:
    ames = fetch_openml("house_prices", as_frame=True)
    df = ames.frame

    categorical_cols = ["MSZoning", "Neighborhood", "BldgType", "HouseStyle"]
    numerical_cols = ["GrLivArea", "TotalBsmtSF", "YearBuilt", "LotArea"]

    X = df[categorical_cols + numerical_cols].copy()
    y = (df["SalePrice"] / 1000.0).astype(np.float32)

    for col in numerical_cols:
        X[col] = X[col].fillna(X[col].median())
    for col in categorical_cols:
        X[col] = X[col].fillna("Missing")

except (OSError, ValueError, RuntimeError):
    # Offline fallback in case network access to OpenML is unavailable
    rng = np.random.default_rng(42)
    n_samples = 1460
    X = pd.DataFrame(
        {
            "MSZoning": rng.choice(["RL", "RM", "C", "FV"], size=n_samples),
            "Neighborhood": rng.choice(
                ["CollgCr", "Veenker", "Crawfor", "NoRidge"], size=n_samples
            ),
            "BldgType": rng.choice(
                ["1Fam", "2fmCon", "Duplex", "Twnhs"], size=n_samples
            ),
            "HouseStyle": rng.choice(["1Story", "2Story", "1.5Fin"], size=n_samples),
            "GrLivArea": rng.uniform(500, 3500, size=n_samples),
            "TotalBsmtSF": rng.uniform(0, 2500, size=n_samples),
            "YearBuilt": rng.integers(1950, 2020, size=n_samples),
            "LotArea": rng.uniform(2000, 25000, size=n_samples),
        }
    )
    y = (
        X["GrLivArea"] * 0.08
        + X["YearBuilt"] * 0.05
        + rng.normal(0, 15, size=n_samples)
    ).astype(np.float32)
    categorical_cols = ["MSZoning", "Neighborhood", "BldgType", "HouseStyle"]
    numerical_cols = ["GrLivArea", "TotalBsmtSF", "YearBuilt", "LotArea"]

print(f"Dataset shape: {X.shape}")
print(f"Categorical features: {categorical_cols}")
print(f"Numerical features: {numerical_cols}")

# %%
# Building the Hybrid Pipeline
# ----------------------------
#
# Our pipeline consists of three phases:
#
# 1. **CPU Preprocessing:** :class:`~sklearn.compose.ColumnTransformer` processes
#    categorical and numerical columns using standard NumPy/pandas objects.
# 2. **Device Transfer:** :class:`~sklearn.preprocessing.FunctionTransformer`
#    converts the transformed array into a PyTorch tensor on the chosen device.
# 3. **Hardware-Accelerated Modeling:**
#    :class:`~sklearn.kernel_approximation.Nystroem` projects samples to an
#    RBF reproducing kernel space, and :class:`~sklearn.linear_model.RidgeCV`
#    optimizes regularization using generalized cross-validation (GCV).

cpu_preprocessor = ColumnTransformer(
    transformers=[
        (
            "cat",
            OneHotEncoder(sparse_output=False, drop="first", handle_unknown="ignore"),
            categorical_cols,
        ),
        (
            "num",
            SplineTransformer(n_knots=5, degree=3),
            numerical_cols,
        ),
    ]
)

device_transfer = FunctionTransformer(to_best_device)

gpu_pipeline = Pipeline(
    [
        ("preprocessor", cpu_preprocessor),
        ("transfer", device_transfer),
        ("nystroem", Nystroem(kernel="rbf", n_components=300, random_state=42)),
        ("regressor", RidgeCV(alphas=np.logspace(-3, 3, 10))),
    ]
)

# %%
# Cross-Validation with Array API Dispatching
# -------------------------------------------
#
# We enable Array API support with ``config_context(array_api_dispatch=True)``.
# Array-API-compatible estimators (:class:`~sklearn.kernel_approximation.Nystroem`
# and :class:`~sklearn.linear_model.RidgeCV`) dispatch linear algebra routines
# directly to the hardware device backing the input tensor.

current_device = get_best_torch_device()
print(f"Executing pipeline on device: {current_device}")

cv = KFold(n_splits=5, shuffle=True, random_state=42)

t0 = time.perf_counter()
with config_context(array_api_dispatch=True):
    cv_results = cross_validate(
        gpu_pipeline,
        X,
        y,
        cv=cv,
        scoring="neg_root_mean_squared_error",
        return_train_score=True,
    )
duration = time.perf_counter() - t0

test_rmse = -cv_results["test_score"]
train_rmse = -cv_results["train_score"]

print(f"5-Fold CV Duration: {duration:.2f} seconds")
print(f"Mean Test RMSE: ${test_rmse.mean():.2f}k +/- ${test_rmse.std():.2f}k")
print(f"Mean Train RMSE: ${train_rmse.mean():.2f}k +/- ${train_rmse.std():.2f}k")

# %%
# Serialization and Cross-Device Pickling
# ---------------------------------------
#
# Because device detection is encapsulated inside the transfer step, the fitted
# pipeline can be pickled and restored on other machines regardless of their
# hardware accelerators.

with config_context(array_api_dispatch=True):
    gpu_pipeline.fit(X, y)
    initial_preds = gpu_pipeline.predict(X[:3])
    print(f"Initial predictions (device={initial_preds.device}):\n{initial_preds}")

# Serialize and deserialize pipeline
pipeline_bytes = pickle.dumps(gpu_pipeline)
restored_pipeline = pickle.loads(pipeline_bytes)

with config_context(array_api_dispatch=True):
    restored_preds = restored_pipeline.predict(X[:3])
    print(f"Restored predictions (device={restored_preds.device}):\n{restored_preds}")
