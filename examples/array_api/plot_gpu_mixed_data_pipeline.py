"""
==============================================
Mixed-Type Data Pipeline with GPU Acceleration
==============================================

Demonstrates how to build a pipeline that processes mixed-type data -- text
and numeric features -- where text processing is inherently CPU-bound but
downstream computation benefits from GPU acceleration.

``TfidfVectorizer`` operates on strings and can only run on CPU. Numeric
features also arrive as NumPy arrays or pandas columns. After
``ColumnTransformer`` combines both feature types into a single dense matrix,
a ``FunctionTransformer`` moves it to GPU with one call to ``torch.asarray``.
Everything downstream -- PCA and logistic regression -- then runs on GPU via
Array API dispatch.

This pattern generalises to any pipeline where some preprocessing steps are
CPU-only (text vectorisation, categorical encoding, etc.) but the
compute-heavy model fitting benefits from GPU acceleration.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys
import time
from functools import partial

os.environ["SCIPY_ARRAY_API"] = "1"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import torch
except ImportError:
    print("This example requires PyTorch. Install it with: pip install torch")
    sys.exit(1)

import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_20newsgroups
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load dataset and build a DataFrame
# ------------------------------------
# We construct a pandas DataFrame with a text column alongside numeric
# metadata features extracted from the raw messages: character count, line
# count, and the fraction of quoted lines (lines starting with ``>``).

newsgroups = fetch_20newsgroups(subset="all")

texts = newsgroups.data
df = pd.DataFrame({
    "text": texts,
    "n_chars": [len(t) for t in texts],
    "n_lines": [t.count("\n") + 1 for t in texts],
    "quote_ratio": [
        sum(1 for line in t.split("\n") if line.startswith(">"))
        / max(t.count("\n") + 1, 1)
        for t in texts
    ],
})
y = newsgroups.target

df_train, df_test, y_train, y_test = train_test_split(
    df, y, test_size=0.2, random_state=42, stratify=y
)

print(f"Train: {len(df_train)}, Test: {len(df_test)}")
print(f"Columns: {list(df.columns)}")

# %%
# Define the preprocessing
# ------------------------
# ``ColumnTransformer`` applies TF-IDF to the text column and standard
# scaling to the numeric columns. ``sparse_threshold=0`` ensures the
# output is always a dense NumPy array, which ``torch.asarray`` can consume.

numeric_features = ["n_chars", "n_lines", "quote_ratio"]


def make_preprocessor():
    return ColumnTransformer(
        [
            ("tfidf", TfidfVectorizer(max_features=5000), "text"),
            ("numeric", StandardScaler(), numeric_features),
        ],
        sparse_threshold=0,
    )


# %%
# CPU pipeline
# ------------
# All steps run on CPU with NumPy arrays.

pipe_cpu = make_pipeline(
    make_preprocessor(),
    PCA(svd_solver="full", n_components=128),
    LogisticRegression(solver="lbfgs", max_iter=1000),
)

t0 = time.perf_counter()
pipe_cpu.fit(df_train, y_train)
y_pred_cpu = pipe_cpu.predict(df_test)
elapsed_cpu = time.perf_counter() - t0

acc_cpu = accuracy_score(y_test, y_pred_cpu)
print(f"CPU pipeline:  {elapsed_cpu:.1f}s, accuracy={acc_cpu:.4f}")

# %%
# GPU pipeline
# ------------
# Identical to the CPU pipeline except for one extra step: a
# ``FunctionTransformer`` that moves the combined feature matrix to GPU.
# ``ColumnTransformer`` (TF-IDF + scaling) still runs on CPU -- text
# vectorisation cannot run on GPU. But PCA and logistic regression now
# operate on CUDA tensors.

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

pipe_gpu = make_pipeline(
    make_preprocessor(),
    FunctionTransformer(partial(torch.asarray, dtype=torch.float32, device=device)),
    PCA(svd_solver="full", n_components=128),
    LogisticRegression(solver="lbfgs", max_iter=1000),
)

t0 = time.perf_counter()
pipe_gpu.fit(df_train, y_train)
y_pred_gpu = pipe_gpu.predict(df_test)
if device == "cuda":
    torch.cuda.synchronize()
elapsed_gpu = time.perf_counter() - t0

y_pred_gpu_np = (
    np.asarray(y_pred_gpu.cpu()) if hasattr(y_pred_gpu, "cpu") else y_pred_gpu
)
acc_gpu = accuracy_score(y_test, y_pred_gpu_np)
print(f"GPU pipeline:  {elapsed_gpu:.1f}s, accuracy={acc_gpu:.4f}")

sklearn.set_config(array_api_dispatch=False)

# %%
# Compare
# -------
# Both pipelines produce the same accuracy. The text vectorisation step
# (TF-IDF) is CPU-bound in both cases, but PCA and logistic regression
# benefit from GPU acceleration. On larger datasets or with more expensive
# downstream models, the GPU advantage grows.
#
# The key takeaway is the pattern: a single ``FunctionTransformer`` line
# bridges CPU-only preprocessing to GPU-accelerated modelling, with no
# changes to any other pipeline step.

fig, ax = plt.subplots(figsize=(7, 4))
labels = ["CPU pipeline\n(all NumPy)", f"GPU pipeline\n(FunctionTransformer → {device.upper()})"]
bars = ax.bar(labels, [elapsed_cpu, elapsed_gpu])
ax.set_ylabel("Wall-clock time (seconds)")
ax.set_title("Mixed-type pipeline: TF-IDF + numeric → PCA → LogisticRegression")
ax.bar_label(bars, fmt="%.1fs")
plt.tight_layout()
plt.show()
