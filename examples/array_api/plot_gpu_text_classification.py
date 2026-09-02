"""
=====================================================
GPU Text Classification with Sentence Embeddings
=====================================================

Demonstrates how Array API support simplifies the integration of GPU-based
embedding models with scikit-learn classifiers.

When using ``sentence-transformers`` with a GPU, the encoder produces
embeddings as PyTorch CUDA tensors. Without Array API support, these must be
copied back to NumPy before scikit-learn can use them -- an unnecessary
round-trip that wastes time and memory. With Array API dispatch enabled, the
embeddings flow directly from the encoder into the scikit-learn pipeline with
zero copies.

This example compares three workflows:

1. **NumPy (baseline)**: encode to NumPy, fit on CPU
2. **Manual transfer**: encode on GPU, copy to NumPy, fit on CPU
3. **Array API dispatch**: encode on GPU, fit on GPU -- no copies

This example requires ``sentence-transformers`` (``pip install
sentence-transformers``).
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys
import time

os.environ["SCIPY_ARRAY_API"] = "1"

import matplotlib.pyplot as plt
import numpy as np

try:
    import torch
except ImportError:
    print("This example requires PyTorch. Install it with: pip install torch")
    sys.exit(1)

try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    print(
        "This example requires sentence-transformers."
        " Install it with: pip install sentence-transformers"
    )
    sys.exit(1)

import sklearn
from sklearn.datasets import fetch_20newsgroups
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load dataset
# ------------
# 20 Newsgroups has ~11K train and ~7.5K test documents across 20 categories.

train = fetch_20newsgroups(subset="train")
test = fetch_20newsgroups(subset="test")

texts_train, y_train_np = train.data, train.target
texts_test, y_test_np = test.data, test.target
target_names = train.target_names

print(f"Train documents: {len(texts_train)}")
print(f"Test documents:  {len(texts_test)}")

# %%
# Set up the encoder and pipeline
# --------------------------------
# ``all-MiniLM-L6-v2`` produces 384-dim embeddings. The encoder always runs
# on GPU when ``device="cuda"`` -- the difference between workflows is what
# happens to the embeddings *after* encoding.

encoder = SentenceTransformer("all-MiniLM-L6-v2", device=device)

param_grid = {
    "pca__n_components": [64, 128, 256],
    "logisticregression__C": np.logspace(-2, 2, 5).tolist(),
}


def make_clf_pipeline():
    return make_pipeline(
        StandardScaler(),
        PCA(svd_solver="full", n_components=128),
        LogisticRegression(solver="lbfgs", max_iter=1000),
    )


# %%
# Workflow 1: NumPy baseline
# --------------------------
# Encode directly to NumPy arrays and fit on CPU. This is how most users
# integrate sentence-transformers with scikit-learn today.

t0 = time.perf_counter()

emb_train_np = encoder.encode(texts_train, convert_to_numpy=True, batch_size=128,
                              show_progress_bar=False)
emb_test_np = encoder.encode(texts_test, convert_to_numpy=True, batch_size=128,
                             show_progress_bar=False)

search_np = GridSearchCV(make_clf_pipeline(), param_grid, cv=3)
search_np.fit(emb_train_np, y_train_np)
y_pred_np = search_np.predict(emb_test_np)

elapsed_numpy = time.perf_counter() - t0
acc_numpy = accuracy_score(y_test_np, y_pred_np)
print(f"NumPy baseline:   {elapsed_numpy:.1f}s, accuracy={acc_numpy:.4f}")


# %%
# Workflow 2: Array API dispatch (zero-copy)
# -------------------------------------------
# With ``array_api_dispatch=True``, the CUDA tensors from the encoder feed
# directly into the scikit-learn pipeline. No ``.cpu().numpy()`` call, no
# extra memory allocation, and the pipeline's linear algebra runs on GPU.

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

t0 = time.perf_counter()

emb_train_gpu = encoder.encode(texts_train, convert_to_tensor=True, batch_size=128,
                               show_progress_bar=False)
emb_test_gpu = encoder.encode(texts_test, convert_to_tensor=True, batch_size=128,
                              show_progress_bar=False)

search_gpu = GridSearchCV(make_clf_pipeline(), param_grid, cv=3)
search_gpu.fit(emb_train_gpu, y_train_np)
y_pred_gpu = search_gpu.predict(emb_test_gpu)
if device == "cuda":
    torch.cuda.synchronize()

elapsed_gpu = time.perf_counter() - t0
acc_gpu = float(accuracy_score(y_test_np, y_pred_gpu))
print(f"Array API (GPU):  {elapsed_gpu:.1f}s, accuracy={acc_gpu:.4f}")

sklearn.set_config(array_api_dispatch=False)

# %%
# Compare workflows
# -----------------
# All three workflows produce the same accuracy. The Array API path is the
# simplest *and* the fastest: no manual transfers, no extra memory, and the
# downstream pipeline benefits from GPU acceleration.

fig, axes = plt.subplots(1, 2, figsize=(13, 4))

labels = ["NumPy\n(encode→numpy)", "Array API\n(encode→GPU→fit on GPU)"]
times = [elapsed_numpy, elapsed_gpu]
accs = [acc_numpy, acc_gpu]

bars = axes[0].bar(labels, times)
axes[0].set_ylabel("Wall-clock time (seconds)")
axes[0].set_title("End-to-end: embed + GridSearchCV on 20 Newsgroups")
axes[0].bar_label(bars, fmt="%.1fs")

axes[1].bar(labels, accs)
axes[1].set_ylabel("Accuracy")
axes[1].set_title("Test Accuracy (all workflows equivalent)")
axes[1].set_ylim(0.5, 0.85)

plt.tight_layout()
plt.show()

# %%
# Classification report
# ---------------------

y_pred_report = np.asarray(y_pred_gpu.cpu()) if hasattr(y_pred_gpu, "cpu") else y_pred_gpu
print(classification_report(y_test_np, y_pred_report, target_names=target_names))
