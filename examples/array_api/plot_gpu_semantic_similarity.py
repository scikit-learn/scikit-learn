"""
======================================================
GPU Semantic Similarity and Near-Duplicate Detection
======================================================

Demonstrates GPU-accelerated pairwise cosine similarity on sentence
embeddings. Computing the full 18.8K x 18.8K similarity matrix is a dense
matrix multiplication -- the ideal GPU workload. Two practical applications
are shown: semantic search (top-K retrieval for a query document) and
near-duplicate detection (pairs above a cosine similarity threshold).

This example requires ``sentence-transformers`` (``pip install
sentence-transformers``).
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys
import time
import textwrap

os.environ["SCIPY_ARRAY_API"] = "1"

import matplotlib.pyplot as plt
import numpy as np

try:
    import torch
except ImportError:
    print("This example requires PyTorch. Install it with: pip install torch")
    sys.exit(0)

try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    print(
        "This example requires sentence-transformers."
        " Install it with: pip install sentence-transformers"
    )
    sys.exit(0)

import sklearn
from sklearn.datasets import fetch_20newsgroups
from sklearn.metrics.pairwise import cosine_similarity

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load dataset and compute embeddings
# ------------------------------------

newsgroups = fetch_20newsgroups(subset="all")
all_texts = newsgroups.data
all_targets = newsgroups.target
target_names = newsgroups.target_names

print(f"Documents: {len(all_texts)}")

encoder = SentenceTransformer("all-MiniLM-L6-v2", device=device)

embeddings_gpu = encoder.encode(
    all_texts, convert_to_tensor=True, batch_size=128, show_progress_bar=True
)
embeddings_np = np.asarray(
    embeddings_gpu.cpu() if hasattr(embeddings_gpu, "cpu") else embeddings_gpu
)
print(f"Embedding shape: {tuple(embeddings_gpu.shape)}")

# %%
# CPU / NumPy baseline: full pairwise similarity
# -----------------------------------------------
# Computing an (n x n) cosine similarity matrix from (n x 384) embeddings
# is dominated by a matrix multiplication.

t0 = time.perf_counter()
sim_cpu = cosine_similarity(embeddings_np)
elapsed_cpu = time.perf_counter() - t0
print(f"CPU cosine_similarity time: {elapsed_cpu:.2f}s")

# %%
# GPU / PyTorch: full pairwise similarity
# ----------------------------------------

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

t0 = time.perf_counter()
sim_gpu = cosine_similarity(embeddings_gpu)
if device == "cuda":
    torch.cuda.synchronize()
elapsed_gpu = time.perf_counter() - t0
print(f"GPU cosine_similarity time: {elapsed_gpu:.2f}s")

sklearn.set_config(array_api_dispatch=False)

# %%
# Compare timings
# ---------------
# With ~19K documents both CPU and GPU finish quickly, but the GPU is
# already an order of magnitude faster. The core operation is a dense
# matrix multiply (n x d) @ (d x n) that scales as O(n² · d). In
# production corpora with 100K--1M+ documents the CPU time grows
# quadratically into minutes or hours, while the GPU keeps up easily.

speedup = elapsed_cpu / elapsed_gpu if elapsed_gpu > 0 else float("inf")

fig, ax = plt.subplots(figsize=(6, 4))
bars = ax.bar(["NumPy / CPU", f"PyTorch / {device.upper()}"],
              [elapsed_cpu, elapsed_gpu])
ax.set_ylabel("Wall-clock time (seconds)")
ax.set_title(
    f"Pairwise cosine similarity ({len(all_texts)} docs)\nSpeedup: {speedup:.1f}x"
)
ax.bar_label(bars, fmt="%.2fs")
plt.tight_layout()
plt.show()

# %%
# Semantic search: top-K retrieval
# --------------------------------
# Given a query document, find the most similar documents in the corpus.

sim_np = np.asarray(sim_gpu.cpu()) if hasattr(sim_gpu, "cpu") else np.asarray(sim_gpu)

query_idx = 42
query_sims = sim_np[query_idx]
top_k = 5
top_indices = np.argsort(query_sims)[::-1][1:top_k + 1]

print(f"Query (index {query_idx}, group: {target_names[all_targets[query_idx]]}):")
print(textwrap.shorten(all_texts[query_idx].replace("\n", " "), width=120))
print()
for rank, idx in enumerate(top_indices, 1):
    print(f"  #{rank} (sim={query_sims[idx]:.3f}, group: {target_names[all_targets[idx]]})")
    print(f"       {textwrap.shorten(all_texts[idx].replace(chr(10), ' '), width=100)}")

# %%
# Near-duplicate detection
# ------------------------
# Pairs with cosine similarity above a threshold are likely near-duplicates
# or very closely related posts.

threshold = 0.90
upper_tri = np.triu(sim_np, k=1)
dup_rows, dup_cols = np.where(upper_tri > threshold)
print(f"\nPairs with cosine similarity > {threshold}: {len(dup_rows)}")

for i in range(min(3, len(dup_rows))):
    r, c = dup_rows[i], dup_cols[i]
    print(f"\n  Pair (sim={sim_np[r, c]:.3f}):")
    print(f"    Doc {r} [{target_names[all_targets[r]]}]: "
          f"{textwrap.shorten(all_texts[r].replace(chr(10), ' '), width=90)}")
    print(f"    Doc {c} [{target_names[all_targets[c]]}]: "
          f"{textwrap.shorten(all_texts[c].replace(chr(10), ' '), width=90)}")

# %%
# Distribution of pairwise similarities
# --------------------------------------
# The histogram of all pairwise similarities shows that most document pairs
# have low similarity, with a long right tail of related and duplicate pairs.

upper_values = sim_np[np.triu_indices_from(sim_np, k=1)]

fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(upper_values, bins=200, density=True, alpha=0.7)
ax.axvline(threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")
ax.set_xlabel("Cosine similarity")
ax.set_ylabel("Density")
ax.set_title("Distribution of Pairwise Cosine Similarities")
ax.legend()
plt.tight_layout()
plt.show()
