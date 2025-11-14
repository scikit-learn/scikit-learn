# -*- coding: utf-8 -*-
"""
KMeans: required shapes and dtypes
==================================

This example shows:
- The required 2D input shape (n_samples, n_features).
- How integer arrays are cast to float.
- The minor numerical differences between float64 and float32.
- The common 1D input error and the correct reshape.

The example prints small, deterministic summaries to keep execution fast
for the documentation build.
"""

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


# Reproducible toy data: (n_samples, n_features)
rng = np.random.RandomState(42)
X_int = rng.randint(0, 100, size=(300, 2))  # integers; will be cast to float

# Standardize to improve clustering stability; keep both float64/float32 views
X64 = StandardScaler().fit_transform(X_int.astype(np.float64))
X32 = X64.astype(np.float32)

km_kwargs = dict(n_clusters=3, n_init="auto", random_state=0)

# Fit on float64
km64 = KMeans(**km_kwargs).fit(X64)
sil64 = silhouette_score(X64, km64.labels_)

# Fit on float32
km32 = KMeans(**km_kwargs).fit(X32)
sil32 = silhouette_score(X32, km32.labels_)

print(f"X64 dtype: {X64.dtype}, inertia: {km64.inertia_:.3f}, silhouette: {sil64:.3f}")
print(f"X32 dtype: {X32.dtype}, inertia: {km32.inertia_:.3f}, silhouette: {sil32:.3f}")

# Demonstrate the common 1D input error and the fix
x_1d = rng.rand(10)  # shape (10,)
try:
    KMeans(n_clusters=2).fit(x_1d)  # expected to raise an informative error
except Exception as e:
    msg = str(e).splitlines()[0]
    print(f"1D input error: {type(e).__name__}: {msg}")

# Correct 2D shape: (n_samples, 1)
x_col = x_1d.reshape(-1, 1)
km1 = KMeans(n_clusters=2, n_init="auto", random_state=0).fit(x_col)
print(f"Reshaped OK; inertia: {km1.inertia_:.3f}, centers shape: {km1.cluster_centers_.shape}")
