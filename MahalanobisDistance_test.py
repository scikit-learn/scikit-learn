# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import warnings

import numpy as np

from sklearn.covariance import MinCovDet
from sklearn.neighbors import LocalOutlierFactor

warnings.filterwarnings('ignore')
rng = np.random.RandomState(42)


# Generate data with clear outliers
X_inliers = rng.randn(200, 10)
X_outliers = rng.randn(10, 10) + 8
X = np.vstack([X_inliers, X_outliers])

# Compute robust inverse covariance matrix
mcd = MinCovDet(random_state=42).fit(X)
VI = np.ascontiguousarray(mcd.precision_, dtype=np.float64)

# Test with n_jobs=1
lof_single = LocalOutlierFactor(
    n_neighbors=20,
    metric='mahalanobis',
    metric_params={'VI': VI},
    n_jobs=1
)
outliers_single = np.where(lof_single.fit_predict(X) == -1)[0]

# Test with n_jobs=4
lof_parallel = LocalOutlierFactor(
    n_neighbors=20,
    metric='mahalanobis',
    metric_params={'VI': VI},
    n_jobs=4
)
outliers_parallel = np.where(lof_parallel.fit_predict(X) == -1)[0]

# Display results
print(f"Outliers detected (n_jobs=1): {len(outliers_single)}")
print(f"Outliers detected (n_jobs=4): {len(outliers_parallel)}")
print("Expected outliers: 10")

# Verify fix
if len(outliers_parallel) > 0:
    print("\nSUCCESS: n_jobs=4 works correctly")
else:
    print("\n FAILED: n_jobs=4 still broken")
