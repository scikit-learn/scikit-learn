"""
==========================================================================
Fitting an Elastic Net with a precomputed Gram Matrix and Weighted Samples
==========================================================================

The following example shows how to precompute the gram matrix
while using weighted samples with an ElasticNet.

If weighted samples are used, the design matrix must be centered and then
rescaled by the sqrt of the weight vector before the gram matrix
is computed. Note: sample_weight vector is also rescaled to sum to
n_samples, see comment in :func:`linear_model.ElasticNet.fit`.

"""

print(__doc__)

import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.datasets import make_regression


X, y = make_regression(n_samples=int(1e5), noise=0.5)

sample_weight = np.random.lognormal(size=y.shape)

w_norm = sample_weight * (y.shape / np.sum(sample_weight))
X_offset = np.average(X, axis=0, weights=w_norm)
X_c = (X - np.average(X, axis=0, weights=w_norm))
X_r = X_c * np.sqrt(w_norm)[:, np.newaxis]
gram = np.dot(X_r.T, X_r)

lm = ElasticNet(alpha=0.01, precompute=gram)

# We must use the centered design matrix for fitting here. If we use the
# original matrix, the preprocessing code will detect that it is uncentered
# and also toss out our gram matrix. However, if we use X_r, the preprocessing
# code will incorrectly rescale the design matrix a second time.
lm.fit(X_c, y, sample_weight=sample_weight)
