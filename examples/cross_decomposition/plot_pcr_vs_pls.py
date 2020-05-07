"""
==================================================================
Principal Component Regression vs Partial Least Squares Regression
==================================================================

This example compares Principal Component Regression (PCR) and Partial Least
Squares Regression (PLS) on a toy dataset.

PCR simply consists in applying :class:`~sklearn.decomposition.PCA` to the
training data, and then train a regressor on the transformed training sample
`X`. in :class:`~sklearn.decomposition.PCA`, the transformation is purely
unsupervized, meaning that no information about the targets are used.
As a result, PCR may perform poorly in some datasets where the target is
correlated with features that have a low variance: if dimensionality
reduction is applied, these features will be dropped in the transformed data.


PLS is both a transformer and a regressor, and it is quite similar to PCR.
The main difference is that the transformation is supervised. Therefore, it
does not suffer from the issue we just mentionned.
"""

print(__doc__)

##############################################################################
# The data
# --------
#
# We start by creating a simple datasets with two features. The target `y` is
# strongly correlated with the second feature, which has a low variance.

import numpy as np
import matplotlib.pyplot as plt

rng = np.random.RandomState(0)
n_samples = 500
f0 = rng.normal(scale=10, size=n_samples)
f1 = rng.normal(scale=1, size=n_samples)
X = np.c_[f0, f1]
y = f1 + rng.normal(scale=.1, size=n_samples)

fig, axes = plt.subplots(1, 3, figsize=(10, 3))

axes[0].scatter(f0, f1)
axes[0].set(ylim=(-10, 10), xlim=(-20, 20), xlabel='f0', ylabel='f1',
            title='training data')
axes[1].scatter(f0, y)
axes[1].set(xlabel='f0', ylabel='y', title='target vs f0')
axes[2].scatter(f1, y)
axes[2].set(xlabel='f1', ylabel='y', title='target vs f1')
plt.tight_layout()
plt.show()

##############################################################################
# Projection on one component and predictive power
# ------------------------------------------------
#
# We now create two regressors: PCR and PLS, and for our illustration purposes
# we set the number of components to 1. For both models, we plot the first
# component against the target. As expected, the PCA transformation of PCR has
# dropped the second feature because it has a low variance, whichs results in
# the first component having a very low predictive power on the target. On the
# other hand, the PLS regressor manages to capture the effect of the second
# feature thanks to its use of target information during the transformation.
#
# We also print the R-squared scores of both estimators, which further confirms
# that PLS is a better alternative than PCR in this case.

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=rng)

pcr = make_pipeline(PCA(n_components=1), LinearRegression())
pcr.fit(X_train, y_train)
pca = pcr.named_steps['pca']  # retrieve the PCA step of the pipeline

pls = PLSRegression(n_components=1)
pls.fit(X_train, y_train)

fig, axes = plt.subplots(1, 2)
axes[0].scatter(pca.transform(X_test), y_test)
axes[0].set(xlabel='first component', ylabel='y', title='PCR / PCA')
axes[1].scatter(pls.transform(X_test), y_test)
axes[1].set(xlabel='first component', ylabel='y', title='PLS')
plt.tight_layout()
plt.show()
print(f"PCR r-squared {pcr.score(X_test, y_test):.3f}")
print(f"PLS r-squared {pls.score(X_test, y_test):.3f}")

##############################################################################
# As a final remark, we note that PCR with 2 components performs well, since
# the second feature wasn't dropped in this case:

pca_2 = make_pipeline(PCA(n_components=2), LinearRegression())
pca_2.fit(X_train, y_train)
print(f"PCR r-squared with 2 components {pca_2.score(X_test, y_test):.3f}")
