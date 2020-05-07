"""
==================================================================
Principal Component Regression vs Partial Least Squares Regression
==================================================================

This example compares `Principal Component Regression
<https://en.wikipedia.org/wiki/Principal_component_regression>`_ (PCR) and
`Partial Least Squares Regression
<https://en.wikipedia.org/wiki/Partial_least_squares_regression>`_ (PLS) on a
toy dataset. Our goal is to illustrate how PLS can outperform PCR when the
target is strongly correlated with some features that have a low variance.

PCR simply consists in applying :class:`~sklearn.decomposition.PCA` to the
training data (possibly performing dimensionality reduction), and then
training a regressor on the transformed training samples. In
:class:`~sklearn.decomposition.PCA`, the transformation is purely
unsupervised, meaning that no information about the targets is used. As a
result, PCR may perform poorly in some datasets where the target is
correlated with features that have a low variance. Indeed, the dimensionality
reduction of PCA tries to only keep the features that have a high variance.
Those that have a low variance will be dropped, and the final regressor will
not be able to leverage these features (or, strictly speaking, linear
combinations thereof).

PLS is both a transformer and a regressor, and it is quite similar to PCR: it
also applies a dimensionality reduction to the samples before applying a
linear regressor to the transformed data. The main difference with PCR is
that the transformation is supervised. Therefore, as we will see in this
example, it does not suffer from the issue we just mentioned.
"""

print(__doc__)

##############################################################################
# The data
# --------
#
# We start by creating a simple dataset with two features. The target `y` is
# strongly correlated with the second feature, which has a low variance. The
# first feature has a higher variance but no predictive power on the target.

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
# component against the target.
#
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

fig, axes = plt.subplots(1, 2, figsize=(10, 3))
axes[0].scatter(pca.transform(X_test), y_test)
axes[0].set(xlabel='first component', ylabel='y', title='PCR / PCA')
axes[1].scatter(pls.transform(X_test), y_test)
axes[1].set(xlabel='first component', ylabel='y', title='PLS')
plt.tight_layout()
plt.show()

##############################################################################
# As expected, the unsupervized PCA transformation of PCR has dropped the
# second feature because it has a low variance, despite it being the most
# predictive feature. This results in the first component having a low
# predictive power on the target. On the other hand, the PLS regressor manages
# to capture the effect of the second feature thanks to its use of target
# information during the transformation: it can recogize that the second
# feature should not be dropped.
#
# We also print the R-squared scores of both estimators, which further confirms
# that PLS is a better alternative than PCR in this case:

print(f"PCR r-squared {pcr.score(X_test, y_test):.3f}")
print(f"PLS r-squared {pls.score(X_test, y_test):.3f}")

##############################################################################
# As a final remark, we note that PCR with 2 components performs well, since
# the second feature wasn't dropped in this case:

pca_2 = make_pipeline(PCA(n_components=2), LinearRegression())
pca_2.fit(X_train, y_train)
print(f"PCR r-squared with 2 components {pca_2.score(X_test, y_test):.3f}")
