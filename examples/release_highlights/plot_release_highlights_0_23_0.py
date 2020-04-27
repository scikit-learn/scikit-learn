"""
========================================
Release Highlights for scikit-learn 0.23
========================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 0.23, which comes
with many bug fixes and new features! We detail below a few of the major
features of this release. For an exhaustive list of all the changes, please
refer to the :ref:`release notes <changes_0_23>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install scikit-learn
"""

##############################################################################
# Generalized Linear Models
# -------------------------
# Long-awaited Generalized Linear Models with non-normal loss functions are now
# available. In particular, three new regressors were implemented:
# :class:`~sklearn.linear_model.PoissonRegressor`,
# :class:`~sklearn.linear_model.GammaRegressor`, and
# :class:`~sklearn.linear_model.TweedieRegressor`. The Poisson regressor can be
# used to model positive integer counts, or relative frequencies. Read more in
# the :ref:`User Guide <Generalized_linear_regression>`.

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import PoissonRegressor

n_samples, n_features = 100, 20
rng = np.random.RandomState(0)
X = rng.randn(n_samples, n_features)
y = ((10 + X[:, 5]) * 10).astype(int)  # easy target correlated with X[:, 5]
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=rng)
reg = PoissonRegressor()
reg.fit(X_train, y_train)
print(reg.score(X_test, y_test))

##############################################################################
# Scalability and stability improvements to KMeans
# ------------------------------------------------
# The :class:`~sklearn.cluster.KMeans` estimator was entirely re-worked, and it
# is now significantly faster and more stable. In addition, the Elkan algorithm
# is now compatible with sparse matrices. The estimator now uses OpenMP based
# parallelism instead of relying on joblib, so the `n_jobs` parameter has no
# effect anymore. For more details on how to control the number of threads,
# please refer to our :ref:`parallelism` notes.
import scipy
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.metrics import completeness_score

X, y = make_blobs(random_state=rng)
X = scipy.sparse.csr_matrix(X)
X_train, X_test, _, y_test = train_test_split(X, y, random_state=rng)
kmeans = KMeans(algorithm='elkan').fit(X_train)
print(completeness_score(kmeans.predict(X_test), y_test))


##############################################################################
# Improvements to the histogram-based Gradient Boosting estimators
# ----------------------------------------------------------------
# Various improvements were made to
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier` and
# :class:`~sklearn.ensemble.HistGradientBoostingRegressor`. These estimators
# now support :ref:`sample weights <sw_hgbdt>`, and a new Poisson loss function
# was implemented. In addition, users can define :ref:`monotonic constraints
# <monotonic_cst_gbdt>` to constrain the predictions based on the variations of
# specific features. Finally, an automatic early-stopping criteria was added:
# early-stopping is enabled by default when the number of samples exceeds 10k.
# The following snippet illustrates the use of the Poisson loss.

from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor

X = rng.randn(n_samples, n_features)
y = ((10 + X[:, 5]) * 10).astype(int)  # easy target correlated with X[:, 5]
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=rng)
reg = HistGradientBoostingRegressor(loss='poisson', learning_rate=.01)
reg.fit(X_train, y_train)
print(reg.score(X_test, y_test))
