# flake8: noqa
"""
=======================================
Release Highlights for scikit-learn 1.1
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.1! Many bug fixes
and improvements were added, as well as some new key features. We detail
below a few of the major features of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <changes_1_1>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Quantile loss in :class:`ensemble.HistGradientBoostingRegressor`
# ----------------------------------------------------------------
# :class:`ensemble.HistGradientBoostingRegressor` can model quantiles with
# `loss="quantile"` and the new parameter `quantile`.
from sklearn.datasets import make_regression
from sklearn.ensemble import HistGradientBoostingRegressor
import pandas as pd

X, y = make_regression(n_features=4, n_informative=8, noise=10, random_state=0)

quantiles = [0.05, 0.5, 0.95]
hist_quantiles = [
    HistGradientBoostingRegressor(loss="quantile", quantile=quantile).fit(X, y)
    for quantile in quantiles
]
predictions = {
    f"quantile_{quantile:0.2f}": hist.predict(X)
    for quantile, hist in zip(quantiles, hist_quantiles)
}
pd.DataFrame(predictions).iloc[:5]

# %%
# `get_feature_names_out` Avaliable in All Transformers
# -----------------------------------------------------
# :term:`get_feature_names_out` is now avaliable in all Transformers. This enables
# :class:`pipeline.Pipeline` to construct the feature names out for more comlex
# pipelines:
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest
from sklearn.datasets import fetch_openml
from sklearn.linear_model import LogisticRegression

X, y = fetch_openml("titanic", version=1, as_frame=True, return_X_y=True)
numeric_features = ["age", "fare"]
numeric_transformer = make_pipeline(SimpleImputer(strategy="median"), StandardScaler())
categorical_features = ["embarked", "pclass"]

preprocessor = ColumnTransformer(
    [
        ("num", numeric_transformer, numeric_features),
        (
            "cat",
            OneHotEncoder(handle_unknown="ignore", sparse=False),
            categorical_features,
        ),
    ],
    verbose_feature_names_out=False,
)
log_reg = make_pipeline(preprocessor, SelectKBest(k=7), LogisticRegression())
log_reg.fit(X, y)


# %%
# Here we slice the pipeline to include all the steps but the last one. The output
# feature names of this pipeline slice is the features inputted into logistic
# regression. These names directly corresponds to the coefficients in the logistic
# regression:
log_reg_input_features = log_reg[:-1].get_feature_names_out()
pd.Series(log_reg[-1].coef_.ravel(), index=log_reg_input_features)


# %%
# Grouping infrequent categories in :class:`OneHotEncoder`
# --------------------------------------------------------
# :class:`OneHotEncoder` supports aggregating infrequent categories into a single
# output for each feature. The parameters to enable the gathering of infrequent
# categories are `min_frequency` and `max_categories`. See the
# :ref:`User Guide <one_hot_encoder_infrequent_categories>` for more details.
from sklearn.preprocessing import OneHotEncoder
import numpy as np

X = np.array(
    [["dog"] * 5 + ["cat"] * 20 + ["rabbit"] * 10 + ["snake"] * 3], dtype=object
).T
enc = OneHotEncoder(min_frequency=6, sparse=False).fit(X)
enc.infrequent_categories_

# %%
# Since dog and snake are infrequent categories, they are grouped together when
# transformed:
encoded = enc.transform(np.array([["dog"], ["snake"], ["cat"], ["rabbit"]]))
pd.DataFrame(encoded, columns=enc.get_feature_names_out())

# %%
# Performance Improvements
# ------------------------
# Reductions on pairwise distances for dense float64 datasets has been refactored
# to better take advantage of parallelism. For example,
# :meth:`neighbors.NearestNeighbors.kneighbors` and
# :meth:`neighbors.NearestNeighbors.radius_neighbors` can respectively be up to ×20 and
# ×5 faster than previously. In summary, the following functions and estimators
# now benefit from improved performances:
#
# - :func:`metrics.pairwise_distances_argmin`
# - :func:`metrics.pairwise_distances_argmin_min`
# - :class:`cluster.AffinityPropagation`
# - :class:`cluster.Birch`
# - :class:`cluster.MeanShift`
# - :class:`cluster.OPTICS`
# - :class:`cluster.SpectralClustering`
# - :func:`feature_selection.mutual_info_regression`
# - :class:`neighbors.KNeighborsClassifier`
# - :class:`neighbors.KNeighborsRegressor`
# - :class:`neighbors.RadiusNeighborsClassifier`
# - :class:`neighbors.RadiusNeighborsRegressor`
# - :class:`neighbors.LocalOutlierFactor`
# - :class:`neighbors.NearestNeighbors`
# - :class:`manifold.Isomap`
# - :class:`manifold.LocallyLinearEmbedding`
# - :class:`manifold.TSNE`
# - :func:`manifold.trustworthiness`
# - :class:`semi_supervised.LabelPropagation`
# - :class:`semi_supervised.LabelSpreading`
#
# The computation of loss has been refactored using Cython resulting in performance
# improvements in the following estimators:
#
# - :class:`linear_model.LogisticRegression`
# - :class:`linear_model.GammaRegressor`
# - :class:`linear_model.PoissonRegressor`
# - :class:`linear_model.TweedieRegressor`
