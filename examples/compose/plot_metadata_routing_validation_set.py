"""
============================================================================
Routing a validation set through a Pipeline to HistGradientBoostingRegressor
============================================================================

.. currentmodule:: sklearn

This example shows how to pass a validation set (`X_val`, `y_val`) for early stopping of
:class:`~ensemble.HistGradientBoostingRegressor` through a :class:`~pipeline.Pipeline`
and a :class:`~model_selection.GridSearchCV`.

Estimators such as :class:`~ensemble.HistGradientBoostingRegressor` can use a validation
set inside :term:`fit` for early stopping. If that set is split off from `X` after
preprocessing fitted on the full training set, information from the validation set can
leak into the transformers (for example target statistics in
:class:`~preprocessing.TargetEncoder`). A safer pattern is to hold out validation data
and then run it through the same pipeline steps as `X`. This was however not possible in
scikit-learn versions < 1.7.

Today we can use :ref:`Metadata routing <metadata_routing>` to transform selected
metadata listed in the `transform_input` parameter of :class:`~pipeline.Pipeline` (by
default `X_val`) until a later step of the pipeline consumes it. With auto-requests
enabled, the gradient boosting estimator requests `X_val` and `y_val` without an
explicit `set_fit_request`.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# We enable metadata routing and auto-requests, load a subset of the Ames housing
# data, and hold out a validation set for early stopping *before* fitting the pipeline:
from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_openml
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import TargetEncoder

set_config(enable_metadata_routing=True, enable_metadata_auto_requests=True)

ames = fetch_openml(data_id=42165, as_frame=True, parser="auto")
categorical_features = ["Neighborhood", "HouseStyle", "ExterQual"]
numerical_features = ["OverallQual", "GrLivArea", "GarageCars"]
X = ames.data[categorical_features + numerical_features]
y = ames.target

X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=0)

# %%
# Now we build a pipeline and search over a small hyperparameter grid. `X_val` is
# transformed by the preprocessor because it is listed in `transform_input` (the
# default). No `set_fit_request` calls are needed for `X_val` / `y_val` when
# auto-requests are enabled:
preprocessor = ColumnTransformer(
    transformers=[
        (
            "cat",
            TargetEncoder(target_type="continuous"),
            categorical_features,
        ),
    ],
    remainder="passthrough",
)
hist = HistGradientBoostingRegressor(
    early_stopping=True,
    max_iter=50,
    random_state=0,
)
pipe = Pipeline([("prep", preprocessor), ("hist", hist)])
search = GridSearchCV(
    pipe,
    param_grid={"hist__max_depth": [2, 3]},
    cv=2,
)
search.fit(X_train, y_train, X_val=X_val, y_val=y_val)
print(search.best_params_)
print(f"best CV score: {search.best_score_:.3f}")
# Number of trees after early stopping based on the routed validation set:
print(f"n_iter_: {search.best_estimator_.named_steps['hist'].n_iter_}")

# %%
# Finally, we disable both configuration flags:
set_config(enable_metadata_routing=False, enable_metadata_auto_requests=False)
