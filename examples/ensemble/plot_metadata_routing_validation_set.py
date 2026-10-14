"""
=================================================================================
Routing a validation set through a Pipeline to HistGradientBoosting
=================================================================================

.. currentmodule:: sklearn

This example shows how to pass a custom validation set (`X_val`, `y_val`) for early
stopping of :class:`~ensemble.HistGradientBoostingRegressor` through a
:class:`~pipeline.Pipeline` and a :class:`~model_selection.GridSearchCV`.

The pipeline encodes categorical features with :class:`~preprocessing.TargetEncoder` and
leaves numerical features unchanged. With metadata routing and auto-requests enabled,
`X_val` and `y_val` are requested by the gradient boosting estimator automatically. The
pipeline's default `transform_input=("X_val",)` transforms `X_val` with the same steps
as `X` before it reaches the final estimator.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# We enable metadata routing and auto-requests, load a small subset of the Ames
# housing data, and hold out a validation set for early stopping:
from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_openml
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import TargetEncoder

set_config(enable_metadata_routing=True, enable_metadata_auto_requests=True)

ames = fetch_openml(name="house_prices", as_frame=True, parser="auto")
categorical_features = ["Neighborhood", "HouseStyle", "ExterQual"]
numerical_features = ["OverallQual", "GrLivArea", "GarageCars"]
X = ames.data[categorical_features + numerical_features]
y = ames.target.astype(float)

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
# Number of trees after early stopping on the routed validation set:
print(f"n_iter_: {search.best_estimator_.named_steps['hist'].n_iter_}")

# %%
# Finally, we disable both configuration flags:
set_config(enable_metadata_routing=False, enable_metadata_auto_requests=False)
