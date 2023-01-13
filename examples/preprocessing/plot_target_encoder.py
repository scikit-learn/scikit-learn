"""
============================================
Comparing Target Encoder with Other Encoders
============================================

.. currentmodule:: sklearn.preprocessing

The :class:`TargetEncoder` uses the value of the target to encode each
categorical feature. In this example, we will compare three different approaches
for handling categorical features: :class:`TargetEncoder`,
:class:`OrdinalEncoder`, and dropping the category.
"""

# %%
# Loading Data from OpenML
# ========================
# First, we load the wine reviews dataset, where the target is the points given
# be a reviewer:
from sklearn.datasets import fetch_openml

wine_reviews = fetch_openml(data_id=42074, as_frame=True, parser="liac-arff")

df = wine_reviews.frame
df.head()

# %%
# For this example, we use the following subset of numerical and categorical
# features in the data. The target are continuous values from 80 to 100:
import numpy as np

numerical_features = ["price"]
categorical_features = [
    "country",
    "province",
    "region_1",
    "region_2",
    "variety",
    "winery",
]

X = df[numerical_features + categorical_features]
y = df["points"]

np.sort(y.unique())

# %%
# Then, we split the dataset into a training and test set.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

print(f"Samples in training set: {len(X_train)}\nSamples in test set: {len(X_test)}")

# %%
# Building and Training Pipelines with Different Encoders
# =======================================================
# Dropping the categorical features
# ---------------------------------
# As a baseline, we construct a pipeline where the categorical features are
# dropped.
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import HistGradientBoostingRegressor

preprocessor = ColumnTransformer(
    [
        ("num", "passthrough", numerical_features),
        ("cat", "drop", categorical_features),
    ]
)

without_categories_pipe = make_pipeline(
    preprocessor, HistGradientBoostingRegressor(random_state=0)
)
without_categories_pipe

# %%
# Here we train and evaluated on the test set to get a baseline metric when
# the categories are dropped:
from sklearn.metrics import mean_squared_error

without_categories_pipe.fit(X_train, y_train)
without_categories_rsme = mean_squared_error(
    y_test, without_categories_pipe.predict(X_test), squared=False
)
print(f"RMSE for dropping categorical features: {without_categories_rsme:.4}")

# %%
# Using the OrdinalEncoder
# ------------------------
# We create a pipeline using the ordinal categorical preprocessor:
from sklearn.preprocessing import OrdinalEncoder

ordinal_preprocessor = ColumnTransformer(
    [
        ("num", "passthrough", numerical_features),
        (
            "cat",
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1),
            categorical_features,
        ),
    ]
)

ordinal_pipe = make_pipeline(
    ordinal_preprocessor, HistGradientBoostingRegressor(random_state=0)
)
ordinal_pipe

# %%
# When we include the categorical features through ordinal encoding the model improves
# when evaluated with the test set:
ordinal_pipe.fit(X_train, y_train)
ordinal_pipe_rmse = mean_squared_error(
    y_test, ordinal_pipe.predict(X_test), squared=False
)
print(f"RMSE with ordinal encoding: {ordinal_pipe_rmse:.4}")

# %%
# Using the TargetEncoder
# -----------------------
# Finally, we create a pipeline with the :class:`TargetEncoder`:
from sklearn.preprocessing import TargetEncoder

target_preprocessor = ColumnTransformer(
    [
        ("num", "passthrough", numerical_features),
        (
            "cat",
            TargetEncoder(target_type="continuous"),
            categorical_features,
        ),
    ]
)

target_pipe = make_pipeline(
    target_preprocessor, HistGradientBoostingRegressor(random_state=0)
)

# %%
# When the model is evalute on the test set, we see that the
# :class:`TargetEncoder` further improves the predictive performance of the
# model. The target encoding provides more information about the target, which
# the regression model at the end of the pipeline can take advantage of.
target_pipe.fit(X_train, y_train)
target_pipe_rmse = mean_squared_error(
    y_test, target_pipe.predict(X_test), squared=False
)
print(f"RMSE with target encoding: {target_pipe_rmse:.4}")

# %%
from sklearn.preprocessing import OneHotEncoder

ohe_preprocessor = ColumnTransformer(
    [
        ("num", "passthrough", numerical_features),
        (
            "cat",
            OneHotEncoder(
                handle_unknown="ignore", max_categories=50, sparse_output=False
            ),
            categorical_features,
        ),
    ]
)

ohe_pipe = make_pipeline(
    ohe_preprocessor, HistGradientBoostingRegressor(random_state=0)
)
ohe_pipe.fit(X_train, y_train)
ohe_pipe_rmse = mean_squared_error(y_test, ohe_pipe.predict(X_test), squared=False)
print(f"RMSE with OHE encoding: {ohe_pipe_rmse:.4}")
