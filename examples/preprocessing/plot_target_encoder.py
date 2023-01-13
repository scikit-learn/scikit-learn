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
from sklearn.pipeline import Pipeline
from sklearn.ensemble import HistGradientBoostingRegressor

prep = ColumnTransformer(
    [
        ("num", "passthrough", numerical_features),
        ("cat", "drop", categorical_features),
    ]
)

without_categories_pipe = Pipeline(
    [("prep", prep), ("hist", HistGradientBoostingRegressor(random_state=0))]
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
# Since the categorical features have missing values, we impute the feature
# with `'sk_missing'` before passing it to the :class:`OrdinalEncoder`.
from sklearn.preprocessing import OrdinalEncoder

categorical_preprocessor = OrdinalEncoder(
    handle_unknown="use_encoded_value", unknown_value=-1
)

# %%
# We modify the original pipeline to use the ordinal categorical preprocessing:
ordinal_pipe = without_categories_pipe.set_params(prep__cat=categorical_preprocessor)
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
# --------------------------------
# Finally, we replace the ordinal encoder with the :class:`TargetEncoder`:
from sklearn.preprocessing import TargetEncoder

target_pipe = ordinal_pipe.set_params(prep__cat=TargetEncoder(target_type="continuous"))
target_pipe

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
