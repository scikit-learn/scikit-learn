"""
=============================
Target Encoder for Regressors
=============================

.. currentmodule:: sklearn.preprocessing

The :class:`TargetRegressorEncoder` uses target statistics conditioned on
the categorical features for encoding. In this example, we will compare
:class:`TargetRegressorEncoder`, :class:`OrdinalEncoder`, and dropping the
category on a wine review dataset.
"""

# %%
# Loading Data from OpenML
# ========================
# First, we load the wine reviews dataset, where the target is the points given
# be a reviewer:
import warnings
from sklearn.datasets import fetch_openml

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    wine_reviews = fetch_openml(data_id=42074, as_frame=True)

df = wine_reviews.frame
df.head()

# %%
# For this example, we use the following subset of numerical and categorical
# features in the data. The categorical features have a cardinality ranging
# from 18 to 14810:
numerical_features = ['price']
categorical_features = ['country', 'province', 'region_1', 'region_2',
                        'variety', 'winery']

X = df[numerical_features + categorical_features]
y = df['points']
X.nunique().sort_values(ascending=False)

# %%
# We split the dataset into a training and test set:
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

print(f"Samples in training set: {len(X_train)}\n"
      f"Samples in test set: {len(X_test)}")

# %%
# Building and Training Pipelines with Different Encoders
# =======================================================
# Dropping the categorical features
# ---------------------------------
# As a basline, we construct a pipeline where the categorical features are
# dropped.
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer

set_config(display='diagram')  # Show HTML representation of pipeline

prep = ColumnTransformer([
    ('num', SimpleImputer(strategy='median'), numerical_features),
    ('cat', 'drop', categorical_features)
])

reg_drop_cats = Pipeline([
    ('prep', prep), ('hist', HistGradientBoostingRegressor())
])
reg_drop_cats

# %%
# Here we train and use the root mean squared error to evalute the baseline
# model:
from sklearn.metrics import mean_squared_error

reg_drop_cats.fit(X_train, y_train)
reg_drop_cats_rmse = mean_squared_error(
    y_test, reg_drop_cats.predict(X_test)
)
print(f"RMSE for dropping categorical features: {reg_drop_cats_rmse:.4}")

# %%
# Using the OrdinalEncoder
# ------------------------
# Since the categorical features have missing values, we impute the feature
# with `'sk_missing'` before passing it to the :class:`OrdinalEncoder`.
# The `categories` parameter is constructed such that there are not unknown
# values during test time
from sklearn.preprocessing import OrdinalEncoder

categories = [
    X[feat].fillna("sk_missing").unique() for feat in categorical_features
]

cat_prep = Pipeline([
    ('imputer', SimpleImputer(strategy='constant', missing_values=None,
                              fill_value='sk_missing')),
    ('encoder', OrdinalEncoder(categories=categories))
])


# %%
# We modify the original pipeline to use the ordinal categorical preprocessing:
reg_ordinal = reg_drop_cats.set_params(prep__cat=cat_prep)
reg_ordinal

# %%
# When we include the categorical features through ordinal encoding the RMSE
# improves:
reg_ordinal.fit(X_train, y_train)
reg_ordinal_rmse = mean_squared_error(
    y_test, reg_ordinal.predict(X_test), squared=False
)
print(f"RMSE with ordinal encoding: {reg_ordinal_rmse:.4}")

# %%
# Using the TargetEncoder
# -----------------------
# Finally, we replace the ordinal encoder with the
# :class:`TargetRegressorEncoder`:
from sklearn.preprocessing import TargetRegressorEncoder

reg_target = reg_ordinal.set_params(
    prep__cat__encoder=TargetRegressorEncoder())
reg_target

# %%
# The :class:`TargetRegressorEncoder` further improves the RMSE:
reg_target.fit(X_train, y_train)
reg_target_rmse = mean_squared_error(
    y_test, reg_target.predict(X_test), squared=False
)
print(f"RMSE with target encoding: {reg_target_rmse:.4}")
