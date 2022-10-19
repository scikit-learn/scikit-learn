"""
================================
Introducing the `set_output` API
================================

.. currentmodule:: sklearn

This example will demonstrate the `set_output` API to configure transformers to
output pandas DataFrames. `set_output` can be configured per estimator by calling
the `set_output` method or globally by setting `set_config(transform_output="pandas")`.
For details, see
`SLEP018 <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep018/proposal.html>`__.
"""  # noqa

# %%
# First, we load the iris dataset as a DataFrame to demonstrate the `set_output` API.
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split

X, y = load_iris(as_frame=True, return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=0)
X_train.head()

# %%
# To configure an estimator such as :class:`preprocessing.StandardScaler` to return
# DataFrames, call `set_output`. This feature requires pandas to be installed.

from sklearn.preprocessing import StandardScaler

scaler = StandardScaler().set_output(transform="pandas")

scaler.fit(X_train)
X_test_scaled = scaler.transform(X_test)
X_test_scaled.head()

# %%
# `set_output` can be called after `fit` to configure `transform` after the fact.
scaler2 = StandardScaler()

scaler2.fit(X_train)
X_test_np = scaler2.transform(X_test)
print(f"Default output type: {type(X_test_np).__name__}")

scaler2.set_output(transform="pandas")
X_test_df = scaler2.transform(X_test)
print(f"Configured pandas output type: {type(X_test_df).__name__}")

# %%
# In a :class:`pipeline.Pipeline`, `set_output` configures all steps to output
# DataFrames.
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectPercentile

clf = make_pipeline(
    StandardScaler(), SelectPercentile(percentile=75), LogisticRegression()
)
clf.set_output(transform="pandas")
clf.fit(X_train, y_train)

# %%
# Each transformer in the pipeline is configured to return DataFrames. This
# means that the final logistic regression step contains the feature names of the input.
clf[-1].feature_names_in_

# %%
# Next we load the titanic dataset to demonstrate `set_output` with
# :class:`compose.ColumnTransformer` and heterogenous data.
from sklearn.datasets import fetch_openml

X, y = fetch_openml(
    "titanic", version=1, as_frame=True, return_X_y=True, parser="pandas"
)
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y)

# %%
# The `set_output` API can be configured globally by using :func:`set_config` and
# setting `transform_output` to `"pandas"`.
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn import set_config

set_config(transform_output="pandas")

num_pipe = make_pipeline(SimpleImputer(), StandardScaler())
ct = ColumnTransformer(
    (
        ("numerical", num_pipe, ["age", "fare"]),
        (
            "categorical",
            OneHotEncoder(
                sparse_output=False, drop="if_binary", handle_unknown="ignore"
            ),
            ["embarked", "sex", "pclass"],
        ),
    ),
    verbose_feature_names_out=False,
)
clf = make_pipeline(ct, SelectPercentile(percentile=50), LogisticRegression())
clf.fit(X_train, y_train)
clf.score(X_test, y_test)

# %%
# With the global configuration, all transformers output DataFrames. This allows us to
# easily plot the logistic regression coefficients with the corresponding feature names.
import pandas as pd

log_reg = clf[-1]
coef = pd.Series(log_reg.coef_.ravel(), index=log_reg.feature_names_in_)
_ = coef.sort_values().plot.barh()

# %%
# This resets `transform_output` to its default value to avoid impacting other
# examples when generating the scikit-learn documentation
set_config(transform_output="default")
