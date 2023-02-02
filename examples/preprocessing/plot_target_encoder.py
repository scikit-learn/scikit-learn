"""
============================================
Comparing Target Encoder with Other Encoders
============================================

.. currentmodule:: sklearn.preprocessing

The :class:`TargetEncoder` uses the value of the target to encode each
categorical feature. In this example, we will compare three different approaches
for handling categorical features: :class:`TargetEncoder`,
:class:`OrdinalEncoder`, :class:`OneHotEncoder` and dropping the category.
For more information about the :class:`TargetEncoder` see the
:ref:`User Guide <target_encoder>`.

.. note::
  :class:`TargetEncoder` uses a cross validation scheme in
  :meth:`~TargetEncoder.fit_transform` to prevent leaking the target during training.
  In :meth:`~TargetEncoder.fit_transform`, the data is split according to the `cv`
  parameter. Categorical encodings are learned from split and used to encode the other
  split. Afterwards, a final categorical encoding is learned from all the data, which
  is then used to encode data during :meth:`~TargetEncoder.transform`. This means that
  `fit(X, y).transform(X)` does not equal `fit_transform(X, y)`.
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
numerical_features = ["price"]
categorical_features = [
    "country",
    "province",
    "region_1",
    "region_2",
    "variety",
    "winery",
]
target_name = "points"

X = df[numerical_features + categorical_features]
y = df[target_name]

_ = y.hist()

# %%
# Training and Evaluating Pipelines with Different Encoders
# =========================================================
# In this section, we will evaluate pipelines with
# :class:`~sklearn.ensemble.HistGradientBoostingRegressor` with different encoding
# strategies. First, we list out the encoders we will be using to preprocess
# the categorical features:
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import TargetEncoder

categorical_preprocessors = [
    ("drop", "drop"),
    ("ordinal", OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1)),
    (
        "one_hot",
        OneHotEncoder(handle_unknown="ignore", max_categories=10, sparse_output=False),
    ),
    ("target", TargetEncoder(target_type="continuous")),
]

# %%
# Next, we evaluate the models using cross validation and record the results:
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import HistGradientBoostingRegressor

results = []
n_cv_folds = 3
max_iter = 20
for name, categorical_preprocessor in categorical_preprocessors:
    preprocessor = ColumnTransformer(
        [
            ("numerical", "passthrough", numerical_features),
            ("categorical", categorical_preprocessor, categorical_features),
        ]
    )
    pipe = make_pipeline(
        preprocessor, HistGradientBoostingRegressor(random_state=0, max_iter=max_iter)
    )
    rmse_scores = -cross_val_score(
        pipe, X, y, scoring="neg_root_mean_squared_error", cv=n_cv_folds
    )
    results.append(
        {
            "preprocessor": name,
            "rmse_mean": rmse_scores.mean(),
            "rmse_std": rmse_scores.std(),
        }
    )

# %%
# Native Categorical Feature Support
# ----------------------------------
# In this section, we build and evaluate a pipeline that uses native categorical
# feature support in :class:`~sklearn.ensemble.HistGradientBoostingRegressor`,
# which only supports up to 255 unique categories. In our dataset, the most of
# the categorical features have more than 255 unique categories:
n_unique_categories = df[categorical_features].nunique().sort_values(ascending=False)

# %%
# To workaround the limitation above, we group the categorical features into
# low cardinality and high cardinality features. The high cardinality features
# will be target encoded and the low cardinality features will use the native
# categorical feature in gradient boosting.
high_cardinality_features = n_unique_categories[n_unique_categories > 255].index
low_cardinality_features = n_unique_categories[n_unique_categories <= 255].index
mixed_encoded_preprocessor = ColumnTransformer(
    [
        ("numerical", "passthrough", numerical_features),
        (
            "high_cardinality",
            TargetEncoder(target_type="continuous"),
            high_cardinality_features,
        ),
        (
            "low_cardinality",
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1),
            low_cardinality_features,
        ),
    ],
    verbose_feature_names_out=False,
)

# The output of the of the preprocessor must be set to pandas so the
# gradient boosting model can detect the low cardinality features.
mixed_encoded_preprocessor.set_output(transform="pandas")
mixed_pipe = make_pipeline(
    mixed_encoded_preprocessor,
    HistGradientBoostingRegressor(
        random_state=0, max_iter=max_iter, categorical_features=low_cardinality_features
    ),
)
mixed_pipe

# %%
# Finally, we evaluate the pipeline using cross validation and record the results:
rmse_scores = -cross_val_score(
    mixed_pipe, X, y, scoring="neg_root_mean_squared_error", cv=n_cv_folds
)
results.append(
    {
        "preprocessor": "mixed_target_encoder",
        "rmse_mean": rmse_scores.mean(),
        "rmse_std": rmse_scores.std(),
    }
)


# %%
# Plotting the Results
# --------------------
# In this section, we display the results for all the encoders. When evaluating the
# predictive performance on the test set, dropping the categories perform the
# worst and the target encoder performs the best. The target encoding provides
# more information about the target, which the gradient boosted model at the end
# of the pipeline can take advantage of.
import matplotlib.pyplot as plt
import pandas as pd

results_df = (
    pd.DataFrame(results)
    .set_index("preprocessor")
    .sort_values("rmse_mean", ascending=False)
)

fig, ax = plt.subplots()
yticks = range(len(results_df))
ax.barh(
    y=yticks,
    width=results_df["rmse_mean"],
    xerr=results_df["rmse_std"],
    height=0.9,
    color=["C0", "C1", "C2", "C3", "C4"],
)
ax.set(
    xlabel="Root Mean Squared Error",
    ylabel="Model",
    yticks=yticks,
    yticklabels=results_df.index,
)
