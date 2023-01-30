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
# Then, we split the dataset into a training and test set.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

print(
    f"Data points in training set: {len(X_train)}\n"
    f"Data points in test set: {len(X_test)}"
)

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
        OneHotEncoder(handle_unknown="ignore", max_categories=50, sparse_output=False),
    ),
    ("target", TargetEncoder(target_type="continuous")),
]

# %%
# Next, we create and fit the pipelines to the training data and evaluate them on the
# test set:
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error

results = []
for name, categorical_preprocessor in categorical_preprocessors:
    preprocessor = ColumnTransformer(
        [
            ("numerical", "passthrough", numerical_features),
            ("categorical", categorical_preprocessor, categorical_features),
        ]
    )
    pipe = make_pipeline(preprocessor, HistGradientBoostingRegressor(random_state=0))
    pipe.fit(X_train, y_train)

    y_pred = pipe.predict(X_test)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    results.append({"preprocessor": name, "root mean squared error": rmse})

# %%
# Finally, we display the results for all the encoders. When evaluating the
# predictive performance on the test set, dropping the categories perform the
# worst and the target encoder performs the best. The target encoding provides
# more information about the target, which the gradient boosted model at the end
# of the pipeline can take advantage of. In this example, we evaluate with a single
# test split to reduce the runtime. We recommend using cross-validation to check
# the significance of the improvement.
import pandas as pd

_ = pd.DataFrame(results).sort_values("root mean squared error", ascending=False).set_index(
    "preprocessor"
).plot(kind="barh")
