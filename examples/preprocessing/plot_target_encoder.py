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

X = df[numerical_features + categorical_features]
y = df["points"]

y.hist()

# %%
# Then, we split the dataset into a training and test set.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

print(f"Samples in training set: {len(X_train)}\nSamples in test set: {len(X_test)}")

# %%
# Training and Evaluating Pipelines with Different Encoders
# =========================================================
# In this section, we will evalute pipelines with
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
        "one_hot_encoder",
        OneHotEncoder(handle_unknown="ignore", max_categories=50, sparse_output=False),
    ),
    ("target", TargetEncoder(target_type="continuous")),
]

# %%
# Next, we create and fit the pipelines to the training data and evalute them on the
# test set:
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error

results = []
for name, categorical_preprocessor in categorical_preprocessors:
    preprocessor = ColumnTransformer(
        [
            ("nummerical", "passthrough", numerical_features),
            ("categorical", categorical_preprocessor, categorical_features),
        ]
    )
    pipe = make_pipeline(preprocessor, HistGradientBoostingRegressor(random_state=0))
    pipe.fit(X_train, y_train)

    rmse = mean_squared_error(y_test, pipe.predict(X_test), squared=False)
    results.append({"preprocessor": name, "root mean squared error": rmse})

# %%
# Finally, we display the results for all the encoders. When evaluting the
# predictive performance on the test set, dropping the categories perform the
# worst and the target encoder performs the best. The target encoding provides
# more information about the target, which the gradient boosted model at the end
# of the pipeline can take advantage of.
import pandas as pd

pd.DataFrame(results).sort_values("root mean squared error")
