"""
ColumnTransformer for Preprocessing
====================================

Demonstrates how to use :class:`~sklearn.compose.ColumnTransformer` to apply
different preprocessing to numeric and categorical features in a pipeline.

This example shows:
- Scaling numeric features with :class:`~sklearn.preprocessing.StandardScaler`
- Encoding categorical features with :class:`~sklearn.preprocessing.OneHotEncoder`
- Combining preprocessing with a :class:`~sklearn.linear_model.LogisticRegression` classifier

For more details, see the :ref:`Column Transformer <column_transformer>` user guide.
"""

# %%
# Imports and Setup
# -----------------

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.datasets import fetch_openml

set_config(display="diagram")

# %%
# Load the Titanic Dataset
# -------------------------
# We'll use the Titanic dataset from OpenML, which contains both numeric
# and categorical features.

X, y = fetch_openml("titanic", version=1, as_frame=True, return_X_y=True, parser="pandas")

# Select relevant columns
num_cols = ["age", "fare"]
cat_cols = ["sex", "pclass", "embarked"]

X = X[num_cols + cat_cols]
print(f"Dataset shape: {X.shape}")
print(f"\nFirst few rows:\n{X.head()}")

# %%
# Define the Preprocessing Pipeline
# ----------------------------------
# Use ColumnTransformer to apply different transformations to different columns.

preprocessor = ColumnTransformer(
    transformers=[
        ("num", StandardScaler(), num_cols),
        ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), cat_cols),
    ]
)

# %%
# Create the Full Pipeline
# -------------------------
# Combine preprocessing with a classifier.

clf = make_pipeline(
    preprocessor,
    LogisticRegression(max_iter=1000, random_state=42)
)

# Display the pipeline structure
clf

# %%
# Train and Evaluate the Model
# -----------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# Fit the pipeline
clf.fit(X_train, y_train)

# Evaluate on test set
train_score = clf.score(X_train, y_train)
test_score = clf.score(X_test, y_test)

print(f"Training accuracy: {train_score:.3f}")
print(f"Test accuracy: {test_score:.3f}")

# %%
# Inspect the Transformed Features
# ---------------------------------
# View the feature names after transformation.

feature_names = (
    num_cols +
    list(preprocessor.named_transformers_["cat"].get_feature_names_out(cat_cols))
)

print(f"\nTotal features after transformation: {len(feature_names)}")
print(f"Feature names: {feature_names}")
