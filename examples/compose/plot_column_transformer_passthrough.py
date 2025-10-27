# examples/compose/plot_column_transformer_passthrough.py
# -*- coding: utf-8 -*-

"""
=============================================================
ColumnTransformer with remainder='passthrough' and Pandas
=============================================================

This example shows how to keep columns untouched with
``remainder='passthrough'`` while transforming others.
The input is a pandas DataFrame – the most common real-world case.
"""

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report

# -------------------------------------------------
# 1. Create a realistic Pandas DataFrame
# -------------------------------------------------
data = pd.DataFrame({
    "age":      [25, 30, 35, 40, 45, 50, 22, 33],
    "salary":   [50000, 60000, 70000, 80000, 90000, 100000, 45000, 55000],
    "city":     ["NY", "LA", "NY", "SF", "LA", "NY", "SF", "LA"],
    "is_senior":[0, 0, 0, 1, 1, 1, 0, 0],
    "target":   [0, 1, 0, 1, 1, 0, 0, 1]
})

X = data.drop("target", axis=1)
y = data["target"]

# -------------------------------------------------
# 2. ColumnTransformer – scale numeric, encode city,
#     passthrough the binary column `is_senior`
# -------------------------------------------------
ct = ColumnTransformer(
    [
        ("scale", StandardScaler(), ["age", "salary"]),
        ("encode", OneHotEncoder(drop="first", sparse_output=False), ["city"]),
    ],
    remainder="passthrough",   # <-- keeps `is_senior` unchanged
)

# -------------------------------------------------
# 3. Full pipeline + LogisticRegression
# -------------------------------------------------
pipe = Pipeline(
    [
        ("transform", ct),
        ("clf", LogisticRegression(max_iter=1000)),
    ]
)

# -------------------------------------------------
# 4. Train / test split & evaluation
# -------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)

pipe.fit(X_train, y_train)
y_pred = pipe.predict(X_test)

print(f"Accuracy: {accuracy_score(y_test, y_pred):.3f}")
print(classification_report(y_test, y_pred))

# -------------------------------------------------
# 5. Visualise the transformed feature matrix
# -------------------------------------------------
transformed = ct.fit_transform(X)
cols = (
    ["age_scaled", "salary_scaled"] +
    [f"city_{c}" for c in ct.named_transformers_["encode"].get_feature_names_out()] +
    ["is_senior"]
)
print("\nTransformed features (first 5 rows):")
print(pd.DataFrame(transformed, columns=cols).head())
