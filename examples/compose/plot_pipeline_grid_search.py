"""
=================================================================
Displaying a Grid Search over a Pipeline
=================================================================

This example constructs a pipeline that does column transformation followed
by a random forest classifier over a grid search.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of `GridSearchCV` over a `Pipeline` with `RandomForest`
###############################################################################
# This section constructs a pipeline and displays its text and visual
# representation.

from sklearn.compose import make_column_selector as selector
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
import numpy as np
import pandas as pd
from sklearn import set_config

numerical_columns_selector = selector(dtype_exclude=object)
categorical_columns_selector = selector(dtype_include=object)

numeric_transformer = Pipeline(
    steps=[
        ("imputation_mean", SimpleImputer(missing_values=np.nan, strategy="mean")),
        ("scaler", StandardScaler()),
    ]
)

categorical_transformer = Pipeline(
    steps=[
        (
            "imputation_constant",
            SimpleImputer(fill_value="missing", strategy="constant"),
        ),
        ("onehot", OneHotEncoder(handle_unknown="ignore")),
    ]
)

data = pd.DataFrame()

numerical_columns = numerical_columns_selector(data)
categorical_columns = categorical_columns_selector(data)

preprocessor = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, numerical_columns),
        ("cat", categorical_transformer, categorical_columns),
    ]
)

pipe = Pipeline(
    steps=[("preprocessor", preprocessor), ("classifier", RandomForestClassifier())]
)

param_grid = {
    "classifier__n_estimators": [200, 500],
    "classifier__max_features": ["auto", "sqrt", "log2"],
    "classifier__max_depth": [4, 5, 6, 7, 8],
    "classifier__criterion": ["gini", "entropy"],
}

grid_search = GridSearchCV(pipe, param_grid=param_grid, n_jobs=1)

# %%
# To view the text pipeline, the default is `display='text'`
set_config(display="text")
pipe

# %%
# To visualize the diagram, change to `display='diagram'`
set_config(display="diagram")
pipe
