"""
=================================================================
Displaying a Basic Pipeline
=================================================================

This example constructs a basic pipeline that does logistic regression.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of a Basic `Pipeline` with `LogisticRegression`
###############################################################################
# This section constructs a pipeline and displays its text and visual
# representation.

from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn import set_config

steps = [("logistic_regression", LogisticRegression())]
pipe = Pipeline(steps)

# The default in Jupyter notebook is `display='text'`
set_config(display="text")
pipe

# %%
# To visualize the diagram, change to `display='diagram'`
set_config(display="diagram")
pipe

"""
========================================================================
Displaying a Pipeline Chaining a Preprocessing Step and Classifier
========================================================================

This example constructs a pipeline that does standard scaling followed by a
logistic regression classifier.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of `Pipeline` and `StandardScaler` and `LogisticRegression`
###############################################################################
# This section constructs a pipeline and displays its text representation.

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import set_config

steps = [
    ("standard_scaler", StandardScaler()),
    ("logistic_regression", LogisticRegression()),
]
pipe = Pipeline(steps)
pipe

# %%
# To visualize the diagram, change `display='diagram'`
set_config(display="diagram")

pipe

"""
========================================================================
Displaying a Pipeline Chaining Multiple Preprocessing Steps & Classifier
========================================================================

This example constructs a pipeline that does polynomial features and standard
scaling followed by a logistic regression classifier.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of `Pipeline` and `PolynomialFeatures`, `StandardScaler`
# and `LogisticRegression`
###############################################################################
# This section constructs a pipeline and displays its text representation.

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from sklearn import set_config

steps = [
    ("polynomial", PolynomialFeatures()),
    ("standard_scaler", StandardScaler()),
    ("logistic_regression", LogisticRegression()),
]
pipe = Pipeline(steps)
pipe

# To visualize the diagram, change to display='diagram'
set_config(display="diagram")
pipe

"""
=================================================================
Displaying a Pipeline Chaining a PCA and SVC
=================================================================

This example constructs a pipeline that does dimensionality reduction followed
by a support vector classifier.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of `Pipeline` and `PCA` and `SVC`
###############################################################################
# This section constructs a pipeline and displays its text representation.

from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn import set_config

steps = [("reduce_dim", PCA()), ("clf", SVC())]
pipe = Pipeline(steps)
pipe

# %%
# To visualize the diagram, change display='diagram'
set_config(display="diagram")

pipe

"""
=================================================================
Displaying a Complex Pipeline Chaining a Column Transformer
=================================================================

This example constructs a pipeline that does column transformation followed
by a logistic regression classifier.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of `Pipeline` and `ColumnTransformer` and `LogisticRegression`
###############################################################################
# This section constructs a pipeline and displays its text representation.

from sklearn.compose import make_column_selector as selector
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
import pandas as pd
from sklearn import set_config

numerical_columns_selector = selector(dtype_exclude=object)
categorical_columns_selector = selector(dtype_include=object)

data = pd.DataFrame()

numerical_columns = numerical_columns_selector(data)
categorical_columns = categorical_columns_selector(data)

categorical_preprocessor = OneHotEncoder(handle_unknown="ignore")
numerical_preprocessor = StandardScaler()

preprocessor = ColumnTransformer(
    [
        ("one-hot-encoder", categorical_preprocessor, categorical_columns),
        ("standard-scaler", numerical_preprocessor, numerical_columns),
    ]
)

pipe = make_pipeline(preprocessor, LogisticRegression(max_iter=500))
pipe

# %%
# To visualize the diagram, change display='diagram'
set_config(display="diagram")

pipe

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
# This section constructs a pipeline and displays its text representation.

from sklearn.compose import make_column_selector as selector
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline
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

pipe

# %%
# To visualize the diagram, change display='diagram'
set_config(display="diagram")

pipe
