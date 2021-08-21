"""
=================================================================
Displaying a Basic Pipeline with Logistic Regression
=================================================================

This example constructs a basic pipeline that does logistic regression.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then output the pipeline object.
"""

# %%
# Illustration of a Basic `Pipeline` with `LogisticRegression`
###############################################################################
# This section constructs a pipeline and displays its text representation.

from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn import set_config

steps = [("logistic_regression", LogisticRegression())]
pipe = Pipeline(steps)
pipe

# %%
# To visualize the diagram, change display='diagram'
set_config(display="diagram")

pipe

"""
========================================================================
Displaying a Pipeline Chaining a Standard Scaler and Logistic Regression
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
# To visualize the diagram, change display='diagram'
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
