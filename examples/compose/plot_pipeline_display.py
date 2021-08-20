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
