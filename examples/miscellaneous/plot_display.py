"""
=================================================================
Displaying a Pipeline chaining a PCA and SVC
=================================================================

This example constructs a pipeline that does dimensionality reduction followed by a
support vector classifier.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then call the pipeline object.
"""

# %%
# Illustration of `Pipeline` and `PCA` and `SVC`
###############################################################################
# This section constructs a pipeline and displays its text representation.

from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn import set_config

steps = [('reduce_dim', PCA()), ('clf', SVC())]
pipe = Pipeline(steps)
pipe

# %%
# To visualize the diagram, change display='diagram'
set_config(display='diagram')

pipe
