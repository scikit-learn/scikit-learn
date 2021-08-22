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

# %%
# To view the text pipeline, the default is `display='text'`
set_config(display="text")
pipe

# %%
# To visualize the diagram, change to `display='diagram'`
set_config(display="diagram")
pipe
