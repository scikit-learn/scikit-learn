#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=================================================================
Pipeline: chaining a SVC and PCA
=================================================================

This example constructs a pipeline that does support vector
classifier followed by dimensionality reduction.

The default configuration for displaying a pipeline is 'text' where
`set_config(display='text')`.  To visualize the diagram in Jupyter Notebook,
use `set_config(display='diagram')` and then call the pipeline object.

# %%
Illustration of ``Pipeline`` and ``SVC`` and ``PCA``
###############################################################################

This section illustrates the use of a ``Pipeline`` with ``SVC`` and ``PCA``
"""

from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn import set_config

estimators = [('reduce_dim', PCA()), ('clf', SVC())]
pipe = Pipeline(estimators)
pipe

# to visualize the diagram, change display='diagram'
set_config(display='diagram')

pipe
