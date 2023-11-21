"""
==========================================================================
Marginal Distribution Plots of Feature Values in Partial Dependence Plots
==========================================================================

Partial dependence plots can be misleading because partial dependence values
are not calculated on actual data points. They are calculated on a grid
representing a range of the feature values. This means that the partial
dependence plots can be misinterpreted in regions where the data are sparse.
[1]_ This problem can be avoided by plotting a histogram
(or bar chart if the feature is categorical) of the actual feature values for
each point in the partial dependence plot.

This example shows how to add visualizations of marginal distributions to
partial dependence plots. The focus of this example is to illustrate how these
extra visualizations can aid with the interpetation of partial dependence
results. Notice that we're setting ``percentiles=(0, 1)`` in these examples
so we visualize the entire range of the feature values.

For more information on partial dependence plots,
see the :ref:`User guide <partial_dependence>`.

.. [1] `Molnar, Christoph. "Interpretable machine learning.
       A Guide for Making Black Box Models Explainable",
       2019. <https://christophm.github.io/interpretable-ml-book/>`_

"""

# %%
# Diabetes dataset preprocessing
# ----------------------------------
#
# We will use diabetes dataset which is a regression task. The goal is to
# predict disease progression based on 10 features. The following cell loads
# the dataset, splits it into training and test sets, and one-hot encodes "sex"
# as a categorical variable.
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.datasets import load_diabetes
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import PartialDependenceDisplay
from sklearn.model_selection import train_test_split

X, y = load_diabetes(as_frame=True, return_X_y=True)
X["sex"] = pd.get_dummies(X["sex"], drop_first=True)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, random_state=0, train_size=0.6
)

# %%
# Modeling
# --------
# Now fit a random forest model to the training data. We're not going to
# tune the model in this tutorial so we can focus on visualizing the features
# alongside the partial dependence plots.
rf = RandomForestRegressor(max_depth=5, random_state=0)
rf.fit(X_train, y_train)

# %%
# Histograms
# ----------
# The first additional plot we will create is a histogram of the feature values
# for blood pressure. This will help us understand the distribution of the
# feature values in the input data.
#
# Set ``marginal_dist=True`` to create a histogram above the partial dependence
# plot.  Note that we're plotting the test partition.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["bp"],
    marginal_dist=True,
    percentiles=(0, 1),
    ax=ax,
)
# %%
# You can control which partial dependence plots get a marginal distribution plot
# by passing a list of booleans to ``marginal_dist``. In the example below,
# we will create a histogram for the two-way partial dependence plot of bmi and
# blood pressure, but not for the one-way partial dependence plot of blood
# pressure.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["bp", ("bmi", "bp")],
    marginal_dist=[False, True],
    percentiles=(0, 1),
    ax=ax,
)

# %%
# Categoricals
# ------------
# There is support for visualizing the marginal distributions of categorical
# features in partial dependence plots by setting ``marginal_dist=True``,
# when ``categorical_features`` is not `None`. A ``bar`` plot is used
# instead of a ``hist`` plot when the feature is a categorical feature.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["sex"],
    marginal_dist=True,
    percentiles=(0, 1),
    categorical_features=["sex"],
    ax=ax,
)

# %%
# Customize Axes
# --------------
# The axes of the marginal distribution plots can be customized by passing a
# dictionary of dictionaries to ``marginal_dist_kw``. The keys of the
# dictionary are the names of the plot types (either 'hist' or 'bar') and the
# values are dictionaries of keyword arguments to pass to the plotting
# function. In the example below, we change decrease transparency of the bars
# in the bar chart and increase the number of bins in the histogram.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["sex", ("bmi", "bp")],
    categorical_features=["sex"],
    marginal_dist=[True, True],
    marginal_dist_kw={"bar": {"alpha": 1}, "hist": {"bins": 20}},
    percentiles=(0, 1),
    ax=ax,
)
