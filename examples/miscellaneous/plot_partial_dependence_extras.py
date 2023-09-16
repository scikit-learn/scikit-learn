"""
=========================================================
Extra Feature Visualizations for Partial Dependence Plots
=========================================================

Partial dependence plots can be misleading because partial dependence values
are not calculated on actual data points. They are calculated on a grid
representing a range of the feature values. This means that the partial
dependence plots can be misinterpreted in regions where the data are sparse.
[1]_ This problem can be avoided by plotting a scatter, boxplot, or histogram
(or bar chart if the feature is categorical) of the actual feature values for
each point in the partial dependence plot.

This example shows how to add extra visualizations to partial dependence plots
to better understand the relationship between the target and the features. The
focus of this example is to illustrate how extra plots can aid in
interpeting partial dependence results.

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
# Set ``extra_plots='hist'`` to create figure above the partial dependence
# plot.  Note that we're plotting the test partition.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf, X_test, ["bp"], extra_plots="hist", ax=ax
)
# %%
# Now we will add a boxplot to the one-way blood pressure plot as well as
# histograms to the blood pressure and bmi two-way plot.
# To achieve this, pass a list containing 'boxplot', 'hist', or 'scatter' to
#  ``extra_plots``. You can also pass ``None```, for the default behavior of no
# extra plots.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf, X_test, ["bp", ("bmi", "bp")], extra_plots=["boxplot", "hist"], ax=ax
)

# %%
# Customize Axes
# --------------
# The axes of the extra plots can be customized by passing a dictionary of
# dictionaries to ``extra_plots_kw``. The keys of the dictionary are the names
# of the extra plots and the values are dictionaries of keyword arguments to
# pass to the plotting function. In the example below, we remove the fill from
# the histogram and change the color of the boxplot.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["bp", ("bmi", "bp")],
    extra_plots=["boxplot", "hist"],
    extra_plots_kw={"hist": {"fill": False}, "boxplot": {"widths": 0.2}},
    ax=ax,
)


# %%
# Scatter Plots
# -------------
# Scatter plots can be used to visualize the relationship between the target
# and the feature values or between the target and two feature values. The
# ``y`` parameter controls what is plotted on the y-axis.
#
# In one-way partial dependence plots, if ``y`` is ``None``, then predictions
# (regression) or predicted probabilities (classification) are
# plotted. You can pass ``y`` values as an argument to plot them against the
# feature values.
#
# In two-way partial dependence plots, ``y`` must be provided as an argument
# to visualize the points.
#
# In the example below, we plot the true target values for one and two-way
# partial dependece plots. We also change the alpha of the scatter to be
# darker. Note, that we only pass 'scatter' as a string and it applies to all
# plots.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["s5", ("bmi", "bp")],
    extra_plots="scatter",
    extra_plots_kw={"scatter": {"alpha": 0.25}},
    y=y_test,
    ax=ax,
)

# %%
# Scatter Plots with predictions
# ------------------------------
# We can also plot the predictions (regression) or predicted probabilities
# (classification) on the y-axis. In the example below, we plot the predictions
# as a scatter behind the partial dependence line. Compare this figure to the
# one above where we plotted the true target values against s5. Note that the
# legend changed from "True Values" to "Predictions".
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["s5", ("bmi", "bp")],
    extra_plots="scatter",
    extra_plots_kw={"scatter": {"alpha": 0.25}},
    ax=ax,
)

# %%
# Categoricals
# ------------
# There is support for visualizing the feature distributions of categorical
# features in partial dependence plots by setting ``extra_plots='hist'``,
# and when ``categorical_features`` is not `None`. A ``bar`` plot is used
# instead of a ``hist`` plot when the feature is a categorical feature.
_, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
display = PartialDependenceDisplay.from_estimator(
    rf,
    X_test,
    ["sex"],
    extra_plots="hist",
    categorical_features=["sex"],
    ax=ax,
)
