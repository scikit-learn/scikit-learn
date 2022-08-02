"""
===============================================================
Partial Dependence and Individual Conditional Expectation Plots
===============================================================

Partial dependence plots show the dependence between the target function [2]_
and a set of features of interest, marginalizing over the values of all other
features (the complement features). Due to the limits of human perception, the
size of the set of features of interest must be small (usually, one or two)
thus they are usually chosen among the most important features.

Similarly, an individual conditional expectation (ICE) plot [3]_
shows the dependence between the target function and a feature of interest.
However, unlike partial dependence plots, which show the average effect of the
features of interest, ICE plots visualize the dependence of the prediction on a
feature for each :term:`sample` separately, with one line per sample.
Only one feature of interest is supported for ICE plots.

This example shows how to obtain partial dependence and ICE plots from a
:class:`~sklearn.neural_network.MLPRegressor` and a
:class:`~sklearn.ensemble.HistGradientBoostingRegressor` trained on the
California housing dataset. The example is taken from [1]_.

.. [1] T. Hastie, R. Tibshirani and J. Friedman, "Elements of Statistical
       Learning Ed. 2", Springer, 2009.

.. [2] For classification you can think of it as the regression score before
       the link function.

.. [3] :arxiv:`Goldstein, A., Kapelner, A., Bleich, J., and Pitkin, E. (2015).
       "Peeking Inside the Black Box: Visualizing Statistical Learning With Plots of
       Individual Conditional Expectation". Journal of Computational and
       Graphical Statistics, 24(1): 44-65 <1309.6392>`

"""

# %%
# California Housing data preprocessing
# -------------------------------------
#
# Center target to avoid gradient boosting init bias: gradient boosting
# with the 'recursion' method does not account for the initial estimator
# (here the average target, by default).

import pandas as pd
from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split

cal_housing = fetch_california_housing()
X = pd.DataFrame(cal_housing.data, columns=cal_housing.feature_names)
y = cal_housing.target

y -= y.mean()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)

# %%
# 1-way partial dependence with different models
# ----------------------------------------------
#
# In this section, we will compute 1-way partial dependence with two different
# machine-learning models: (i) a multi-layer perceptron and (ii) a
# gradient-boosting. With these two models, we illustrate how to compute and
# interpret both partial dependence plot (PDP) and individual conditional
# expectation (ICE).
#
# Multi-layer perceptron
# ......................
#
# Let's fit a :class:`~sklearn.neural_network.MLPRegressor` and compute
# single-variable partial dependence plots.

from time import time
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import QuantileTransformer
from sklearn.neural_network import MLPRegressor

print("Training MLPRegressor...")
tic = time()
est = make_pipeline(
    QuantileTransformer(),
    MLPRegressor(
        hidden_layer_sizes=(30, 15),
        learning_rate_init=0.01,
        early_stopping=True,
        random_state=0,
    ),
)
est.fit(X_train, y_train)
print(f"done in {time() - tic:.3f}s")
print(f"Test R2 score: {est.score(X_test, y_test):.2f}")

# %%
# We configured a pipeline to scale the numerical input features and tuned the
# neural network size and learning rate to get a reasonable compromise between
# training time and predictive performance on a test set.
#
# Importantly, this tabular dataset has very different dynamic ranges for its
# features. Neural networks tend to be very sensitive to features with varying
# scales and forgetting to preprocess the numeric feature would lead to a very
# poor model.
#
# It would be possible to get even higher predictive performance with a larger
# neural network but the training would also be significantly more expensive.
#
# Note that it is important to check that the model is accurate enough on a
# test set before plotting the partial dependence since there would be little
# use in explaining the impact of a given feature on the prediction function of
# a poor model.
#
# We will plot the partial dependence, both individual (ICE) and averaged one
# (PDP). We limit to only 50 ICE curves to not overcrowd the plot.

from sklearn.inspection import PartialDependenceDisplay

common_params = {
    "subsample": 50,
    "n_jobs": 2,
    "grid_resolution": 20,
    "centered": True,
    "random_state": 0,
}

print("Computing partial dependence plots...")
tic = time()
display = PartialDependenceDisplay.from_estimator(
    est,
    X_train,
    features=["MedInc", "AveOccup", "HouseAge", "AveRooms"],
    kind="both",
    **common_params,
)
print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "Partial dependence of house value on non-location features\n"
    "for the California housing dataset, with MLPRegressor"
)
display.figure_.subplots_adjust(hspace=0.3)

# %%
# Gradient boosting
# .................
#
# Let's now fit a :class:`~sklearn.ensemble.HistGradientBoostingRegressor` and
# compute the partial dependence on the same features.

from sklearn.ensemble import HistGradientBoostingRegressor

print("Training HistGradientBoostingRegressor...")
tic = time()
est = HistGradientBoostingRegressor(random_state=0)
est.fit(X_train, y_train)
print(f"done in {time() - tic:.3f}s")
print(f"Test R2 score: {est.score(X_test, y_test):.2f}")

# %%
# Here, we used the default hyperparameters for the gradient boosting model
# without any preprocessing as tree-based models are naturally robust to
# monotonic transformations of numerical features.
#
# Note that on this tabular dataset, Gradient Boosting Machines are both
# significantly faster to train and more accurate than neural networks. It is
# also significantly cheaper to tune their hyperparameters (the defaults tend
# to work well while this is not often the case for neural networks).
#
# We will plot the partial dependence, both individual (ICE) and averaged one
# (PDP). We limit to only 50 ICE curves to not overcrowd the plot.

print("Computing partial dependence plots...")
tic = time()
display = PartialDependenceDisplay.from_estimator(
    est,
    X_train,
    features=["MedInc", "AveOccup", "HouseAge", "AveRooms"],
    kind="both",
    **common_params,
)
print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "Partial dependence of house value on non-location features\n"
    "for the California housing dataset, with Gradient Boosting"
)
display.figure_.subplots_adjust(wspace=0.4, hspace=0.3)

# %%
# Analysis of the plots
# .....................
#
# We can clearly see on the PDPs (dashed orange line) that the median house price
# shows a linear relationship with the median income (top left) and that the
# house price drops when the average occupants per household increases (top
# middle). The top right plot shows that the house age in a district does not
# have a strong influence on the (median) house price; so does the average
# rooms per household.
#
# The ICE curves (light blue lines) complement the analysis: we can see that
# there are some exceptions (which are better highlighted with the option
# `centered=True`), where the house price remains constant with respect to
# median income and average occupants variations.
# On the other hand, while the house age (top right) does not have a strong
# influence on the median house price on average, there seems to be a number
# of exceptions where the house price increases when
# between the ages 15-25. Similar exceptions can be observed for the average
# number of rooms (bottom left). Therefore, ICE plots show some individual
# effect which are attenuated by taking the averages.
#
# In all plots, the tick marks on the x-axis represent the deciles of the
# feature values in the training data.
#
# We also observe that :class:`~sklearn.neural_network.MLPRegressor` has much
# smoother predictions than
# :class:`~sklearn.ensemble.HistGradientBoostingRegressor`.
#
# However, it is worth noting that we are creating potential meaningless
# synthetic samples if features are correlated.

# %%
# 2D interaction plots
# --------------------
#
# PDPs with two features of interest enable us to visualize interactions among
# them. However, ICEs cannot be plotted in an easy manner and thus interpreted.
# Another consideration is linked to the performance to compute the PDPs. With
# the tree-based algorithm, when only PDPs are requested, they can be computed
# on an efficient way using the `'recursion'` method.
import matplotlib.pyplot as plt

print("Computing partial dependence plots...")
tic = time()
_, ax = plt.subplots(ncols=3, figsize=(9, 4))

# Note that we could have called the method `from_estimator` three times and
# provide one feature, one kind of plot, and one axis for each call.
display = PartialDependenceDisplay.from_estimator(
    est,
    X_train,
    features=["AveOccup", "HouseAge", ("AveOccup", "HouseAge")],
    kind=["both", "both", "average"],
    ax=ax,
    **common_params,
)

print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "Partial dependence of house value on non-location features\n"
    "for the California housing dataset, with Gradient Boosting"
)
display.figure_.subplots_adjust(wspace=0.4, hspace=0.3)

# %%
# The two-way partial dependence plot shows the dependence of median house
# price on joint values of house age and average occupants per household. We
# can clearly see an interaction between the two features: for an average
# occupancy greater than two, the house price is nearly independent of the
# house age, whereas for values less than two there is a strong dependence on
# age.
#
# 3D interaction plots
# --------------------
#
# Let's make the same partial dependence plot for the 2 features interaction,
# this time in 3 dimensions.
import numpy as np

# unused but required import for doing 3d projections with matplotlib < 3.2
import mpl_toolkits.mplot3d  # noqa: F401

from sklearn.inspection import partial_dependence

fig = plt.figure()

features = ("AveOccup", "HouseAge")
pdp = partial_dependence(
    est, X_train, features=features, kind="average", grid_resolution=10
)
XX, YY = np.meshgrid(pdp["values"][0], pdp["values"][1])
Z = pdp.average[0].T
ax = fig.add_subplot(projection="3d")
fig.add_axes(ax)

surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor="k")
ax.set_xlabel(features[0])
ax.set_ylabel(features[1])
ax.set_zlabel("Partial dependence")
# pretty init view
ax.view_init(elev=22, azim=122)
plt.colorbar(surf)
plt.suptitle(
    "Partial dependence of house value on median\n"
    "age and average occupancy, with Gradient Boosting"
)
plt.subplots_adjust(top=0.9)
plt.show()
