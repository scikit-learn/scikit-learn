"""
==============================================
Individual Conditional Expectation (ICE) Plots
==============================================

Similar to a partial dependence (PD) plot, an individual conditional
expectation (ICE) plot [1]_ shows the dependence between the target function
[2]_ and a 'target' feature, marginalizing over the values of all other
features (the complement features). However, unlike PD plots which show the
average effect of the 'target' features, ICE plots visualizes the dependence
of the prediction on a feature for each instance separately, with one line
per instance.

This example shows how to obtain ICE plots from a
:class:`~sklearn.neural_network.MLPRegressor` trained on the
California housing dataset.

The target variables for ICE plot are: median income (`MedInc`), average
occupants per household (`AvgOccup`), median house age (`HouseAge`), and
average rooms per household (`AveRooms`).

.. [1] Goldstein, A., Kapelner, A., Bleich, J., and Pitkin, E., Peeking Inside
       the Black Box: Visualizing Statistical Learning With Plots of
       Individual Conditional Expectation. (2015) Journal of Computational and
       Graphical Statistics, 24(1): 44-65 (https://arxiv.org/abs/1309.6392)

.. [2] For classification you can think of it as the regression score before
       the link function.
"""
print(__doc__)

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import QuantileTransformer
from sklearn.pipeline import make_pipeline

from sklearn.inspection import plot_individual_conditional_expectation
from sklearn.inspection import plot_partial_dependence
from sklearn.neural_network import MLPRegressor
from sklearn.datasets import fetch_california_housing


##############################################################################
# Loading California Housing data
# -------------------------------

cal_housing = fetch_california_housing()
X = pd.DataFrame(cal_housing.data, columns=cal_housing.feature_names)
y = cal_housing.target

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,
                                                    random_state=0)

##############################################################################
# Configuring a training pipeline
# -------------------------------
#
# Let's configuring a training pipeline fit a MLPRegressor:

print("Training MLPRegressor...")
est = make_pipeline(QuantileTransformer(),
                    MLPRegressor(hidden_layer_sizes=(50, 50),
                                 learning_rate_init=0.01,
                                 early_stopping=True))
est.fit(X_train, y_train)
print("Test R2 score: {:.2f}".format(est.score(X_test, y_test)))


# We configured a pipeline to scale the numerical input features and tuned the
# neural network size and learning rate to get a reasonable compromise between
# training time and predictive performance on a test set.
#
# Note that it is important to check that the model is accurate enough on a
# test set before plotting the ICE plots since there would be little use in
# explaining the impact of a given feature on the prediction function of a
# poor model.


##############################################################################
# ICE computation
# ---------------
#
# Let's now compute the ICE plots for the neural network:
#
# One of the main limitation of ICE plots is that the plot get overcrowded
# if many ICE curves are drawn. Due to this we will use a random sample of
# 50 instances for the ICE plot.

X_train_sample = X_train.sample(50, random_state=0)

features = ['MedInc', 'AveOccup', 'HouseAge', 'AveRooms']

print('Computing ICE plots...')
plot_individual_conditional_expectation(est, X_train_sample, features,
                                        n_jobs=3, grid_resolution=20,
                                        centre=False, n_cols=2,
                                        line_kw={'linewidth': 0.5})
fig = plt.gcf()
fig.suptitle('ICE of house value on non-location features')
fig.subplots_adjust(hspace=0.3)


# As the ICE curves do not have the same origin, sometimes it can be hard to
# see the differences between ICE curves. The ICE curves can be made to have
# same origin with ``centre`` parameter set to True.


##############################################################################
# Centered ICE computation
# ------------------------
#
# Let's now compute the ICE plots with ``centre`` parameter set to True:

print('Computing centered ICE plots...')
plot_individual_conditional_expectation(est, X_train_sample, features,
                                        n_jobs=3, grid_resolution=20,
                                        centre=True, n_cols=2,
                                        line_kw={'linewidth': 0.5})
fig = plt.gcf()
fig.suptitle('Centered ICE of house value on non-location features')
fig.subplots_adjust(hspace=0.3)


##############################################################################
# Partial Dependence computation
# ------------------------------
#
# In ICE plots it might not be easy to see the average effect of the 'target'
# variable. Hence, it is recommended to use ICE plot with partial dependency
# plots.
#
# Let's compute single-variable partial dependence plots

print('Computing partial dependence plots...')
plot_partial_dependence(est, X_train, features, n_jobs=3, grid_resolution=20,
                        n_cols=2)
fig = plt.gcf()
fig.suptitle('Partial dependence of house value on non-location features')
fig.subplots_adjust(hspace=0.3)

plt.show()

##############################################################################
# Analysis of the plots
# ---------------------
#
# From the PD plot, we can see that the median house price increases with the
# median income (top left) and that the median house price drops when the
# average occupants per household increases (top right). However, from ICE
# plots we can see that there are some exceptions, where the house price
# remain constant with median income and average occupants.
# On the other hand, while the house age (bottom left) does not have a strong
# influence on the median house price on average, there seems to be a number
# of exceptions  where the house price increase when between the ages 15-25.
# Similar exceptions can be observed for average number of rooms (bottom
# right).
