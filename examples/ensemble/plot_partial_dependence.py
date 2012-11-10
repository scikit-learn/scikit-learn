"""
========================
Partial Dependence Plots
========================

Partial dependence plots show the dependence between the target function [1]_
and a set of 'target' features, marginalizing over the
values of all other features (the complement features). Due to the limits
of human perception the size of the target feature set must be small (usually,
one or two) thus the target features are usually chosen among the most
important features (see :meth:`ensemble.GradientBoostingRegressor.feature_importances_`).

This example shows how to obtain partial dependence plots from a
:class:`ensemble.GradientBoostingRegressor` trained on the California
housing dataset. The example is taken from [HTF2009]_.

The first plot shows four one-dimensional partial dependence plots. The target
variables are: median income (`MedInc`), avg. occupants per household (`AvgOccup`),
median house age (`HouseAge`), and avg. rooms per household (`AveRooms`).
We can clearly see that the median house price shows a linear relationship
with the median income (top left) and that the house price drops when the
avg. occupants per household increases (top right).
The lower left plot shows that the house age in a district does not have
a strong influence on the (median) house price.
The red tick marks on the x-axis represent the deciles of the feature values
in the trainig data.

Partial dependence plots with two target features enable us to visualize
interactions among them. The two-variable partial dependence plot shows the
dependence of median house price on joint values of house age and avg.
occupants per household. We can clearly see an interaction between the
two features:
For an avg. occupancy greather than two, the house price is nearly independent
of the house age, whereas for values less than two there is a strong dependence
on age.

.. [HTF2009] T. Hastie, R. Tibshirani and J. Friedman, "Elements of Statistical Learning Ed. 2", Springer, 2009.

.. [1] For classification you can think of it as the regression score before the link function.
"""
print __doc__

import numpy as np
import pylab as pl
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
from scipy.stats.mstats import mquantiles

from sklearn.cross_validation import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import gradient_boosting
from sklearn.datasets.california_housing import fetch_california_housing

# fetch California housing dataset
cal_housing = fetch_california_housing()

# split 80/20 train-test
X_train, X_test, y_train, y_test = train_test_split(cal_housing.data,
                                                    cal_housing.target,
                                                    test_size=0.2,
                                                    random_state=1)
names = cal_housing.feature_names

print('_' * 80)
print("Training GBRT...")
clf = GradientBoostingRegressor(n_estimators=800, max_depth=4,
                                learning_rate=0.1, loss='huber',
                                random_state=1)
clf.fit(X_train, y_train)
print("done.")

print('_' * 80)
print('One-way partial dependence plots')
print
fig = pl.figure()
sub_plots = []

for i, fx in enumerate([0, 5, 1, 2]):
    name = names[fx]
    target_feature = np.array([fx], dtype=np.int32)

    ax = pl.subplot(2, 2, i + 1)

    # plot partial dependence
    pdp, (axis,) = gradient_boosting.partial_dependence(clf, target_feature,
                                                        X=X_train)
    ax.plot(axis, pdp.ravel(), 'g-')

    # plot data deciles
    deciles = mquantiles(X_train[:, fx], prob=np.arange(0.1, 1.0, 0.1))
    trans = matplotlib.transforms.blended_transform_factory(ax.transData,
                                                            ax.transAxes)
    ax.vlines(deciles, 0.0, 0.05, transform=trans, color='red')

    pl.xlabel(name)
    pl.ylabel('Partial Dependence')

    sub_plots.append(ax)

# set common ylim
y_min = min((ax.get_ylim()[0] for ax in sub_plots))
y_max = max((ax.get_ylim()[1] for ax in sub_plots))
for ax in sub_plots:
    ax.set_ylim((y_min, y_max))

fig.suptitle('Partial dependence of house value on nonlocation features ' \
             'for the California housing dataset')

print('_' * 80)
print('Two-way partial dependence plot')
print
fig = pl.figure()


target_feature = np.array([1, 5], dtype=np.int32)
pdp, (x_axis, y_axis) = gradient_boosting.partial_dependence(clf,
                                                             target_feature,
                                                             X=X_train,
                                                             grid_resolution=50)
XX, YY = np.meshgrid(x_axis, y_axis)
Z = pdp.T.reshape(XX.shape).T
ax = Axes3D(fig)
surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=pl.cm.BuPu)
ax.set_xlabel(names[target_feature[0]])
ax.set_ylabel(names[target_feature[1]])
ax.set_zlabel('Partial dependence')
#  pretty init view
ax.view_init(elev=22, azim=122)
pl.colorbar(surf)
pl.suptitle('Partial dependence of house value on median age and ' \
            'average occupancy')
pl.show()
