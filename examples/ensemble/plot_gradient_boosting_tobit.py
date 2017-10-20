"""
=============================================
Grabit: Gradient Boosting for the Tobit model
=============================================

Demonstrate the Grabit model (Gradient Boosting with a Tobit loss fuction).

This example fits a Gradient Boosting model with a Tobit loss function
to a simulated data set.

Feature importance plots and partial dependence plots are shown.
"""

# Author: Fabio Sigrist <fabio.sigrist@gmail.com>
#
# License: BSD 3 clause


print(__doc__)

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble.partial_dependence import plot_partial_dependence
from sklearn.ensemble.partial_dependence import partial_dependence
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


"""
Simulate data
"""
n = 10000
X = np.random.rand(n, 4)
X = (X - 0.5) * 2
x1t = np.abs(X[:, 0]) ** 1.5
x2t = np.abs(X[:, 1]) ** 1.5
r2 = (x1t ** 2 + x2t ** 2) ** 0.5
e = np.random.normal(scale=1, size=n)
y = (4 * np.cos(np.pi * 2 * 1.1 * r2) + 4 * X[:, 2] + e)

# Censoring: 66% of the data is censored (33% lower and 33% upper censoring)
yc = y.copy()
yl = np.percentile(y, q=33)
yu = np.percentile(y, q=66)
yc[y >= yu] = yu
yc[y <= yl] = yl

"""
Define and learn model
"""
model = GradientBoostingRegressor(loss='tobit', yl=yl, yu=yu)
model.fit(X, yc)


"""
Variable importance plots
"""
feature_names = np.array(['V1', 'V2', 'V3', 'V4'])
feature_importance = model.feature_importances_
# make importances relative to max importance
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
plt.subplot(1, 2, 2)
plt.barh(pos, feature_importance[sorted_idx], align='center')
plt.yticks(pos, feature_names[sorted_idx])
plt.xlabel('Relative Importance')
plt.title('Variable Importance')
plt.show()


"""
Partial dependence plots
"""
# Univariate partial dependence plots
features = [0, 1, 2, 3]
names = [0, 1, 2, 3]
fig, axs = plot_partial_dependence(model, X, features,
                                   feature_names=names,
                                   n_jobs=3, grid_resolution=50)
fig.suptitle('Partial dependence plots')
plt.subplots_adjust(top=0.9)

# 3d plot for bivariate partial dependence plot
target_feature = (0, 1)
pdp, axes = partial_dependence(model, target_feature,
                               X=X, grid_resolution=50)
x1, x2 = np.meshgrid(axes[0], axes[1])
yh = pdp[0].reshape(list(map(np.size, axes))).T

fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(x1, x2, yh, rstride=1, cstride=1,
                       cmap=plt.cm.BuPu, edgecolor='k')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('y')
#  pretty init view
ax.view_init(elev=45, azim=120)
plt.colorbar(surf)
plt.suptitle('Partial dependence plot of first two variables')
plt.subplots_adjust(top=0.9)
plt.show()
