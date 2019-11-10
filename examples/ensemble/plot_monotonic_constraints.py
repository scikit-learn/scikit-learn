"""
Example reproducing
https://xgboost.readthedocs.io/en/latest//tutorials/monotonic.html

Mostly for my beloved reviewers.
"""

from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import plot_partial_dependence
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.RandomState(0)

n_samples = 5000
f_0 = rng.random(size=n_samples)  # positive correlation with y
f_1 = rng.random(size=n_samples)  # negative correslation with y
X = np.c_[f_0, f_1]
noise = rng.normal(loc=0.0, scale=0.01, size=n_samples)
y = (5 * f_0 + np.sin(10 * np.pi * f_0) -
     5 * f_1 - np.cos(10 * np.pi * f_1) +
     noise)

fig, (ax0, ax1) = plt.subplots(ncols=2)


gbdt = HistGradientBoostingRegressor()
gbdt.fit(X, y)
unconstrained_disp0 = plot_partial_dependence(
    gbdt, X, features=[0], line_kw={'linewidth': 4, 'label': 'unconstrained'},
    ax=ax0)
unconstrained_disp1 = plot_partial_dependence(
    gbdt, X, features=[1], line_kw={'linewidth': 4, 'label': 'unconstrained'},
    ax=ax1)

gbdt = HistGradientBoostingRegressor()
gbdt._monotonic_cst = [1, -1]  # INC, DEC
gbdt.fit(X, y)
constrained_disp0 = plot_partial_dependence(
    gbdt, X, features=[0], line_kw={'linewidth': 4, 'label': 'constrained'},
    ax=unconstrained_disp0.axes_)
constrained_disp1 = plot_partial_dependence(
    gbdt, X, features=[1], line_kw={'linewidth': 4, 'label': 'constrained'},
    ax=unconstrained_disp1.axes_)

constrained_disp0.axes_[0, 0].plot(X[:, 0], y, 'o', alpha=.6, zorder=-1)
constrained_disp0.axes_[0, 0].set_ylim(-6, 6)
constrained_disp1.axes_[0, 0].plot(X[:, 1], y, 'o', alpha=.6, zorder=-1)
constrained_disp1.axes_[0, 0].set_ylim(-6, 6)

plt.legend()
fig.suptitle("Monotonic constraints illustration")
plt.show()
