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
f_0 = rng.rand(n_samples)  # positive correlation with y
f_1 = rng.rand(n_samples)  # negative correlation with y
X = np.c_[f_0, f_1]
print(X.shape)
noise = rng.normal(loc=0.0, scale=0.01, size=n_samples)
y = (5 * f_0 + np.sin(10 * np.pi * f_0) -
     5 * f_1 - np.cos(10 * np.pi * f_1) +
     noise)

fig, ax = plt.subplots()


gbdt = HistGradientBoostingRegressor()
gbdt.fit(X, y)
disp = plot_partial_dependence(
    gbdt, X, features=[0, 1],
    line_kw={'linewidth': 4, 'label': 'unconstrained'},
    ax=ax)

gbdt = HistGradientBoostingRegressor()
gbdt._monotonic_cst = [1, -1]  # INC, DEC
gbdt.fit(X, y)
plot_partial_dependence(
    gbdt, X, features=[0, 1],
    feature_names=('First feature\nPositive constraint',
                   'Second feature\nNegtive constraint'),
    line_kw={'linewidth': 4, 'label': 'constrained'},
    ax=disp.axes_)

for f_idx in (0, 1):
    disp.axes_[0, f_idx].plot(X[:, f_idx], y, 'o', alpha=.3, zorder=-1)
    disp.axes_[0, f_idx].set_ylim(-6, 6)

plt.legend()
fig.suptitle("Monotonic constraints illustration")

plt.show()
