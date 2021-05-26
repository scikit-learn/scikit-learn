"""
===================
Quantile regression
===================
This example illustrates how quantile regression can predict non-trivial
conditional quantiles.

The left figure shows the case when the error distribution is normal,
but has non-constant variance, i.e. with heteroscedasticity.

The right figure shows an example of an asymmetric error distribution,
namely the Pareto distribution.
"""
print(__doc__)
# Authors: David Dale <dale.david@mail.ru>
#          Christian Lorentzen <lorentzen.ch@gmail.com>
# License: BSD 3 clause
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import QuantileRegressor, LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_val_score


def plot_points_highlighted(x, y, model_low, model_high, ax):
    """Plot points with highlighting."""
    mask = y <= model_low.predict(X)
    ax.scatter(x[mask], y[mask], c="k", marker="x")
    mask = y > model_high.predict(X)
    ax.scatter(x[mask], y[mask], c="k", marker="x")
    mask = (y > model_low.predict(X)) & (y <= model_high.predict(X))
    ax.scatter(x[mask], y[mask], c="k")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5), sharey=True)

rng = np.random.RandomState(42)
x = np.linspace(0, 10, 100)
X = x[:, np.newaxis]
y = 10 + 0.5 * x + rng.normal(loc=0, scale=0.5 + 0.5 * x, size=x.shape[0])
y_mean = 10 + 0.5 * x
ax1.plot(x, y_mean, "k--")

quantiles = [0.05, 0.5, 0.95]
models = []
for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, alpha=0)
    qr.fit(X, y)
    ax1.plot(x, qr.predict(X))
    models.append(qr)

plot_points_highlighted(x, y, models[0], models[2], ax1)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("Quantiles of heteroscedastic Normal distributed target")
ax1.legend(["true mean"] + quantiles)


a = 5
y = 10 + 0.5 * x + 10 * (rng.pareto(a, size=x.shape[0]) - 1 / (a - 1))
ax2.plot(x, y_mean, "k--")

models = []
for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, alpha=0)
    qr.fit(X, y)
    ax2.plot([0, 10], qr.predict([[0], [10]]))
    models.append(qr)

plot_points_highlighted(x, y, models[0], models[2], ax2)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("Quantiles of asymmetric Pareto distributed target")
ax2.legend(["true mean"] + quantiles, loc="lower right")
ax2.yaxis.set_tick_params(labelbottom=True)

plt.show()

# %%
# Note that both targets have the same mean value, indicated by the dashed
# black line. As the Normal distribution is symmetric, mean and median are
# identical and the predicted 0.5 quantile almost hits the true mean.
# In the Pareto case, the difference between predicted median and true mean
# is evident. We also marked the points below the 0.05 and above 0.95
# predicted quantiles by small crosses. You might count them and consider
# that we have 100 samples in total.
#
# The second part of the example shows that LinearRegression minimizes MSE
# in order to predict the mean, while QuantileRegressor with `quantile=0.5`
# minimizes MAE in order to predict the median. Both do their own job well.

models = [LinearRegression(), QuantileRegressor(alpha=0)]
names = ["OLS", "Quantile"]

print("# In-sample performance")
for model_name, model in zip(names, models):
    print(model_name + ":")
    model.fit(X, y)
    mae = mean_absolute_error(model.predict(X), y)
    rmse = np.sqrt(mean_squared_error(model.predict(X), y))
    print(f"MAE = {mae:.4}  RMSE = {rmse:.4}")
print("\n# Cross-validated performance")
for model_name, model in zip(names, models):
    print(model_name + ":")
    mae = -cross_val_score(model, X, y, cv=3,
                           scoring="neg_mean_absolute_error").mean()
    rmse = np.sqrt(-cross_val_score(model, X, y, cv=3,
                                    scoring="neg_mean_squared_error").mean())
    print(f"MAE = {mae:.4}  RMSE = {rmse:.4}")
