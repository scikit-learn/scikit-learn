"""
==================================================
Using Winsorizer to deal with outliers in features
==================================================

Many estimation techniques are highly sensitive to the presence of outliers in
the data. Winsorization can be an effective way of dealing with outliers.

This example demonstrates the use of :class:`preprocessing.Winsorization` to
clip outliers in feature data. Note that :class:`preprocessing.Winsorization`
follows the ``Transformer`` API so it can be used as part of a preprocessing
:class:`sklearn.pipeline.Pipeline`). However, because it follows the
``Transformer`` API, it does not allow for winsorization of target variables.

Ordinary least squares regression (:class:`linear_model.LinearRegression`) is
particularly sensitive to outliers. Winsorization can be used as a
preprocessing step to improve the robustness of a linear model fit using
this method. One alternative way to deal with outliers when fitting linear
models to use :class:`linear_model.HuborRegressor`. A downside of doing this
is the fit procedure for :class:`linear_model.HuborRegressor` can be more
computationally intensive because it uses BFGS, an iterative method, to solve
a nonlinear optimization problem. This is especially the case when the number
of samples in the data set is much larger than the number of features.
"""

print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Winsorizer

# Generate toy data.
rng = np.random.RandomState(0)
X, y = make_regression(n_samples=50, n_features=1, random_state=0, noise=1.0,
                       bias=100.0)

# Add some strong outliers to the dataset.
X_outliers = rng.normal(0, 2.0, size=(5, 1))
X[-5:, :] = X_outliers

# import pdb; pdb.set_trace()
plt.plot(X, y, 'b.')

raw_ols = LinearRegression()
raw_ols.fit(X, y)

wins_ols = make_pipeline(Winsorizer(0.10), LinearRegression())
wins_ols.fit(X, y)

x = np.linspace(X.min(), X.max(), 7)

ols1 = raw_ols
ols2 = wins_ols.named_steps['linearregression']

fit1 = ols1.coef_ * x + ols1.intercept_
fit2 = ols2.coef_ * x + ols2.intercept_

plt.plot(x, fit1, label='OLS regression')
plt.plot(x, fit2, label='OLS regression with winsorization')

plt.title('Linear regression with and without winsorization')
plt.xlabel('X')
plt.ylabel('y')
plt.legend(loc=0)
plt.show()
