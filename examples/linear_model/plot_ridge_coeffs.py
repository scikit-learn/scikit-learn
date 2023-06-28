"""
=========================================================
Ridge coefficients as a function of the L2 regularization
=========================================================

A model that overfits learns the training data too well, capturing both the
underlying patterns and the noise in the data. However, when applied to unseen
data, the learned associations may not hold. We normally detect this when we
apply our trained predictions to the test data and see the performance drop
significantly compared to the training data.

One way to overcome overfitting is through regularization. Regularization
penalizes large training weights, forcing the model to shrink certain
coefficients. Regularization reduces a model's reliance on specific information
obtained from the training samples.

This example illustrates how L2 regularization in a
:class:`~sklearn.linear_model.Ridge` regression affects a model's performance by
adding a penalty term to the Loss that increases with the coefficients (β).

The regularized loss function is given by:

:math:`RegularizedLoss = Loss(X, y, β) + Penalty(β)`

where:

- :math:`X` is the input data,

- :math:`y` is the target variable,

- :math:`β` is the vector of coefficients or weights associated with the
  features.

In Ridge regularization, the penalty term is proportional to the sum of the
squares of the coefficients. In effect, it discourages any single coefficient βi
from taking an excessively large value. This indirectly influences the
individual coefficients towards smaller values, forcing the model to shrink
certain coefficients. Thus, Ridge regularization is promoting model simplicity
or feature selection to prevent overfitting.
"""

# Author: Kornel Kielczewski -- <kornel.k@plusnet.pl>

# %%
# Creating Dataset
# ----------------
# We make a toy dataset, that's suitable to detect regression. Please note that
# in this example the data is non-noisy, hence it is possible to extract the
# exact coefficients.
from sklearn.datasets import make_regression

X, y, w = make_regression(n_samples=100, n_features=10, coef=True, random_state=1)
# Obtain the true coefficients
print(f"The true coefficient of this regression problem are:\n{w}")

# %%
# Training Ridge Regressor
# ------------------------
# We use :class:`~sklearn.linear_model.Ridge`, a linear model with L2
# regularization. The model parameter `alpha` is a positive constant that
# multiplies the penalty term, controlling the regularization strength.
import numpy as np

from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error

clf = Ridge()

# Generate values for `alpha` that are evenly distributed on a logarithmic scale
alphas = np.logspace(-3, 4, 200)
coefs = []
errors_coefs = []

# Train the model with different regularisation strengths
for a in alphas:
    clf.set_params(alpha=a).fit(X, y)
    coefs.append(clf.coef_)
    errors_coefs.append(mean_squared_error(clf.coef_, w))

# %%
# Plotting Results
# ----------------
# We now plot the 10 different regularized coefficients as a function of the
# regularization parameter `alpha` where each color represents a different
# coefficient.
#
# On the right-hand-side, we plot how the errors of the coefficients from the
# estimator change as a function of regularization.
import matplotlib.pyplot as plt
import pandas as pd

alphas = pd.Index(alphas, name="alpha")
coefs = pd.DataFrame(coefs, index=alphas, columns=[f"Feature {i}" for i in range(10)])
errors = pd.Series(errors_coefs, index=alphas, name="Mean squared error")

fig, axs = plt.subplots(1, 2, figsize=(20, 6))

coefs.plot(
    ax=axs[0],
    logx=True,
    ylabel="Ridge coefficient values",
    title="Ridge coefficients as a function of the regularization strength",
)
_ = errors.plot(
    ax=axs[1],
    logx=True,
    ylabel="Mean squared error",
    title="Coefficient error as a function of the regularization strength",
)

# %%
# Interpretation
# --------------
# The plot on the left-hand side shows how the regularization strength (`alpha`)
# affects the Ridge regression coefficients. Smaller values of `alpha` (weak
# regularization), allow the coefficients to closely resemble the true
# coefficients (`w`) used to generate the dataset. This is because no additional
# noise was added to our artificial data set. As `alpha` increases, the
# coefficients shrink towards zero, gradually reducing the formerly more
# significant features' impact.
#
# The right-hand side plot shows the `mean squared error` between the
# coefficients found by the model and the true coefficients (`w`). It provides a
# measure that relates to how exact our ridge model is in comparison to the true
# generative model. A low error means that it found coefficients closer to the
# ones of the true generative model. In this case, since our toy dataset was
# non-noisy, we can see that the least regularized model retrieves coefficients
# closest to the true coefficients (`w`) (error is close to 0).
#
# These plots show the trade-off between model simplicity (smaller coefficients,
# feature selection) and model performance. When `alpha` is small, the model
# captures the intricate details of the training data, potentially achieving
# higher performance (if there is little noise) or worse performance, (if there
# is a lot of noise). As `alpha` increases, the model becomes simpler by
# shrinking the coefficients, which can lead to a sacrifice in performance (if
# there is little noise) or to a better performance on unseen data (if there is
# a lot of noise).
#
# In our example with a noise-free toy dataset, since regularization aims to
# reduce the impact of noise, setting `alpha = 0` best recovers the true
# coefficients. However, in real-world scenarios where data typically contains
# noise, a careful selection of the regularization parameter `alpha` is
# necessary to balance model simplicity and model performance on unseen data.

# %%
