"""
=========================================================
Ridge coefficients as a function of the L2 Regularization
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

The regularized loss function is given by: :math:`\\mathcal{L}(X, y, \\beta) = \\|
y - X \\beta \\|^{2}_{2} + \\alpha \\| \\beta \\|^{2}_{2}`

where :math:`X` is the input data, :math:`y` is the target variable,
:math:`\\beta` is the vector of coefficients associated with the features, and
:math:`\\alpha` is the regularization strength.

The regularized loss function aims to balance the trade-off between accurately
predicting the training set and to prevent overfitting.

In this regularized loss, the left-hand side (:math:`\\|y - X\\beta\\|^{2}_{2}`)
measures the squared difference between the actual target variable, :math:`y`,
and the predicted values. Minimizing this term alone could lead to overfitting,
as the model may become too complex and sensitive to noise in the training data.

To address overfitting, Ridge regularization adds a constraint, called a penalty
term, (:math:`\\alpha \\| \\beta\\|^{2}_{2}`) to the loss function. This penalty
term is the sum of the squares of the model's coefficients, multiplied by the
regularization strength :math:`\\alpha`. By introducing this constraint, Ridge
regularization discourages any single coefficient :math:`\\beta_{i}` from taking
an excessively large value and encourages smaller and more evenly distributed
coefficients. Higher values of :math:`\\alpha` force the coefficients towards
zero. However, an excessively high :math:`\\alpha` can result in an underfit
model that fails to capture important patterns in the data.

Therefore, the regularized loss function combines the prediction accuracy term
and the penalty term. By adjusting the regularization strength, practitioners
can fine-tune the degree of constraint imposed on the weights, training a model
capable of generalizing well to unseen data while avoiding overfitting.

For a comparison between Ridge and Huber Regularization see
:doc:`~examples.linear_model.plot_huber_vs_ridge`
"""

# Author: Kornel Kielczewski -- <kornel.k@plusnet.pl>

# %%
# Purpose of this example
# -----------------------
# For the purpose of showing how Ridge regularization works, we will create two
# data sets: non-noisy and noisy. Then we will train them on a range of
# regularization strengths (:math:`\alpha`) and plot how the trained
# coefficients and the mean squared error between those and the original values
# behave as functions of the regularization strength.

# %%
# Ridge Regularization on non-noisy data set
# ------------------------------------------

# %%
# Creating non-noisy data set
#############################
# We make a toy data set, that's suitable to detect regression. Please note that
# in this example the data is non-noisy, hence it is possible to extract the
# exact coefficients.
from sklearn.datasets import make_regression

X, y, w = make_regression(n_samples=100, n_features=10, coef=True, random_state=1)

# Obtain the true coefficients
print(f"The true coefficient of this regression problem are:\n{w}")

# %%
# Training Ridge Regressor
##########################
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
# Note: In a real setting you should use :class:`~sklearn.linear_model.RidgeCV`
# for training a Ridge regression with built-in cross-validation, but for the
# purpose of demonstration we will keep it simple here.

# %%
# Plotting trained Coefficients and Mean Squared Errors
#######################################################
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

# Interpreting the plots
########################
# The plot on the left-hand side shows how the regularization strength (`alpha`)
# affects the Ridge regression coefficients. Smaller values of `alpha` (weak
# regularization), allow the coefficients to closely resemble the true
# coefficients (`w`) used to generate the data set. This is because no additional
# noise was added to our artificial data set. As `alpha` increases, the
# coefficients shrink towards zero, gradually reducing the impact of the
# features that were formerly more significant.
#
# The right-hand side plot shows the `mean squared error` between the
# coefficients found by the model and the true coefficients (`w`). It provides a
# measure that relates to how exact our ridge model is in comparison to the true
# generative model. A low error means that it found coefficients closer to the
# ones of the true generative model. In this case, since our toy data set was
# non-noisy, we can see that the least regularized model retrieves coefficients
# closest to the true coefficients (`w`) (error is close to 0).

# %%
# Ridge Regularization on noisy data set
# --------------------------------------

# %%
# Creating noisy data set
#########################
# As a more realistic example, we generate an artificial data set with 100
# samples and 10 features. Out of these features, 3 are informative and
# contribute to the regression, while the remaining 7 features do not have any
# effect on the target variable (their true coefficients are 0). To introduce
# variability and mimic real-world scenarios, we add Gaussian noise with a
# standard deviation of 100 to each coefficient.
X, y, w = make_regression(
    n_samples=100, n_features=10, n_informative=3, noise=100, coef=True, random_state=1
)

# Obtain the true coefficients
print(f"The true coefficient of this regression problem are:\n{w}")


# %%
# Training Ridge Regressor
##########################
# Just as in the first example, we use :class:`~sklearn.linear_model.Ridge`, a
# linear model with L2 regularization. The model parameter `alpha` is a positive
# constant that multiplies the penalty term, controlling the regularization
# strength.
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
# Note: In a real setting you should use :class:`~sklearn.linear_model.RidgeCV`
# for training a Ridge regression with built-in cross-validation, but for the
# purpose of demonstration we will keep it simple here.

# %%
# Plotting trained Coefficients and Mean Squared Errors
#######################################################
# Like in the first example, we now plot the 10 different regularized
# coefficients as a function of the regularization parameter `alpha` where each
# color represents a different coefficient.
#
# On the right-hand-side, we plot how the errors of the coefficients from the
# estimator change as a function of regularization.
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
# Interpreting the plots
########################
# The plot on the left-hand side shows how the regularization strength (`alpha`)
# affects the Ridge regression coefficients. Smaller values of `alpha` (weak
# regularization), don't affect the trained coefficients very much. Comparing
# those with the true coefficients `w` from the noisy data set however reveals
# that the model hardly captures the true weights. As `alpha` increases, the
# coefficients shrink towards zero.
#
# The right-hand side plot shows the `mean squared error` between the
# coefficients found by the model and the true coefficients (`w`). It provides a
# measure that relates to how exact our ridge model is in comparison to the true
# generative model. A low error means that it found coefficients closer to the
# ones of the true generative model. In this case, we can see that the curve has
# a minimum at `alpha = 10 ** 1.7` (about 50), where the `mean squared error` is
# lowest. This means, that this particular regularization strength recovers the
# true coefficients best. A higher regularization however would cause the model
# to underfit.

# %%
# Conclusion
# ----------
# The four plots visualize the trade-off between fitting the training data too
# closely and fitting it too loosely in order to archive the best performance in
# regularized models. When `alpha` is small, the model captures the intricate
# details of the training data, whether those were caused by noise or by actual
# information. As `alpha` increases, the highest coefficients shrink more
# rapidly, rendering their corresponding features less influential in the
# training process. This can enhance the model's ability to generalize to unseen
# data, but it also poses the risk of loosing performance if the regularization
# becomes too strong.
#
# In real-world scenarios where data typically includes noise, selecting an
# appropriate `alpha` value becomes crucial in striking a balance between an
# overfitting and an underfitting model.

# %%
