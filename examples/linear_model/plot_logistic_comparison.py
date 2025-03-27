"""
=============================================================
Logistic Regression: From Binary to Multiclass Classification
=============================================================

This example illustrates the concept of logistic regression, progressing from
binary classification to multiclass approaches. It combines and extends the
concepts from two separate examples:

- The binary classification example that demonstrates the logistic function
- The multiclass example comparing multinomial and one-vs-rest approaches

The example starts with visualizing how logistic regression performs binary
classification using the logistic (sigmoid) function to model probabilities.
It then extends to the multiclass case, comparing multinomial and one-vs-rest
approaches by visualizing their decision boundaries and hyperplanes.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Part 1: Binary Logistic Regression
# ----------------------------------
#
# We begin with a basic illustration of binary logistic regression.
# We generate a simple synthetic dataset with two classes and visualize
# how logistic regression models the probability of class membership
# using the logistic (sigmoid) function.

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import expit

from sklearn.datasets import make_blobs
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.multiclass import OneVsRestClassifier

# %%
# Generate a binary classification dataset
# ---------------------------------------
#
# We create a simple synthetic dataset with two classes separated by a threshold.
xmin, xmax = -5, 5
n_samples = 100
np.random.seed(0)
X_binary = np.random.normal(size=n_samples)
y_binary = (X_binary > 0).astype(float)
X_binary[X_binary > 0] *= 4
X_binary += 0.3 * np.random.normal(size=n_samples)

X_binary = X_binary[:, np.newaxis]

# %%
# Fit the binary classifier
# -------------------------
#
# We train a logistic regression model and also a linear regression model
# to illustrate the difference between them.
clf_binary = LogisticRegression(C=1e5)
clf_binary.fit(X_binary, y_binary)

# %%
# Visualize the binary model
# -------------------------
#
# Plot the data points and both models to see how logistic regression applies
# the sigmoid function to model probabilities.
plt.figure(figsize=(8, 6))
plt.scatter(X_binary.ravel(), y_binary, label="Training data", color="black", zorder=20)
X_test = np.linspace(-5, 10, 300)

# Plot the logistic regression model (sigmoid curve)
loss = expit(X_test * clf_binary.coef_ + clf_binary.intercept_).ravel()
plt.plot(X_test, loss, label="Logistic Regression Model", color="red", linewidth=3)

# Plot linear regression for comparison
ols = LinearRegression()
ols.fit(X_binary, y_binary)
plt.plot(
    X_test,
    ols.coef_ * X_test + ols.intercept_,
    label="Linear Regression Model",
    linewidth=1,
)
plt.axhline(0.5, color=".5", label="Decision Boundary")

plt.ylabel("Probability / Target Value")
plt.xlabel("Feature Value")
plt.xticks(range(-5, 10))
plt.yticks([0, 0.5, 1])
plt.ylim(-0.25, 1.25)
plt.xlim(-4, 10)
plt.legend(loc="lower right")
plt.title("Binary Logistic Regression")
plt.tight_layout()

# %%
# The plot above illustrates several key concepts:
#
# - Logistic regression uses the sigmoid function to transform linear predictions
#   into probabilities between 0 and 1.
# - The decision boundary is at probability 0.5, where the sigmoid curve crosses
#   the horizontal line.
# - Unlike linear regression, logistic regression's output is constrained between
#   0 and 1, making it suitable for probability estimation.
#
# Part 2: Multiclass Logistic Regression
# --------------------------------------
#
# Now we extend to the multiclass case, where we have more than two classes.
# We'll compare two approaches:
#
# - Multinomial logistic regression: handles all classes simultaneously
# - One-vs-Rest (OvR) logistic regression: trains a binary classifier for each class

# %%
# Generate a multiclass dataset
# ----------------------------
#
# We create a synthetic dataset with three classes.
centers = [[-5, 0], [0, 1.5], [5, -1]]
X, y = make_blobs(n_samples=1_000, centers=centers, random_state=40)
transformation = [[0.4, 0.2], [-0.4, 1.2]]
X = np.dot(X, transformation)

fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(X[:, 0], X[:, 1], c=y, edgecolor="black")
ax.set(title="Multiclass Dataset", xlabel="Feature 1", ylabel="Feature 2")
ax.legend(*scatter.legend_elements(), title="Classes")
plt.tight_layout()

# %%
# Train multiclass classifiers
# ---------------------------
#
# We train both multinomial and one-vs-rest logistic regression classifiers.
logistic_regression_multinomial = LogisticRegression().fit(X, y)
logistic_regression_ovr = OneVsRestClassifier(LogisticRegression()).fit(X, y)

accuracy_multinomial = logistic_regression_multinomial.score(X, y)
accuracy_ovr = logistic_regression_ovr.score(X, y)

# %%
# Visualize decision boundaries
# ----------------------------
#
# We compare the decision boundaries produced by both methods.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

for model, title, ax in [
    (
        logistic_regression_multinomial,
        f"Multinomial Logistic Regression\n(Accuracy: {accuracy_multinomial:.3f})",
        ax1,
    ),
    (
        logistic_regression_ovr,
        f"One-vs-Rest Logistic Regression\n(Accuracy: {accuracy_ovr:.3f})",
        ax2,
    ),
]:
    DecisionBoundaryDisplay.from_estimator(
        model,
        X,
        ax=ax,
        response_method="predict",
        alpha=0.8,
    )
    scatter = ax.scatter(X[:, 0], X[:, 1], c=y, edgecolor="k")
    legend = ax.legend(*scatter.legend_elements(), title="Classes")
    ax.add_artist(legend)
    ax.set_title(title)

plt.tight_layout()

# %%
# The decision boundaries show how the two approaches differ:
#
# - Multinomial logistic regression considers all classes simultaneously
# - One-vs-rest trains a separate binary classifier for each class
#
# These differences in approach can lead to different decision boundaries,
# especially in regions where classes overlap.


# %%
# Visualize hyperplanes
# --------------------
#
# To better understand the differences, we visualize the hyperplanes that represent
# the boundaries where the probability for a class equals 0.5.
def plot_hyperplanes(classifier, X, ax):
    xmin, xmax = X[:, 0].min(), X[:, 0].max()
    ymin, ymax = X[:, 1].min(), X[:, 1].max()
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))

    if isinstance(classifier, OneVsRestClassifier):
        coef = np.concatenate([est.coef_ for est in classifier.estimators_])
        intercept = np.concatenate([est.intercept_ for est in classifier.estimators_])
    else:
        coef = classifier.coef_
        intercept = classifier.intercept_

    for i in range(coef.shape[0]):
        w = coef[i]
        a = -w[0] / w[1]
        xx = np.linspace(xmin, xmax)
        yy = a * xx - (intercept[i]) / w[1]
        ax.plot(xx, yy, "--", linewidth=3, label=f"Class {i}")

    return ax.get_legend_handles_labels()


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

for model, title, ax in [
    (
        logistic_regression_multinomial,
        "Multinomial Logistic Regression Hyperplanes",
        ax1,
    ),
    (logistic_regression_ovr, "One-vs-Rest Logistic Regression Hyperplanes", ax2),
]:
    hyperplane_handles, hyperplane_labels = plot_hyperplanes(model, X, ax)
    scatter = ax.scatter(X[:, 0], X[:, 1], c=y, edgecolor="k")
    scatter_handles, scatter_labels = scatter.legend_elements()

    all_handles = hyperplane_handles + scatter_handles
    all_labels = hyperplane_labels + scatter_labels

    ax.legend(all_handles, all_labels, title="Classes")
    ax.set_title(title)

plt.tight_layout()
plt.show()

# %%
# Key differences and takeaways:
# -----------------------------
#
# 1. Binary vs. Multiclass:
#    - Binary logistic regression models probabilities using the sigmoid function
#    - Multiclass logistic regression extends this to multiple classes
#
# 2. Multinomial vs. One-vs-Rest:
#    - Multinomial optimizes all class probabilities simultaneously
#    - One-vs-Rest trains independent binary classifiers
#    - Multinomial often produces more consistent decision boundaries
#
# 3. Practical considerations:
#    - For binary problems, both approaches are equivalent
#    - For multiclass, multinomial is generally preferred as it produces
#      better-calibrated probabilities and more interpretable results
#    - The choice may depend on the specific problem and computational constraints
