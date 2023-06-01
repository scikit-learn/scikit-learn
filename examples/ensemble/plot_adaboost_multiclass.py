"""
=====================================
Multi-class AdaBoosted Decision Trees
=====================================

This example shows how boosting can improve prediction accuracy on a multi-class
problem by reproducing Figure 1 of [1].

- [1] J. Zhu, H. Zou, S. Rosset, T. Hastie, "Multi-class AdaBoost", 2009.
"""

# Noel Dawe <noel.dawe@gmail.com>
# License: BSD 3 clause

# %%
# The classification dataset is constructed by taking a ten-dimensional standard
# normal distribution and defining three classes separated by nested concentric
# ten-dimensional spheres such that roughly equal numbers of samples are in each
# class (quantiles of the :math:`\chi^2` distribution).
from sklearn.datasets import make_gaussian_quantiles

X, y = make_gaussian_quantiles(
    n_samples=13000, n_features=10, n_classes=3, random_state=1
)

# %%
# We're splitting the data in train and test set. 10000 samples are used for
# training, 3000 for test.
n_split = 3000

X_train, X_test = X[:n_split], X[n_split:]
y_train, y_test = y[:n_split], y[n_split:]

# %%
# Next, we compare the performance of the `SAMME` and `SAMME.R` algorithms. Note
# that in :class:`ensemble.AdaBoostClassifier` the default value for the
# `algorithm` parameter is `SAMME.R`. `SAMME.R` uses the probability estimates
# to update the additive model, while `SAMME` uses the class labels predicted by
# the base estimators (:class:`tree.DecisionTreeClassifier`) to update the
# sample weights. They both combine the estimators' predictions through weighted
# voting.
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

bdt_real = AdaBoostClassifier(
    DecisionTreeClassifier(max_depth=2),
    n_estimators=300,
    learning_rate=1,
    algorithm="SAMME.R",
)

bdt_discrete = AdaBoostClassifier(
    DecisionTreeClassifier(max_depth=2),
    n_estimators=300,
    learning_rate=1.5,
    algorithm="SAMME",
)

bdt_real.fit(X_train, y_train)
bdt_discrete.fit(X_train, y_train)

# %%
# We're calculating the accuracy score of the predicted values against the
# actual ones. As the example will illustrate, the `SAMME.R` algorithm typically
# converges faster than `SAMME`, achieving a lower test error with fewer
# boosting iterations.
from sklearn.metrics import accuracy_score

real_test_errors = []
discrete_test_errors = []

for real_test_predict, discrete_test_predict in zip(
    bdt_real.staged_predict(X_test), bdt_discrete.staged_predict(X_test)
):
    real_test_errors.append(1.0 - accuracy_score(real_test_predict, y_test))
    discrete_test_errors.append(1.0 - accuracy_score(discrete_test_predict, y_test))

# %%
# We want to extract the `estimator_errors_` and the `estimator_weights_`
# from the fitted AdaBoostClassifiers in order to plot them. Unlike `SAMME`, the
# `SAMME.R` algorithm does not calculate the `estimator_weights_`, because the
# weights are determined based on the probabilities predicted by the base
# estimators.
#
# Since boosting might terminate early, leaving some of the `n_estimators`
# estimators untrained, we want to make sure to only keep attribtes from trained
# Decision Trees. So, we crop those attribute's arrays to the actual number of
# trained trees:
n_trees_discrete = len(bdt_discrete)
n_trees_real = len(bdt_real)

discrete_estimator_errors = bdt_discrete.estimator_errors_[:n_trees_discrete]
real_estimator_errors = bdt_real.estimator_errors_[:n_trees_real]
discrete_estimator_weights = bdt_discrete.estimator_weights_[:n_trees_discrete]


# %%
# The following subplots provide insights into the performance, convergence, and
# behavior of the `SAMME` and `SAMME.R` algorithms during the boosting process.
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 5))

plt.subplot(131)
plt.title("Errors on the test set")
plt.plot(range(1, n_trees_discrete + 1), discrete_test_errors, c="black", label="SAMME")
plt.plot(
    range(1, n_trees_real + 1),
    real_test_errors,
    c="black",
    linestyle="dashed",
    label="SAMME.R",
)
plt.legend()
plt.ylim(0.18, 0.62)
plt.ylabel("Test Error")
plt.xlabel("Number of Trees")

plt.subplot(132)
plt.title("Classification error while boosting")
plt.plot(
    range(1, n_trees_discrete + 1),
    discrete_estimator_errors,
    "b",
    label="SAMME",
    alpha=0.5,
)
plt.plot(
    range(1, n_trees_real + 1), real_estimator_errors, "r", label="SAMME.R", alpha=0.5
)
plt.legend()
plt.ylabel("Error")
plt.xlabel("Number of Trees")
plt.ylim((0.2, max(real_estimator_errors.max(), discrete_estimator_errors.max()) * 1.2))
plt.xlim((-20, len(bdt_discrete) + 20))

plt.subplot(133)
plt.title("Estimator weights for SAMME algorithm")
plt.plot(range(1, n_trees_discrete + 1), discrete_estimator_weights, "b", label="SAMME")
plt.legend()
plt.ylabel("Weight")
plt.xlabel("Number of Trees")
plt.ylim((0, discrete_estimator_weights.max() * 1.2))
plt.xlim((-20, n_trees_discrete + 20))

# prevent overlapping y-axis labels
plt.subplots_adjust(wspace=0.25)
plt.show()

# %%
# The left plot shows the errors (1 - accuracy) of each algorithm on the test
# set after each boosting iteration. In the middle plot, the classification
# errors of the individual estimators during the boosting process are shown. The
# subplot to the right represents the weights assigned to each estimator
# (Decision Tree) in the `SAMME` algorithm. All trees have a weight of 1 in the
# `SAMME.R` algorithm, which is therefore not shown.
