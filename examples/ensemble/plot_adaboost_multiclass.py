"""
=====================================
Multi-class AdaBoosted Decision Trees
=====================================

This example shows how boosting can improve prediction accuracy on a multi-class
problem. It reproduces the experiment depicted by Figure 1 in [1].

.. topic:: References:
    .. [1] J. Zhu, H. Zou, S. Rosset, T. Hastie,
       `"Multi-class AdaBoost"
       <https://hastie.su.domains/Papers/samme.pdf>`_,
       2009.
"""

# Noel Dawe <noel.dawe@gmail.com>
# License: BSD 3 clause

# %%
# Creating the dataset
# --------------------
# The classification dataset is constructed by taking a ten-dimensional standard
# normal distribution and defining three classes separated by nested concentric
# ten-dimensional spheres such that roughly equal numbers of samples are in each
# class (quantiles of the :math:`\chi^2` distribution).
from sklearn.datasets import make_gaussian_quantiles

X, y = make_gaussian_quantiles(
    n_samples=13_000, n_features=10, n_classes=3, random_state=1
)

# %%
# We split the dataset into 2 sets: 10,000 samples are used for training and
# 3,000 for test.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=10_000, shuffle=False
)

# %%
# Training `AdaBoostClassifier` with two different algorithms
# -----------------------------------------------------------
# Next, we compare the performance of the `SAMME` and `SAMME.R` algorithms. Note
# that in :class:`~sklearn.ensemble.AdaBoostClassifier` the default value for
# the `algorithm` parameter is `SAMME.R`. `SAMME.R` uses the probability
# estimates to update the additive model, while `SAMME` uses the class labels
# predicted by the base estimators (i.e.
# :class:`~sklearn.tree.DecisionTreeClassifier`) to update the sample weights.
# They both combine the estimators' predictions through weighted voting.
#
# Here, we define the weak learner as a
# :class:`~sklearn.tree.DecisionTreeClassifier` and set the maximum number of
# leaves to 4. In a real setting, this parameter should be tuned. We set it to a
# rather low value to limit the runtime of the example.
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
import time

weak_learner = DecisionTreeClassifier(max_depth=1, max_leaf_nodes=4)

start = time.time()
samme_r = AdaBoostClassifier(
    estimator=weak_learner,
    n_estimators=300,
    algorithm="SAMME.R",
    random_state=42,
).fit(X_train, y_train)
end = time.time()
print(f"Training time of SAMME.R algorithm: {end - start:.3f}s")

start = time.time()
samme = AdaBoostClassifier(
    estimator=weak_learner,
    n_estimators=300,
    algorithm="SAMME",
    random_state=42,
).fit(X_train, y_train)
end = time.time()
print(f"Training time of SAMME algorithm: {end - start:.3f}s")

# %%
# Comparing the Results
# ---------------------
# We calculate the misclassification error of the additive models at each
# boosting iteration to assess their performance. To obtain the predictions at
# each iteration, we use the
# :func:`~sklearn.esamble._weight_boosting.staged_predict` method, which returns
# a generator containing the predictions of all the weak learners as the
# training progresses. The final prediction is determined by majority voting
# among the weak learners. We compare these predictions to the actual values in
# `y_test` by calculating the difference. This process allows us to generate the
# data needed for plotting a learning curve."

from sklearn.metrics import accuracy_score


def misclassification_error(y_true, y_pred):
    return 1 - accuracy_score(y_true, y_pred)


samme_r_errors = [
    misclassification_error(y_test, y_pred) for y_pred in samme_r.staged_predict(X_test)
]
samme_errors = [
    misclassification_error(y_test, y_pred) for y_pred in samme.staged_predict(X_test)
]

# %%
# For performance review, we also want to extract the `estimator_errors_` and
# the `estimator_weights_` from the fitted AdaBoostClassifiers.
# `estimator_errors_` is an array of weighted errors of each individual weak
# learner immediately after it has been trained at each boosting iteration.
# `estimator_weights_` are only calculated in the `SAMME` algorithm, because in
# `SAMME.R` the estimator weights are updated using the predicted class
# probabilities instead of the misclassification rate
#
# Since boosting might terminate early, leaving some of the `n_estimators`
# estimators untrained, we want to make sure to only keep attribtes from trained
# Decision Trees. So, we crop those attribute's arrays to the actual number of
# trained trees:
n_trees_discrete = len(samme)
n_trees_real = len(samme_r)

samme_estimator_errors = samme.estimator_errors_[:n_trees_discrete]
samme_r_estimator_errors = samme_r.estimator_errors_[:n_trees_real]
samme_r_estimator_weights = samme.estimator_weights_[:n_trees_discrete]


# %%
# The following subplots provide insights into the performance, convergence, and
# behavior of the `SAMME` and `SAMME.R` algorithms during the boosting process.
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 5))

plt.subplot(131)
plt.title("Classification errors on the test set")
plt.plot(range(1, n_trees_discrete + 1), samme_errors, c="black", label="SAMME")
plt.plot(
    range(1, n_trees_real + 1),
    samme_r_errors,
    c="black",
    linestyle="dashed",
    label="SAMME.R",
)
plt.legend()
plt.ylim(0.18, 0.62)
plt.ylabel("Test Error")
plt.xlabel("Iterations")

plt.subplot(132)
plt.title("Classification errors of individual weak learners \n on the training set")
plt.plot(
    range(1, n_trees_discrete + 1),
    samme_estimator_errors,
    "b",
    label="SAMME",
    alpha=0.5,
)
plt.plot(
    range(1, n_trees_real + 1),
    samme_r_estimator_errors,
    "r",
    label="SAMME.R",
    alpha=0.5,
)
plt.legend()
plt.ylabel("Error")
plt.xlabel("Iterations")
plt.ylim((0.2, max(samme_r_estimator_errors.max(), samme_estimator_errors.max()) * 1.2))
plt.xlim((-20, len(samme) + 20))

plt.subplot(133)
plt.title("Estimator weights of individual weak learners \n for SAMME algorithm")
plt.plot(range(1, n_trees_discrete + 1), samme_r_estimator_weights, "b", label="SAMME")
plt.legend()
plt.ylabel("Weight")
plt.xlabel("Iterations")
plt.ylim((0, samme_r_estimator_weights.max() * 1.2))
plt.xlim((-20, n_trees_discrete + 20))

# prevent overlapping y-axis labels
plt.subplots_adjust(wspace=0.25)
plt.show()

# %%
# The left plot shows the majoritarian classification errors (1 - accuracy)
# of each algorithm on the test set after each boosting iteration. In the middle
# plot, the classification errors of the individual estimators during the
# training are shown. The subplot to the right represents the weights assigned
# to each estimator in the `SAMME` algorithm. For the `SAMME.R`
# algorithm all trees have a weight of 1, which is therefore not plotted.
