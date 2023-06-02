"""
=====================================
Multi-class AdaBoosted Decision Trees
=====================================

This example shows how boosting can improve prediction accuracy on a multi-class
problem. It reproduces the experiment depicted by Figure 1 in [1]_.

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

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=10_000)

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
# leaves to 8. In a real setting, this parameter should be tuned. We set it to a
# rather low value to limit the runtime of the example.
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
import time

weak_learner = DecisionTreeClassifier(max_depth=1, max_leaf_nodes=8)
n_estimators = 600

start = time.time()
samme_r = AdaBoostClassifier(
    estimator=weak_learner,
    n_estimators=n_estimators,
    algorithm="SAMME.R",
    random_state=42,
).fit(X_train, y_train)
end = time.time()
print(f"Training time of SAMME.R algorithm: {end - start:.3f}s")

start = time.time()
samme = AdaBoostClassifier(
    estimator=weak_learner,
    n_estimators=n_estimators,
    algorithm="SAMME",
    random_state=42,
).fit(X_train, y_train)
end = time.time()
print(f"Training time of SAMME algorithm: {end - start:.3f}s")

# %%
# Analysis
# --------
# Convergence of the additive AdaBoostClassifier
# **********************************************
# We calculate the misclassification error of the additive models at each
# boosting iteration to assess their performance. To obtain the predictions at
# each iteration, we use the
# :func:`~sklearn.esamble._weight_boosting.staged_predict` method, which returns
# a generator containing the predictions of all the weak learners during the
# training progresses. The final prediction is determined by majority vote
# among the weak learners. We compare these predictions to the actual values in
# `y_test` by calculating the difference. This process allows us to generate the
# data needed for plotting a learning curve.

from sklearn.metrics import accuracy_score


def misclassification_error(y_true, y_pred):
    return 1 - accuracy_score(y_true, y_pred)


import pandas as pd

boosting_errors = pd.DataFrame(
    {
        "Number of trees": range(1, n_estimators + 1),
        "SAMME.R": [
            misclassification_error(y_test, y_pred)
            for y_pred in samme_r.staged_predict(X_test)
        ],
        "SAMME": [
            misclassification_error(y_test, y_pred)
            for y_pred in samme.staged_predict(X_test)
        ],
    }
).set_index("Number of trees")
ax = boosting_errors.plot()
ax.set_ylabel("Misclassification error on test set")
ax.set_title("Convergence of AdaBoost algorithms")

# %%
# The plot shows the majoritarian classification errors (1 - accuracy) of each
# algorithm on the test set after each boosting iteration.
#
# The misclassification error jitters more with the `SAMME` algorithm than with
# `SAMME.R`, because it uses the discrete outputs of the weak learners to train
# the boosted model.
#
# The convergence of AdaBoostClassifier is mainly influenced by the
# `learning_rate`, the number of weak learners used (`n_estimators`), and the
# expressivity of the weak learners (e.g. `max_leaf_nodes`).

# %%
# Errors of the Weak Learners
# ***************************
# We also want to extract the `estimator_errors_` and
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
# trained trees (using `[:len(samme)]`):

weak_learners_info = pd.DataFrame(
    {
        "Number of trees": range(1, n_estimators + 1),
        "Errors": samme.estimator_errors_[: len(samme)],
        "Weights": samme.estimator_weights_[: len(samme)],
    }
).set_index("Number of trees")

axs = weak_learners_info.plot(
    subplots=True, layout=(1, 2), figsize=(10, 4), legend=False, color="tab:blue"
)
axs[0, 0].set_ylabel("Train error")
axs[0, 0].set_title("Weak learner's training error")
axs[0, 1].set_ylabel("Weight")
axs[0, 1].set_title("Weak learner's weight")
fig = axs[0, 0].get_figure()
fig.suptitle("Weak learner's errors and weights for SAMME")
fig.tight_layout()

# %%
# In the left plot, the classification errors of the individual estimators
# during the training are shown. The plot to the right represents the weights
# assigned to each estimator in the `SAMME` algorithm.
