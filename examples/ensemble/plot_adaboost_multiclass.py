"""
=====================================
Multi-class AdaBoosted Decision Trees
=====================================

This example shows how boosting can improve the prediction accuracy on a
multi-class problem. It reproduces a similar experiment as depicted by Figure 1
in [1]_.

.. topic:: References:

    .. [1] :doi:`J. Zhu, H. Zou, S. Rosset, T. Hastie, "Multi-class adaboost."
           Statistics and its Interface 2.3 (2009): 349-360.
           <10.4310/SII.2009.v2.n3.a8>`

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
    n_samples=2_000, n_features=10, n_classes=3, random_state=1
)

# %%
# We split the dataset into 2 sets: 70 percent of the samples are used for
# training and the remaining 30 percent for testing.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=0.7, random_state=42
)

# %%
# Training `AdaBoostClassifier`
# -----------------------------
# We train the estimator, that utilizes boosting to improve classification
# accuracy. Boosting is a method designed to train weak learners (base
# estimators) that learn from their predecessor's mistakes.
#
# Here, we define the weak learner as a
# :class:`~sklearn.tree.DecisionTreeClassifier` and set the maximum number of
# leaves to 8. In a real setting, this parameter should be tuned. We set it to a
# rather low value to limit the runtime of the example.
#
# The `SAMME` algorithm build into the
# :class:`~sklearn.ensemble.AdaBoostClassifier` then uses the class labels
# predicted by the base estimators to update the sample weights. The base
# estimators' predictions are combined through weighted voting.
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier

weak_learner = DecisionTreeClassifier(max_leaf_nodes=8)
n_estimators = 300

adaboost_clf = AdaBoostClassifier(
    estimator=weak_learner,
    n_estimators=n_estimators,
    algorithm="SAMME",
    random_state=42,
).fit(X_train, y_train)
# %%
# Analysis
# --------
# Convergence of the AdaBoostClassifier
# *************************************
# To demonstrate how boosting can increase the accuracy, we compare the
# `misclassification_error` of the :class:`~sklearn.tree.DecisionTreeClassifier`,
# as a baseline score, with the `misclassification_error` of the boosted Trees.
from sklearn.metrics import accuracy_score


def misclassification_error(y_true, y_pred):
    return 1 - accuracy_score(y_true, y_pred)


weak_learners_misclassification_error = misclassification_error(
    y_test, weak_learner.fit(X_train, y_train).predict(X_test)
)

print(
    "DecisionTreeClassifier's misclassification_error:"
    f" {weak_learners_misclassification_error}"
)

# %%
# After training the :class:`~sklearn.tree.DecisionTreeClassifier` model, the
# achieved error, which is less than 0.66, surpasses the expected value that
# would have been obtained by random guessing alone.
#
# Now, we calculate the misclassification error of the additive model at each
# boosting iteration to assess it's performance. To obtain the predictions at
# each iteration, we use the
# :func:`~sklearn.ensemble._weight_boosting.staged_predict` method which returns
# a generator containing the predictions of all the weak learners during the
# training progresses. The final prediction is determined by majority vote among
# the weak learners. We compare these predictions to the actual values in
# `y_test` by calculating the difference. This process allows us to generate the
# data needed for plotting a learning curve.
import pandas as pd

boosting_errors = pd.DataFrame(
    {
        "Number of trees": range(1, n_estimators + 1),
        "AdaBoost": [
            misclassification_error(y_test, y_pred)
            for y_pred in adaboost_clf.staged_predict(X_test)
        ],
    }
).set_index("Number of trees")
ax = boosting_errors.plot()
ax.set_ylabel("Misclassification error on test set")
ax.set_title("Convergence of AdaBoost algorithm")

# %%
# The plot shows the classification error, i.e. (1 - accuracy), on the test
# set after each boosting iteration. We can see that the error of the boosted
# trees converges at approximately 0.3 after 50 trees has been trained, thus
# archiving much higher accuracy than training a single tree.
#
# The misclassification error jitters because it uses the discrete outputs of
# the weak learners to train the boosted model.
#
# The convergence of :class:`~sklearn.ensemble.AdaBoostClassifier` is mainly
# influenced by the `learning_rate`, the number of weak learners used
# (`n_estimators`), and the expressivity of the weak learners (e.g.
# `max_leaf_nodes`).

# %%
# Errors of the Weak Learners
# ***************************
# We also want to extract the `estimator_errors_` and the `estimator_weights_`
# from the fitted :class:`~sklearn.ensemble.AdaBoostClassifier`.
# `estimator_errors_` is an array of weighted errors of each individual weak
# learner immediately after it has been trained at each boosting iteration,
# while the `estimator_weights_` represents the corresponding weights.
#
# Since boosting might terminate early, leaving some of the `n_estimators`
# estimators untrained, we want to make sure to only keep attribtes from trained
# Decision Trees. So, we crop those attribute's arrays to the actual number of
# trained trees (using `[:len(adaboost_clf)]`):

weak_learners_info = pd.DataFrame(
    {
        "Number of trees": range(1, n_estimators + 1),
        "Errors": adaboost_clf.estimator_errors_[: len(adaboost_clf)],
        "Weights": adaboost_clf.estimator_weights_[: len(adaboost_clf)],
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
fig.suptitle("Weak learner's errors and weights for the AdaBoostClassifier")
fig.tight_layout()

# %%
# In the left plot, the classification errors of the individual estimators
# during the training are shown. The plot to the right represents the weights
# assigned to each base estimator.
