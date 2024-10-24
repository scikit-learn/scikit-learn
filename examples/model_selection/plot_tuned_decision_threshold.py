"""
======================================================
Post-hoc tuning the cut-off point of decision function
======================================================

Once a binary classifier is trained, the :term:`predict` method outputs class label
predictions corresponding to a thresholding of either the :term:`decision_function` or
the :term:`predict_proba` output. The default threshold is defined as a posterior
probability estimate of 0.5 or a decision score of 0.0. However, this default strategy
may not be optimal for the task at hand.

This example shows how to use the
:class:`~sklearn.model_selection.TunedThresholdClassifierCV` to tune the decision
threshold, depending on a metric of interest.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# The diabetes dataset
# --------------------
#
# To illustrate the tuning of the decision threshold, we will use the diabetes dataset.
# This dataset is available on OpenML: https://www.openml.org/d/37. We use the
# :func:`~sklearn.datasets.fetch_openml` function to fetch this dataset.
from sklearn.datasets import fetch_openml

diabetes = fetch_openml(data_id=37, as_frame=True, parser="pandas")
data, target = diabetes.data, diabetes.target

# %%
# We look at the target to understand the type of problem we are dealing with.
target.value_counts()

# %%
# We can see that we are dealing with a binary classification problem. Since the
# labels are not encoded as 0 and 1, we make it explicit that we consider the class
# labeled "tested_negative" as the negative class (which is also the most frequent)
# and the class labeled "tested_positive" the positive as the positive class:
neg_label, pos_label = target.value_counts().index

# %%
# We can also observe that this binary problem is slightly imbalanced where we have
# around twice more samples from the negative class than from the positive class. When
# it comes to evaluation, we should consider this aspect to interpret the results.
#
# Our vanilla classifier
# ----------------------
#
# We define a basic predictive model composed of a scaler followed by a logistic
# regression classifier.
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

model = make_pipeline(StandardScaler(), LogisticRegression())
model

# %%
# We evaluate our model using cross-validation. We use the accuracy and the balanced
# accuracy to report the performance of our model. The balanced accuracy is a metric
# that is less sensitive to class imbalance and will allow us to put the accuracy
# score in perspective.
#
# Cross-validation allows us to study the variance of the decision threshold across
# different splits of the data. However, the dataset is rather small and it would be
# detrimental to use more than 5 folds to evaluate the dispersion. Therefore, we use
# a :class:`~sklearn.model_selection.RepeatedStratifiedKFold` where we apply several
# repetitions of 5-fold cross-validation.
import pandas as pd

from sklearn.metrics import make_scorer, recall_score
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate

scoring = {
    "accuracy": "accuracy",
    "balanced_accuracy": "balanced_accuracy",
    "recall": make_scorer(recall_score, pos_label=pos_label),
}

cv_scores = [
    "train_accuracy",
    "test_accuracy",
    "train_balanced_accuracy",
    "test_balanced_accuracy",
    "train_recall",
    "test_recall",
]

cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)


def run_model_cv_results(model, cv):
    cv_results_model = pd.DataFrame(
        cross_validate(
            model,
            data,
            target,
            scoring=scoring,
            cv=cv,
            return_train_score=True,
            return_estimator=True,
        )
    )
    return cv_results_model


cv_results_vanilla_model = run_model_cv_results(model, cv)
cv_results_vanilla_model[cv_scores].aggregate(["mean", "std"]).T


# %%
# Our predictive model succeeds to grasp the relationship between the data and the
# target. The training and testing scores are close to each other, meaning that our
# predictive model is not overfitting. We can also observe that the balanced accuracy is
# lower than the accuracy, due to the class imbalance previously mentioned.
#
# For this classifier, we let the decision threshold, used convert the probability of
# the positive class into a class prediction, to its default value: 0.5. However, this
# threshold might not be optimal. If our interest is to maximize the balanced accuracy,
# we should select another threshold that would maximize this metric.
#
# The :class:`~sklearn.model_selection.TunedThresholdClassifierCV` meta-estimator allows
# to tune the decision threshold of a classifier given a metric of interest.
#
# Tuning the decision threshold
# -----------------------------
#
# We create a :class:`~sklearn.model_selection.TunedThresholdClassifierCV` and
# configure it to maximize the balanced accuracy. We evaluate the model using the same
# cross-validation strategy as previously.
from sklearn.model_selection import TunedThresholdClassifierCV

tuned_model = TunedThresholdClassifierCV(estimator=model, scoring="balanced_accuracy")
cv_results_tuned_model = run_model_cv_results(tuned_model, cv)
cv_results_tuned_model[cv_scores].aggregate(["mean", "std"]).T

# %%
# In comparison with the vanilla model, we observe that the balanced accuracy score
# increased. Of course, it comes at the cost of a lower accuracy score. It means that
# our model is now more sensitive to the positive class but makes more mistakes on the
# negative class.
#
# However, it is important to note that this tuned predictive model is internally the
# same model as the vanilla model: they have the same fitted coefficients. We will
# see shortly that the models differ in their decision threshold.
import matplotlib.pyplot as plt

vanilla_model_coef = pd.DataFrame(
    [est[-1].coef_.ravel() for est in cv_results_vanilla_model["estimator"]],
    columns=diabetes.feature_names,
)
tuned_model_coef = pd.DataFrame(
    [est.estimator_[-1].coef_.ravel() for est in cv_results_tuned_model["estimator"]],
    columns=diabetes.feature_names,
)

fig, ax = plt.subplots(ncols=2, figsize=(12, 4), sharex=True, sharey=True)
vanilla_model_coef.boxplot(ax=ax[0])
ax[0].set_ylabel("Coefficient value")
ax[0].set_title("Vanilla model")
tuned_model_coef.boxplot(ax=ax[1])
ax[1].set_title("Tuned model")
_ = fig.suptitle("Coefficients of the predictive models")

# %%
# To illustrate the flexibility in the choice of metrics, we make another classifier
# with the :class:`~sklearn.model_selection.TunedThresholdClassifierCV` class.
# Earlier we saw that maximizing the balanced accuracy came at the cost of decreased
# performance on the negative class. Suppose we want to maximize balanced accuracy,
# but still maintain some performance on the negative class, say we want our true
# negative rate (TNR) to be at least 75%. This can be done by creating a custom
# scoring function, as shown below.
import numpy as np

from sklearn.metrics import confusion_matrix, make_scorer


def balanced_accuracy_min_tnr(y_true, y_pred, bound):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    tpr = tp / (tp + fn)
    tnr = tn / (tn + fp)

    if tnr < bound:
        return -np.inf
    return (tpr + tnr) / 2  # Balanced accuracy


constrained_balanced_accuracy = make_scorer(balanced_accuracy_min_tnr, bound=0.75)

constrained_metric_model = TunedThresholdClassifierCV(
    estimator=model, scoring=constrained_balanced_accuracy
)

cv_results_constrained_metric_model = run_model_cv_results(constrained_metric_model, cv)

cv_results_constrained_metric_model[cv_scores].aggregate(["mean", "std"]).T
# %%
# As before, we observe an increase in balanced accuracy and a decrease in our accuracy
# recall scores. However, as a result of insisting that the TNR to be at least 0.75
# (by giving a score of `-np.inf` to any model that has a TNR below that bound) we
# do not see as much of an improvement in recall as our previous model. Later on we
# will plot the corresponding point of this classifier on an ROC curve to better
# understand this constraint. For now we observe the decision thresholds selected in
# our :class:`~sklearn.model_selection.TunedThresholdClassifierCV` models.
fig, ax = plt.subplots(figsize=(6, 5))


def process_decision_thresholds(cv_results):
    decision_threshold = pd.Series(
        [est.best_threshold_ for est in cv_results["estimator"]]
    )

    decision_threshold.plot.kde(ax=ax)
    ax.axvline(decision_threshold.mean(), color="k", linestyle="--")
    return decision_threshold.mean()


mean_threshold = process_decision_thresholds(cv_results_tuned_model)
mean_threshold_constrained = process_decision_thresholds(
    cv_results_constrained_metric_model
)

_ = ax.legend(
    [
        "Tuned Model",
        f"Mean Decision Threshold : {mean_threshold:.2f}",
        "Tuned Model (Constrained Metric)",
        f"Mean Decision Threshold : {mean_threshold_constrained:.2f}",
    ]
)

# %%
# On average, a decision threshold around 0.32 maximizes the balanced accuracy,
# and a decision threshold of 0.43 maximizes the constrained balanced accuracy.
# These are different from the default decision threshold of 0.5. Thus tuning the
# decision threshold is particularly important when the output of the predictive
# model is used to make decisions. Besides, the metric used to tune the decision
# threshold should be chosen carefully. Here, we used the balanced accuracy but it
# might not be the most appropriate metric for the problem at hand. The choice of the
# "right" metric is usually problem-dependent and might require some domain knowledge.
# Refer to the example entitled,
# :ref:`sphx_glr_auto_examples_model_selection_plot_cost_sensitive_learning.py`,
# for more details.


# We conclude this section by plotting the mean ROC curve and the points on the curve
# corresponding to our classifiers. To do the latter, we will compute the FPR and
# TPR on each test set passed to our classifiers. TPR is the same as recall,
# which we already have. FPR can be obtained from the TPR and balanced accuracy,
# recalling that `:math:`FPR=1 - TNR` and
# :math:`\text{Balanced Accuracy}\,=\frac{TPR+TNR}{2}`
def compute_point_roc(cv_model_results):
    ba = cv_model_results["test_balanced_accuracy"].mean()
    tpr = cv_model_results["test_recall"].mean()
    return (1 - (2 * ba - tpr), tpr)


roc_points = {
    "vanilla": compute_point_roc(cv_results_vanilla_model),
    "max_balanced_accuracy": compute_point_roc(cv_results_tuned_model),
    "constrained_balanced_accuracy": compute_point_roc(
        cv_results_constrained_metric_model
    ),
}
# %%
# We now compute the test ROC on the test set from each train, test, splited generated
# from our :class:`~sklearn.model_selection.RepeatedStratifiedKFold` instance. We saw
# earlier the underlying model of all classifiers  are the same, so it suffices to
# compute each ROC from the vanilla model. We then average these to obtain the
# mean ROC.

import warnings

from sklearn.metrics import RocCurveDisplay

test_indices = [test_index for _, test_index in cv.split(data, target)]

fig, ax = plt.subplots(figsize=(6, 6))
tprs = []

mean_fpr = np.linspace(0, 1, 100)

for i in range(len(test_indices)):
    test_idx = test_indices[i]
    estimator = cv_results_vanilla_model["estimator"][i]

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="No artists with labels found to put in legend",
            category=UserWarning,
        )

        viz = RocCurveDisplay.from_estimator(
            estimator,
            data.iloc[test_idx],
            target.iloc[test_idx],
            alpha=0.2,
            pos_label=pos_label,
            lw=1,
            ax=ax,
            plot_chance_level=(i == len(test_indices) - 1),
            label=None,
        )

        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0

ax.plot(mean_fpr, mean_tpr, color="b", label=r"Mean ROC", lw=2, alpha=0.8)

ax.plot(
    roc_points["vanilla"][0],
    roc_points["vanilla"][1],
    marker="o",
    markersize=10,
    color="orange",
    label="Cut-off point at probability 0.5",
)

ax.plot(
    roc_points["max_balanced_accuracy"][0],
    roc_points["max_balanced_accuracy"][1],
    marker=">",
    markersize=10,
    color="blue",
    label=f"Cut-off point at probability {mean_threshold:.2f}",
)

ax.plot(
    roc_points["constrained_balanced_accuracy"][0],
    roc_points["constrained_balanced_accuracy"][1],
    marker="<",
    markersize=10,
    color="green",
    label=f"Cut-off point at probability {mean_threshold_constrained:.2f}",
)

ax.legend(loc="lower right")
_ = ax.set(
    xlabel="False Positive Rate", ylabel="True Positive Rate", title="ROC Curves"
)
# %%
# We note that the point corresponding to the classifier maximizing the balanced
# accuracy is furthest 'up the curve'. Of all three points, we can see this one has
# the highest TNR (recall that :math:`FPR=1 - TNR`). Also, notice that the point
# corresponding to the classifier with the constrained metric is between the points
# corresponding to the vanilla classifier and the classifier maximizing balanced
# accuracy with no constraints. Indeed, the classifier with the constrained metric
# is discouraged from thresholds achieving a false positive rate higher than
# 1 - 0.75 = 0.25, since the TNR is constrained to be at least 0.75.
