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
# We evaluate our model using cross-validation. We use the precision
# and recall to report the performance of our model.

# Cross-validation allows us to study the variance of the decision threshold across
# different splits of the data. However, the dataset is rather small and it would be
# detrimental to use more than 5 folds to evaluate the dispersion. Therefore, we use
# a :class:`~sklearn.model_selection.RepeatedStratifiedKFold` where we apply several
# repetitions of 5-fold cross-validation.
import pandas as pd

from sklearn.metrics import make_scorer, precision_score, recall_score
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate

scoring = {
    "recall": make_scorer(recall_score, pos_label=pos_label),
    "precision": make_scorer(precision_score, pos_label=pos_label),
}

cv_scores = [
    "train_recall",
    "test_recall",
    "train_precision",
    "test_precision",
]

# TODO change later to 5, 10 repeats
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=42)

params = {
    "scoring": scoring,
    "cv": cv,
    "return_train_score": True,
    "return_estimator": True,
}

cv_results_vanilla_model = pd.DataFrame(cross_validate(model, data, target, **params))
cv_results_vanilla_model[cv_scores].aggregate(["mean", "std"]).T
# %%
# Our predictive model succeeds to grasp the relationship between the data and the
# target. The training and testing scores are close to each other, meaning that our
# predictive model is not overfitting.
#
# For this classifier, we let the decision threshold, used convert the probability of
# the positive class into a class prediction, remain at its default value: 0.5.
# However, this threshold might not be optimal. For instance, if we are interested in
# maximizing the f1-score, which is defined as the harmonic mean between
# precision and recall, we probably need a different decision threshold.
#
# The :class:`~sklearn.model_selection.TunedThresholdClassifierCV` meta-estimator allows
# to tune the decision threshold of a classifier given a metric of interest.
#
# Tuning the decision threshold
# -----------------------------
#
# We create a :class:`~sklearn.model_selection.TunedThresholdClassifierCV` and
# configure it to maximize the f1-score, using the
# `func:f1_score` function from the :mod:`sklearn.metrics` module.
# We evaluate the model using the same cross-validation strategy as previously.
from sklearn.metrics import f1_score
from sklearn.model_selection import TunedThresholdClassifierCV

tuned_model = TunedThresholdClassifierCV(
    estimator=model, scoring=make_scorer(f1_score, pos_label=pos_label)
)
cv_results_tuned_model = pd.DataFrame(
    cross_validate(tuned_model, data, target, **params)
)

cv_results_tuned_model[cv_scores].aggregate(["mean", "std"]).T

# %%
# It is important to note that this tuned predictive model is internally the
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
# Only the decision threshold of each model was changed during the cross-validation.
fig, ax = plt.subplots(figsize=(6, 5))

decision_threshold = pd.Series(
    [est.best_threshold_ for est in cv_results_tuned_model["estimator"]]
)

decision_threshold.plot.kde(ax=ax)
ax.axvline(
    decision_threshold.mean(),
    color="k",
    linestyle="--",
    label=f"Mean decision threshold: {decision_threshold.mean():.2f}",
)

ax.legend(loc="lower right")
ax.set_xlabel("Decision threshold")
_ = ax.set_title(
    "Distribution of the decision threshold \nacross different cross-validation folds"
)
# %%
# On average, a decision threshold around 0.30 maximizes the f1-score,
# which is lower than the default decision threshold of 0.5. To understand the
# tradeoffs better, we view the corresponding points of these thresholds
# on a precision-recall curve. For ease of presentation, we will consider only
# one of the test sets in our cross-validation split, so that only one precision-recall
# curve is considered. We will also be defining a function to generate this plot,
# as later on we will want to plot the precision-recall curve again with an
# additional point, corresponding to another decision threshold.

from sklearn.metrics import PrecisionRecallDisplay

cv_models = {
    "Vanilla": cv_results_vanilla_model,
    "Max F1": cv_results_tuned_model,
}


def plot_pr_curve():

    fig, ax = plt.subplots(figsize=(6, 5))

    test_indices = [test_index for _, test_index in cv.split(data, target)]
    idx = 0  # Consider only the first test set

    display = PrecisionRecallDisplay.from_estimator(
        cv_results_vanilla_model["estimator"][idx],
        data.iloc[test_indices[idx]],
        target.iloc[test_indices[idx]],
        name="precision-recall",
        plot_chance_level=1,
        ax=ax,
    )

    display.ax_.set_title("Precision-Recall curve")
    symbols = ["o", "<", ">"]
    colors = ["orange", "blue", "green"]

    for i, model in enumerate(cv_models):
        ax.plot(
            cv_models[model]["test_recall"][idx],
            cv_models[model]["test_precision"][0],
            marker=symbols[i],
            markersize=10,
            color=colors[i],
            label=f"{model} model",
        )
    ax.legend()


plot_pr_curve()

# %%
# From the plot above, we see that the new model achieves
# higher recall and lower precision when compared
# to the vanilla model.
# This is intuitive when we think about the corresponding decision thresholds.
# By decreasing the threshold from 0.5
# to 0.3, we allow our classifier to label more samples as positive,
# which also increases the number of false positives.
#
# The metric used to tune the decision threshold should be chosen carefully.
# The choice of the "right" metric to optimize for is usually problem-dependent,
# and might require some domain knowledge. Here, we used the f1-score, which
# usually achieves a good balance between precision and recall. However, this
# may not be most appropriate for our situation. In classifying patients as diabetic,
# we may decide that recall is more important than precision, since the cost of
# leaving diabetes untreated can be very harfmul to the patient. Having said that,
# we may still want to maintain a baseline level of precision, say 60%.
#
# This can be done as before, with the
# :class:`~sklearn.model_selection.TunedThresholdClassifierCV` class,
# and a custom scoring function.

import numpy as np

from sklearn.metrics import confusion_matrix


def constrained_recall(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(
        y_true, y_pred, labels=[neg_label, pos_label]
    ).ravel()
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)

    if precision < 0.6:
        return -np.inf
    return recall


constrained_recall_model = TunedThresholdClassifierCV(
    estimator=model, scoring=make_scorer(constrained_recall)
)

cv_results_constrained_recall_model = pd.DataFrame(
    cross_validate(constrained_recall_model, data, target, **params)
)

cv_results_constrained_recall_model[cv_scores].aggregate(["mean", "std"]).T
# %%
# We conclude this section by plotting the precision-recall curve a second time,
# but with an additional point that corresponds to our new model which maximizes
# recall under a constraint of maintaining 60% precision.
cv_models["Max Recall Constrained"] = cv_results_constrained_recall_model
plot_pr_curve()
# %%
# We can see now the tradeoffs between precision and recall across all our models,
# on the given test set. The precision of our new model is lower than the
# vanilla model, but higher than the model maximizing the f1-score.
# On the other hand, our new model has higher recall than the vanilla model, but
# lower recall than the f1-score model. By combining custom scoring functions and the
# :class:`~sklearn.model_selection.TunedThresholdClassifierCV` class, the user gains
# flexibility when deciding what metric to optimize for. The decision threshold can be
# carefully modified by choosing a metric and maximizing it.
