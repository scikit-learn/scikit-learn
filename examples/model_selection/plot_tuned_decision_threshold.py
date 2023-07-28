"""
==================================================
Post-tuning the cut-off point of decision function
==================================================

Once a classifier is trained, the output of the :term:`predict` method output class
label predictions corresponding to a thresholding of either the :term:`decision
function` or the :term:`predict_proba` output. For a binary classifier, the default
threshold is defined as a posterior probability estimate of 0.5 or a decision score of
0.0. However, this default strategy may not be optimal for the task at hand.

This example shows how to use the
:class:`~sklearn.model_selection.TunedThresholdClassifier` to tune the decision function
threshold, depending on a metric of interest as well as under a specific constraint.
"""

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
# We see that we are dealing with a binary classification problem. Since the labels are
# not encoded as 0 and 1, we will store which label we considered the negative class
# and which one we considered the positive class: "tested_negative" will be considered
# the negative class and "tested_positive" the positive class.
#
# We also observed that this binary problem is slightly imbalanced where we have around
# twice more samples from the negative class than from the positive class. When it
# comes to evaluation, we should consider this aspect to interpret the results.
neg_label, pos_label = target.value_counts().index

# %%
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
import pandas as pd

from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate

scoring = ["accuracy", "balanced_accuracy"]
cv_scores = [
    "train_accuracy",
    "test_accuracy",
    "train_balanced_accuracy",
    "test_balanced_accuracy",
]
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
cv_results_vanilla_model = pd.DataFrame(
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
cv_results_vanilla_model[cv_scores].aggregate(["mean", "std"]).T

# %%
# Our predictive model succeed to grasp relationship between the data and the target.
# The training and testing scores are close to each other, meaning that our predictive
# model is not overfitting. We also observe that the balanced accuracy is lower than
# the accuracy, due to the class imbalanced previously mentioned.
#
# For this classifier, we used a decision threshold of 0.5 to convert the probability
# of the positive class into a class prediction. However, this threshold might not be
# optimal. If our interest is to maximize the balanced accuracy, we should select
# another threshold that would maximize this metric.
#
# The :class:`~sklearn.model_selection.TunedThresholdClassifier` allows to tune the
# decision threshold of a classifier given a metric of interest.
#
# Tuning the decision threshold
# -----------------------------
#
# We create a :class:`~sklearn.model_selection.TunedThresholdClassifier` and we
# configure it to maximize the balanced accuracy. We evaluate the model using the same
# cross-validation strategy as previously.
from sklearn.model_selection import TunedThresholdClassifier

tuned_model = TunedThresholdClassifier(
    estimator=model, objective_metric="balanced_accuracy"
)
cv_results_tuned_model = pd.DataFrame(
    cross_validate(
        tuned_model,
        data,
        target,
        scoring=scoring,
        cv=cv,
        return_train_score=True,
        return_estimator=True,
    )
)
cv_results_tuned_model[cv_scores].aggregate(["mean", "std"]).T

# %%
# In comparison with the vanilla model, we observe that the balanced accuracy score
# increased. Of course, it comes at the cost of a lower accuracy score. It means that
# our model is now more sensitive to the positive class but makes more mistakes on the
# negative class.
#
# However, it is important to note that this tuned predictive model is internally the
# same model as the vanilla model.
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
decision_threshold = pd.Series(
    [est.decision_threshold_ for est in cv_results_tuned_model["estimator"]],
)
ax = decision_threshold.plot.kde()
ax.axvline(
    decision_threshold.mean(),
    color="k",
    linestyle="--",
    label=f"Mean decision threshold: {decision_threshold.mean():.2f}",
)
ax.set_xlabel("Decision threshold")
ax.legend(loc="upper right")
_ = ax.set_title(
    "Distribution of the decision threshold \nacross different cross-validation folds"
)

# %%
# In average, a decision threshold around 0.32 is maximizing the balanced accuracy. It
# is thus different from the default decision threshold of 0.5. Tuning the decision
# threshold is thus particularly important when the output of the predictive model
# is used to make decisions. Besides, the metric used to tune the decision threshold
# should be chosen carefully. Here, we used the balanced accuracy but it might not be
# the most appropriate metric for the problem at hand. The choice of the "right" metric
# is usually problem-dependent and might require some domain knowledge. Refer to the
# example entitled,
# :ref:`sphx_glr_auto_examples_model_selection_plot_tuned_threshold_classifier.py`,
# for more details.
#
# Tuning the decision threshold under constraint
# ----------------------------------------------
#
# In some cases, we do not want only to maximize a specific metric but instead maximize
# a metric while satisfying a constraint on another metric. In the current example, we
# could imagine that the decision of our predictive model will be reviewed by a medical
# doctor. In this case, this doctor will only accept a ratio of false positive.
# Therefore, we are interesting at maximizing the true positive rate while having a
# a false positive rate lower than a given threshold.
#
# The :class:`~sklearn.model_selection.TunedThresholdClassifier` allows to tune the
# decision threshold with such specification. We should how to proceed using a single
# cross-validation split to display the Receiver Operating Characteristic (ROC) curves
# to get intuition on the problem.
#
# First, we split the data into a training and testing set.

# %%
from sklearn.model_selection import train_test_split

data_train, data_test, target_train, target_test = train_test_split(
    data, target, random_state=42
)

# %%
# Now, we will train both the vanilla and tuned model on the training set. We recall
# that the tuned model is internally maximizing the balanced accuracy for the moment.
model.fit(data_train, target_train)
tuned_model.fit(data_train, target_train)

# %%
# To show the benefit on optimizing a metric under constraint, we will evaluate the
# models using the ROC curves statistics: the true positive rate (TPR) and the false
# positive rate (FPR).
#
# The FPR is not defined in scikit-learn and we define it below:
from sklearn.metrics import confusion_matrix, make_scorer, recall_score


def fpr_score(y, y_pred, neg_label, pos_label):
    cm = confusion_matrix(y, y_pred, labels=[neg_label, pos_label])
    tn, fp, _, _ = cm.ravel()
    tnr = tn / (tn + fp)
    return 1 - tnr


tpr_score = recall_score  # TPR and recall are the same metric
scoring = {
    "fpr": make_scorer(fpr_score, neg_label=neg_label, pos_label=pos_label),
    "tpr": make_scorer(tpr_score, pos_label=pos_label),
}

# %%
# Now, we plot the ROC curve of both models and the FPR and TPR statistics for the
# decision thresholds of both models.
from sklearn.metrics import RocCurveDisplay

disp = RocCurveDisplay.from_estimator(
    model, data_test, target_test, name="Vanilla model", linestyle="--", alpha=0.5
)
RocCurveDisplay.from_estimator(
    tuned_model,
    data_test,
    target_test,
    name="Tuned model",
    linestyle="-.",
    alpha=0.5,
    ax=disp.ax_,
)
disp.ax_.plot(
    scoring["fpr"](model, data_test, target_test),
    scoring["tpr"](model, data_test, target_test),
    marker="o",
    markersize=10,
    color="tab:blue",
    label="Default cut-off point at a probability of 0.5",
)
disp.ax_.plot(
    scoring["fpr"](tuned_model, data_test, target_test),
    scoring["tpr"](tuned_model, data_test, target_test),
    marker=">",
    markersize=10,
    color="tab:orange",
    label=f"Cut-off point at probability of {tuned_model.decision_threshold_:.2f}",
)
disp.ax_.legend()
_ = disp.ax_.set_title("ROC curves")

# %%
# We observe that both models have the same ROC curves. This is expected since the tuned
# model is only a post-processing step of the vanilla model. The tuning is only
# changing the decision threshold threshold as displayed by the markers blue and orange.
# To optimize the balanced accuracy, the tuned model moved the decision threshold is
# moved from 0.5 to 0.22. By shifting this point, we increase the FPR while increasing
# the TPR: in short we make more false positive but also more true positive. This is
# exactly what we concluded in the previous section when looking at the balanced
# accuracy score.
#
# However, this decision threshold might not be acceptable for our medical doctor. He
# might be instead interested to have a low FPR, let say lower than 5%. For this level
# of FPR, he would like our predictive model to maximize the TPR.
#
# The :class:`~sklearn.model_selection.TunedThresholdClassifier` allows to specify such
# constraint by providing the name of the metric and the constraint value. Here, we use
# `max_tpr_at_tnr_constraint` which is exactly what we want. Since the true negative
# rate (TNR) is equal to 1 - FPR, we can rewrite the constraint value as
# `1 - 0.05 = 0.95`.

# %%
constraint_value = 0.95
tuned_model.set_params(
    objective_metric="max_tpr_at_tnr_constraint",
    constraint_value=constraint_value,
    pos_label=pos_label,
)
tuned_model.fit(data_train, target_train)

# %%
# Now, we can plot the ROC curves and analyse the results.
import matplotlib.pyplot as plt

_, axs = plt.subplots(ncols=2, figsize=(12, 5))

disp = RocCurveDisplay(
    fpr=1 - tuned_model.objective_scores_[0],
    tpr=tuned_model.objective_scores_[1],
    estimator_name="ROC of the tuned model",
    pos_label=pos_label,
)
axs[0].plot(
    1 - tuned_model.objective_score_[0],
    tuned_model.objective_score_[1],
    marker="o",
    markersize=10,
    color="tab:blue",
    label=f"Cut-off point at probability of {tuned_model.decision_threshold_:.2f}",
)
axs[0].axvline(
    1 - constraint_value, 0, 1, color="tab:blue", linestyle="--", label="FPR constraint"
)
axs[0].set_title("Average ROC curve for the tuned model\nacross CV folds")
RocCurveDisplay.from_estimator(
    model,
    data_test,
    target_test,
    name="Vanilla model",
    linestyle="--",
    alpha=0.5,
    ax=axs[1],
)
RocCurveDisplay.from_estimator(
    tuned_model,
    data_test,
    target_test,
    name="Tuned model",
    linestyle="-.",
    alpha=0.5,
    ax=axs[1],
)
axs[1].plot(
    scoring["fpr"](model, data_test, target_test),
    scoring["tpr"](model, data_test, target_test),
    marker="o",
    markersize=10,
    color="tab:blue",
    label="Default cut-off point at a probability of 0.5",
)
axs[1].plot(
    1 - tuned_model.objective_score_[0],
    tuned_model.objective_score_[1],
    marker="^",
    markersize=10,
    color="tab:orange",
    label=f"Cut-off point at probability of {tuned_model.decision_threshold_:.2f}",
)
axs[1].legend()
axs[1].set_title("ROC curves")
disp.plot(ax=axs[0])

# %%
# We start with the right-hand side plot. It depicts the ROC curves as in the previous
# section. We observe that the control point of the tuned model moved to a low FPR
# that was defined by our constraint. To achieve this low FPR, the decision threshold
# was moved to a probability of 0.72.
#
# The left-hand side plot shows the averaged ROC curve on the internal validation set
# across the different cross-validation folds. This curve is used to define the decision
# threshold. The vertical dashed line represents the FPR constraint that we defined.
# The decision threshold corresponds to the maximum TPR on the left of this dashed line
# and is represented by a blue marker.
#
# An important point to note is that the decision threshold is defined on averaged
# statistics on an internal validation set. It means that the constraint is respected
# on the train/validation dataset but not necessarily on the test set, in case the
# statistical performance of the model differ from the train/validation set to the test
# set (i.e. overfitting).
