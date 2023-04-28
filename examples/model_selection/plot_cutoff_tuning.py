"""
==================================================
Post-tuning the cut-off point of decision function
==================================================

Once a classifier is trained, the output of the :term:`predict` method output hard
predictions corresponding to a thresholding of either the :term:`decision function`
or the :term:`predict_proba` output. For binary classifier, the default threshold is
defined as a probability score of 0.5 or a decision function value of 0.0.

However, it happens that this default strategy is not optimized for the task at hand.
Here, we use the "Statlog" German credit dataset [1]_ to illustrate such an application.
In this dataset, the task is to predict whether a person is "good" or "bad" credit.
In addition, cost-matrix is provided where it is defined that classifying a "bad"
credit as "good" is 5 times more costly than classifying a "good" credit as "bad".

We use the :class:`~sklearn.model_selection.CutOffClassifier` to tune the
cut-off point of the decision function to minimize the business cost provided to us.

.. topic:: References
   .. [1] "Statlog (German Credit Data) Data Set", UCI Machine Learning Repository,
        https://archive.ics.uci.edu/ml/datasets/Statlog+%28German+Credit+Data%29
"""

# %%
# "Statlog" German credit dataset
# -------------------------------
#
# We fetch the German credit dataset from OpenML.
import sklearn
from sklearn.datasets import fetch_openml

sklearn.set_config(transform_output="pandas")
german_credit = fetch_openml(data_id=31, as_frame=True, parser="pandas")
X, y = german_credit.data, german_credit.target

# %%
# We check the feature types available in `X`.
X.info()

# %%
X.head()

# %%
# Many features are categorical and usually string-encoded. We need to encode
# these categories when we develop ou predictive model. Let's check the targets.
y.value_counts()

# %%
# A first observation is that the dataset is imbalanced. We would need to be careful
# when evaluating our predictive model and use a family of metrics that are adapted
# to this setting.
#
# In addition, we observe that the target is string-encoded. Some metrics
# (e.g. precision and recall) require to provide the label of interest also called
# the "positive label". Here, we define that our goal is to predict whether or not
# a sample is a "bad" credit.
pos_label = "bad"

# %%
# To carry our analysis, we split our dataset using a single stratified split.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=0)

# %%
# We are ready to design our predictive model and the associated evaluation strategy.
#
# Evaluation metrics
# ------------------
#
# In this section, we define a set of metrics that we use later. To see
# the effect of tuning the cut-off point, we evaluate the predictive model using
# the Receiver Operating Characteristic (ROC) curve and the Precision-Recall curve.
# The values reported on these plots are therefore the true positive rate (TPR) and
# the false positive rate (FPR) for the ROC curve and the precision and recall for the
# Precision-Recall curve.
#
# From these four metrics, scikit-learn does not provide a scorer for the FPR. We
# therefore need to define a small custom function to compute it.
import numpy as np
from sklearn.metrics import confusion_matrix, make_scorer, precision_score, recall_score


def fpr_score(y, y_pred, **kwargs):
    cm = confusion_matrix(y, y_pred)
    classes = np.unique(y)
    pos_label = kwargs.get("pos_label", classes[-1])
    pos_label_idx = np.searchsorted(classes, pos_label)
    if pos_label_idx == 0:
        cm = cm[::-1, ::-1]
    tn, fp, _, _ = cm.ravel()
    tnr = tn / (tn + fp)
    return 1 - tnr


# %%
# As previously stated, the "positive label" is not defined as the value "1" and calling
# some of the metrics with this non-standard value raise an error. We need to
# provide the indication of the "positive label" to the metrics.
#
# We therefore need to define a scikit-learn scorer using
# :func:`~sklearn.metrics.make_scorer` where the information is passed. We store all
# created scorer in a dictionary. To use them, we need to pass the fitted model and
# the data and target on which we want to evaluate the predictive model.
tpr_score = recall_score  # TPR and recall are the same metric
scoring = {
    "precision": make_scorer(precision_score, pos_label=pos_label),
    "recall": make_scorer(recall_score, pos_label=pos_label),
    "fpr": make_scorer(fpr_score, pos_label=pos_label),
    "tpr": make_scorer(tpr_score, pos_label=pos_label),
}

# %%
# In addition, the original research [1]_ defines a business metric. They provide a
# cost-matrix where they define that predicting a "bad" credit as "good" is 5 times more
# costly than the opposite. We define a dictionary containing this information and a
# score function that computes the cost.
cost_gain_matrix = {"tp": 0, "tn": 0, "fp": -1, "fn": -5}


def gain_cost_score(y, y_pred, **kwargs):
    cm = confusion_matrix(y, y_pred)
    classes = np.unique(y)
    pos_label = kwargs.get("pos_label", classes[-1])
    pos_label_idx = np.searchsorted(classes, pos_label)
    if pos_label_idx == 0:
        cm = cm[::-1, ::-1]
    costs_and_gain = np.array(
        [
            [kwargs["cost_gain_matrix"]["tn"], kwargs["cost_gain_matrix"]["fp"]],
            [kwargs["cost_gain_matrix"]["fn"], kwargs["cost_gain_matrix"]["tp"]],
        ]
    )
    return (costs_and_gain * cm).sum()


scoring["cost_gain"] = make_scorer(
    gain_cost_score, pos_label=pos_label, cost_gain_matrix=cost_gain_matrix
)
# %%
# Vanilla predictive model
# ------------------------
#
# Design of the predictive model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# In this section, we design our predictive model. We use a gradient boosting classifier
# using :class:`~sklearn.ensemble.HistGradientBoostingClassifier`. However, we need to
# encode the categorical features with numerical values and we therefore use an
# :class:`~sklearn.preprocessing.OrdinalEncoder` to do so. The numerical features are
# kept as-is. To recover the name of the categorical columns, we use the helper function
# :func:`~sklearn.compose.make_column_selector` and the fact that the categorical
# features are stored as `category` dtype.
from sklearn.compose import ColumnTransformer, make_column_selector as selector
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OrdinalEncoder

categorical_columns = selector(dtype_include="category")(X)
numerical_columns = selector(dtype_exclude="category")(X)

preprocessor = ColumnTransformer(
    [
        (
            "categorical",
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1),
            categorical_columns,
        ),
        ("numerical", "passthrough", numerical_columns),
    ],
    verbose_feature_names_out=False,
)
model = Pipeline(
    [
        ("preprocessor", preprocessor),
        (
            "classifier",
            HistGradientBoostingClassifier(
                categorical_features=categorical_columns, random_state=0
            ),
        ),
    ]
)

model.fit(X_train, y_train)

# %%
# Evaluation of the predictive model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# We evaluate the performance of our predictive model using the ROC and Precision-Recall
# curves.
import matplotlib.pyplot as plt
from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

PrecisionRecallDisplay.from_estimator(
    model, X_test, y_test, pos_label=pos_label, ax=axs[0], name="GBDT"
)
axs[0].plot(
    scoring["recall"](model, X_test, y_test),
    scoring["precision"](model, X_test, y_test),
    marker="o",
    markersize=10,
    color="tab:blue",
    label="Default cut-off point at a probability of 0.5",
)
axs[0].set_title("Precision-Recall curve")
axs[0].legend()

RocCurveDisplay.from_estimator(
    model,
    X_test,
    y_test,
    pos_label=pos_label,
    ax=axs[1],
    name="GBDT",
    plot_chance_level=True,
)
axs[1].plot(
    scoring["fpr"](model, X_test, y_test),
    scoring["tpr"](model, X_test, y_test),
    marker="o",
    markersize=10,
    color="tab:blue",
    label="Default cut-off point at a probability of 0.5",
)
axs[1].set_title("ROC curve")
axs[1].legend()
_ = fig.suptitle("Evaluation of the vanilla GBDT model")

# %%
# We recall that these curves corresponds show the statistical performance of the
# predictive model for different cut-off points. For the Precision-Recall curve, the
# reported metrics are the precision and recall and for the ROC curve, the reported
# metrics are the TPR (same as recall) and FPR.
#
# Here, the different cut-off points corresponds at different level of probability
# scores ranging between 0 and 1. By default, `model.predict` uses a cut-off point at
# a probability of 0.5. Thus, the metrics for this cut-off point is reported with the
# blue dot on the curves: it corresponds to the statistical performance of the model
# when using `model.predict`.
#
# However, we recall that the original aim was to minimize the cost (or maximize the
# gain) by the business metric. We can compute the value of the business metric:
scoring["cost_gain"](model, X_test, y_test)

# %%
# At this stage, we don't know if some other cut-off points could lead to a greater
# gain. To be able to decide, if another cut-off point is better, we would need to
# compute the cost-gain using the business metric for all possible cut-off points and
# choose the optimal one. This strategy is quite tedious to implement and the
# :class:`~sklearn.metrics.CutOffClassifier` class is here to help us. It automatically
# compute the cost-gain for all possible cut-off points and choose the optimal one.
#
# Tuning the cut-off point
# ------------------------
#
# We use :class:`~sklearn.metrics.CutOffClassifier` to tune the cut-off point. We need
# to provide the business metric to optimize as well as the positive label. Internally,
# the optimum cut-off point is chosen such that it maximized the business metric
# via cross-validation. By default a 5-fold stratified cross-validation is used.
from sklearn.model_selection import CutOffClassifier

model_tuned = CutOffClassifier(
    estimator=model,
    pos_label=pos_label,
    objective_metric=cost_gain_matrix,
)
model_tuned.fit(X_train, y_train)

# %%
# We plot the ROC and Precision-Recall curves for the vanilla model and the tuned model.
# Also we plot the cut-off points that would be used by each model.
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21, 6))

linestyles = ("dashed", "dotted")
markerstyles = ("o", ">")
colors = ("tab:blue", "tab:orange")
names = ("Vanilla GBDT", "Tuned GBDT")
for idx, (est, linestyle, marker, color, name) in enumerate(
    zip((model, model_tuned), linestyles, markerstyles, colors, names)
):
    decision_threshold = getattr(est, "decision_threshold_", 0.5)
    PrecisionRecallDisplay.from_estimator(
        est,
        X_test,
        y_test,
        pos_label=pos_label,
        linestyle=linestyle,
        color=color,
        ax=axs[0],
        name=name,
    )
    axs[0].plot(
        scoring["recall"](est, X_test, y_test),
        scoring["precision"](est, X_test, y_test),
        marker,
        markersize=10,
        color=color,
        label=f"Cut-off point at probability of {decision_threshold:.2f}",
    )
    RocCurveDisplay.from_estimator(
        est,
        X_test,
        y_test,
        pos_label=pos_label,
        linestyle=linestyle,
        color=color,
        ax=axs[1],
        name=name,
        plot_chance_level=idx == 1,
    )
    axs[1].plot(
        scoring["fpr"](est, X_test, y_test),
        scoring["tpr"](est, X_test, y_test),
        marker,
        markersize=10,
        color=color,
        label=f"Cut-off point at probability of {decision_threshold:.2f}",
    )

axs[0].set_title("Precision-Recall curve")
axs[0].legend()
axs[1].set_title("ROC curve")
axs[1].legend()

axs[2].plot(
    model_tuned.decision_thresholds_, model_tuned.objective_scores_, color="tab:orange"
)
axs[2].plot(
    model_tuned.decision_threshold_,
    model_tuned.objective_score_,
    "o",
    markersize=10,
    color="tab:orange",
    label="Optimal cut-off point for the business metric",
)
axs[2].legend()
axs[2].set_xlabel("Decision threshold (probability)")
axs[2].set_ylabel("Objective score (using cost-matrix)")
axs[2].set_title("Objective score as a function of the decision threshold")

_ = fig.suptitle("Comparison of the cut-off point for the vanilla and tuned GBDT model")

# %%
# The first remark is that both classifiers have exactly the same ROC and
# Precision-Recall curves. It is expected because by default, the classifier is fitted
# on the same training data. In a later section, we discuss more in details the
# available options regarding model refitting and cross-validation.
#
# The second remark is that the cut-off points of the vanilla and tuned model are
# different. To understand why the tuned model has chosen this cut-off point, we can
# look at the right-hand side plot that plots the objective score that is our exactly
# the same as our business metric. We see that the optimum threshold corresponds to the
# maximum of the objective score.
#
# We can now check if choosing this cut-off point leads to a better score on the testing
# set:
scoring["cost_gain"](model_tuned, X_test, y_test)

# %%
# We observe that the decision generalized on the testing set leading to a better
# business score.
#
# Consideration regarding model refitting and cross-validation
# ------------------------------------------------------------
#
# In the above experiment, we use the default parameter of the
# :class:`~sklearn.model_selection.CutOffClassifier`. By default, the cut-off point is
# tuned using a 5-fold stratified cross-validation. Also, the underlying predictive
# model is refitted on the entire training data once the cut-off point is chosen.
#
# These two strategies can be changed by providing the `refit` and `cv` parameters.
# For instance, one could provide a fitted `estimator` and set `cv="prefit"`, in which
# case the cut-off point is found on the entire dataset provided at fitting time.
# Also, the underlying classifier is not be refitted. Here, we can try to do such
# experiment.
model.fit(X_train, y_train)
model_tuned.set_params(cv="prefit").fit(X_train, y_train)


# %%
# Then, we evaluate our model with the same approach as before:
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21, 6))

linestyles = ("dashed", "dotted")
markerstyles = ("o", ">")
colors = ("tab:blue", "tab:orange")
names = ("Vanilla GBDT", "Tuned GBDT")
for idx, (est, linestyle, marker, color, name) in enumerate(
    zip((model, model_tuned), linestyles, markerstyles, colors, names)
):
    decision_threshold = getattr(est, "decision_threshold_", 0.5)
    PrecisionRecallDisplay.from_estimator(
        est,
        X_test,
        y_test,
        pos_label=pos_label,
        linestyle=linestyle,
        color=color,
        ax=axs[0],
        name=name,
    )
    axs[0].plot(
        scoring["recall"](est, X_test, y_test),
        scoring["precision"](est, X_test, y_test),
        marker,
        markersize=10,
        color=color,
        label=f"Cut-off point at probability of {decision_threshold:.2f}",
    )
    RocCurveDisplay.from_estimator(
        est,
        X_test,
        y_test,
        pos_label=pos_label,
        linestyle=linestyle,
        color=color,
        ax=axs[1],
        name=name,
        plot_chance_level=idx == 1,
    )
    axs[1].plot(
        scoring["fpr"](est, X_test, y_test),
        scoring["tpr"](est, X_test, y_test),
        marker,
        markersize=10,
        color=color,
        label=f"Cut-off point at probability of {decision_threshold:.2f}",
    )

axs[0].set_title("Precision-Recall curve")
axs[0].legend()
axs[1].set_title("ROC curve")
axs[1].legend()

axs[2].plot(
    model_tuned.decision_thresholds_, model_tuned.objective_scores_, color="tab:orange"
)
axs[2].plot(
    model_tuned.decision_threshold_,
    model_tuned.objective_score_,
    "o",
    markersize=10,
    color="tab:orange",
    label="Optimal cut-off point for the business metric",
)
axs[2].legend()
axs[2].set_xlabel("Decision threshold (probability)")
axs[2].set_ylabel("Objective score (using cost-matrix)")
axs[2].set_title("Objective score as a function of the decision threshold")

_ = fig.suptitle("Tuned GBDT model without refitting and using the entire dataset")

# %%
# We observe the that the optimum cut-off point is different than in the previous
# experiment. If we look at the right-hand side plot, we observe that the objective
# score has large plateau with a minimum cost (around 0). This behavior is symptomatic
# of an overfitting. Because we disable cross-validation, we tuned the cut-off point on
# the same set as the model was trained on, and this is the reason for the observed
# overfitting.
#
# This option should therefore be used with caution. One needs to make sure that the
# data providing at fitting time to the
# :class:`~sklearn.model_selection.CutOffClassifier` is not the same as the data used to
# train the underlying classifier. This could happen sometimes when the idea is just
# to tune the predictive model on a completely new validation set without a costly
# complete refit.
#
# In the case that cross-validation is too costly, a potential alternative is to use
# a single train-test split by providing a floating number in range `[0, 1]` to the
# `cv` parameter. It splits the data into a training and testing set. Let's
# explore this option:
model_tuned.set_params(cv=0.75).fit(X_train, y_train)

# %%
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21, 6))

linestyles = ("dashed", "dotted")
markerstyles = ("o", ">")
colors = ("tab:blue", "tab:orange")
names = ("Vanilla GBDT", "Tuned GBDT")
for idx, (est, linestyle, marker, color, name) in enumerate(
    zip((model, model_tuned), linestyles, markerstyles, colors, names)
):
    decision_threshold = getattr(est, "decision_threshold_", 0.5)
    PrecisionRecallDisplay.from_estimator(
        est,
        X_test,
        y_test,
        pos_label=pos_label,
        linestyle=linestyle,
        color=color,
        ax=axs[0],
        name=name,
    )
    axs[0].plot(
        scoring["recall"](est, X_test, y_test),
        scoring["precision"](est, X_test, y_test),
        marker,
        markersize=10,
        color=color,
        label=f"Cut-off point at probability of {decision_threshold:.2f}",
    )
    RocCurveDisplay.from_estimator(
        est,
        X_test,
        y_test,
        pos_label=pos_label,
        linestyle=linestyle,
        color=color,
        ax=axs[1],
        name=name,
        plot_chance_level=idx == 1,
    )
    axs[1].plot(
        scoring["fpr"](est, X_test, y_test),
        scoring["tpr"](est, X_test, y_test),
        marker,
        markersize=10,
        color=color,
        label=f"Cut-off point at probability of {decision_threshold:.2f}",
    )

axs[0].set_title("Precision-Recall curve")
axs[0].legend()
axs[1].set_title("ROC curve")
axs[1].legend()

axs[2].plot(
    model_tuned.decision_thresholds_, model_tuned.objective_scores_, color="tab:orange"
)
axs[2].plot(
    model_tuned.decision_threshold_,
    model_tuned.objective_score_,
    "o",
    markersize=10,
    color="tab:orange",
    label="Optimal cut-off point for the business metric",
)
axs[2].legend()
axs[2].set_xlabel("Decision threshold (probability)")
axs[2].set_ylabel("Objective score (using cost-matrix)")
axs[2].set_title("Objective score as a function of the decision threshold")

_ = fig.suptitle("Tuned GBDT model without refitting and using the entire dataset")

# %%
# Regarding the cut-off point, we observe that the optimum is similar to the multiple
# repeated cross-validation case. However, be aware that a single split does not account
# for the variability of the fit/predict process and thus we are enable to know if there
# is any variance in the cut-off point. The repeated cross-validation averages this
# effect.
#
# The second observation is about the ROC and Precision-Recall curves the tuned model.
# We observe that the curve is now different from the vanilla model. It is not
# surprising because we train the underlying classifier on a subset of the data provided
# at `fit` while keeping a validation set to tune the cut-off point.
