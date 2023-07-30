"""
===============================================================
Post-tuning decision threshold based on cost-sensitive learning
===============================================================

Once a classifier is trained, the output of the :term:`predict` method output class
label predictions corresponding to a thresholding of either the :term:`decision
function` or the :term:`predict_proba` output. For a binary classifier, the default
threshold is defined as a posterior probability estimate of 0.5 or a decision score of
0.0.

However, this default strategy may not be optimal for the task at hand.
Here, we use the "Statlog" German credit dataset [1]_ to illustrate a use case.
In this dataset, the task is to predict whether a person has a "good" or "bad" credit.
In addition, a cost-matrix is provided that specifies the cost of
misclassification. Specifically, misclassifying a "bad" credit as "good" is five
times more costly than misclassifying a "good" credit as "bad".

We use the :class:`~sklearn.model_selection.TunedThresholdClassifier` to select the
cut-off point of the decision function that minimizes the provided business
cost.

.. topic:: References

    .. [1] "Statlog (German Credit Data) Data Set", UCI Machine Learning Repository,
       `Link
       <https://archive.ics.uci.edu/ml/datasets/Statlog+%28German+Credit+Data%29>`_.

    .. [2] `Charles Elkan, "The Foundations of Cost-Sensitive Learning",
       International joint conference on artificial intelligence.
       Vol. 17. No. 1. Lawrence Erlbaum Associates Ltd, 2001.
       <https://cseweb.ucsd.edu//~elkan/rescale.pdf>`_
"""

# %%
# Cost-sensitive learning with constant gains and costs
# -----------------------------------------------------
#
# In this first section, we illustrate the use of the
# :class:`~sklearn.model_selection.TunedThresholdClassifier` in a setting of
# cost-sensitive learning when the gains and costs associated to each entry of the
# confusion matrix are constant. We use the problematic presented in [2]_ using the
# "Statlog" German credit dataset [1]_.
#
# "Statlog" German credit dataset
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
# Many features are categorical and usually string-encoded. We need to encode
# these categories when we develop our predictive model. Let's check the targets.
y.value_counts()

# %%
# Another observation is that the dataset is imbalanced. We would need to be careful
# when evaluating our predictive model and use a family of metrics that are adapted
# to this setting.
#
# In addition, we observe that the target is string-encoded. Some metrics
# (e.g. precision and recall) require to provide the label of interest also called
# the "positive label". Here, we define that our goal is to predict whether or not
# a sample is a "bad" credit.
pos_label, neg_label = "bad", "good"

# %%
# To carry our analysis, we split our dataset using a single stratified split.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=0)

# %%
# We are ready to design our predictive model and the associated evaluation strategy.
#
# Evaluation metrics
# ^^^^^^^^^^^^^^^^^^
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
from sklearn.metrics import confusion_matrix


def fpr_score(y, y_pred, neg_label, pos_label):
    cm = confusion_matrix(y, y_pred, labels=[neg_label, pos_label])
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
# the custom scorers in a dictionary. To use them, we need to pass the fitted model,
# the data and the target on which we want to evaluate the predictive model.
from sklearn.metrics import make_scorer, precision_score, recall_score

tpr_score = recall_score  # TPR and recall are the same metric
scoring = {
    "precision": make_scorer(precision_score, pos_label=pos_label),
    "recall": make_scorer(recall_score, pos_label=pos_label),
    "fpr": make_scorer(fpr_score, neg_label=neg_label, pos_label=pos_label),
    "tpr": make_scorer(tpr_score, pos_label=pos_label),
}

# %%
# In addition, the original research [1]_ defines a business metric. They provide a
# cost-matrix which encodes that predicting a "bad" credit as "good" is 5 times more
# costly than the opposite. We define a python function that weight the confusion
# matrix and return the overall cost.
import numpy as np


def gain_cost_score(y, y_pred, neg_label, pos_label):
    cm = confusion_matrix(y, y_pred, labels=[neg_label, pos_label])
    cost_matrix = np.array([[0, -1], [-5, 0]])
    return np.sum(cm * cost_matrix)


scoring["cost_gain"] = make_scorer(
    gain_cost_score, neg_label=neg_label, pos_label=pos_label
)
# %%
# Vanilla predictive model
# ^^^^^^^^^^^^^^^^^^^^^^^^
#
# We first design our predictive model consisting of a
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier`. We encode the categorical
# features with an :class:`~sklearn.preprocessing.OrdinalEncoder` but the numerical
# features are kept as they are. To identify the categorical columns, we use the helper
# function :func:`~sklearn.compose.make_column_selector` and the fact that the
# categorical features are stored as `category` dtype.
from sklearn.compose import ColumnTransformer
from sklearn.compose import make_column_selector as selector
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
# We recall that these curves give insights on the statistical performance of the
# predictive model for different cut-off points. For the Precision-Recall curve, the
# reported metrics are the precision and recall and for the ROC curve, the reported
# metrics are the TPR (same as recall) and FPR.
#
# Here, the different cut-off points correspond to different levels of posterior
# probability estimates ranging between 0 and 1. By default, `model.predict` uses a
# cut-off point at a probability estimate of 0.5. The metrics for such cut-off point are
# reported with the blue dot on the curves: it corresponds to the statistical
# performance of the model when using `model.predict`.
#
# However, we recall that the original aim was to minimize the cost (or maximize the
# gain) by the business metric. We can compute the value of the business metric:
print(f"Business defined metric: {scoring['cost_gain'](model, X_test, y_test)}")

# %%
# At this stage we don't know if any other cut-off can lead to a greater gain. To find
# the optimal one, we need to compute the cost-gain using the business metric for all
# possible cut-off points and choose the best. This strategy can be quite tedious to
# implement by hand, but the :class:`~sklearn.model_selection.TunedThresholdClassifier`
# class is here to help us. It automatically computes the cost-gain for all possible
# cut-off points and optimizes for the `objective_metric`.
#
# .. _cost_sensitive_learning_example:
#
# Tuning the cut-off point
# ^^^^^^^^^^^^^^^^^^^^^^^^
#
# We use :class:`~sklearn.model_selection.TunedThresholdClassifier` to tune the cut-off
# point. We need to provide the business metric to optimize as well as the
# positive label. Internally, the optimum cut-off point is chosen such that it
# maximizes the business metric via cross-validation. By default a 5-fold
# stratified cross-validation is used.
from sklearn.model_selection import TunedThresholdClassifier

tuned_model = TunedThresholdClassifier(
    estimator=model,
    pos_label=pos_label,
    objective_metric=scoring["cost_gain"],
)
tuned_model.fit(X_train, y_train)

# %%
# We plot the ROC and Precision-Recall curves for the vanilla model and the tuned model.
# Also we plot the cut-off points that would be used by each model.
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21, 6))

linestyles = ("dashed", "dotted")
markerstyles = ("o", ">")
colors = ("tab:blue", "tab:orange")
names = ("Vanilla GBDT", "Tuned GBDT")
for idx, (est, linestyle, marker, color, name) in enumerate(
    zip((model, tuned_model), linestyles, markerstyles, colors, names)
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
    tuned_model.decision_thresholds_, tuned_model.objective_scores_, color="tab:orange"
)
axs[2].plot(
    tuned_model.decision_threshold_,
    tuned_model.objective_score_,
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
# on the same training data. In a later section, we discuss more in detail the
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
print(f"Business defined metric: {scoring['cost_gain'](tuned_model, X_test, y_test)}")

# %%
# We observe that the decision generalized on the testing set leading to a better
# business score.
#
# .. _tunedthresholdclassifier_no_cv:
#
# Consideration regarding model refitting and cross-validation
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# In the above experiment, we use the default setting of the
# :class:`~sklearn.model_selection.TunedThresholdClassifier`. In particular, the cut-off
# point is tuned using a 5-fold stratified cross-validation. Also, the
# underlying predictive model is refitted on the entire training data once the
# cut-off point is chosen.
#
# These two strategies can be changed by providing the `refit` and `cv` parameters.
# For instance, one could provide a fitted `estimator` and set `cv="prefit"`, in which
# case the cut-off point is found on the entire dataset provided at fitting time.
# Also, the underlying classifier is not be refitted. Here, we can try to do such
# experiment.
model.fit(X_train, y_train)
tuned_model.set_params(cv="prefit").fit(X_train, y_train)


# %%
# Then, we evaluate our model with the same approach as before:
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21, 6))

linestyles = ("dashed", "dotted")
markerstyles = ("o", ">")
colors = ("tab:blue", "tab:orange")
names = ("Vanilla GBDT", "Tuned GBDT")
for idx, (est, linestyle, marker, color, name) in enumerate(
    zip((model, tuned_model), linestyles, markerstyles, colors, names)
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
    tuned_model.decision_thresholds_, tuned_model.objective_scores_, color="tab:orange"
)
axs[2].plot(
    tuned_model.decision_threshold_,
    tuned_model.objective_score_,
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
# :class:`~sklearn.model_selection.TunedThresholdClassifier` is not the same as the data
# used to train the underlying classifier. This could happen sometimes when the idea is
# just to tune the predictive model on a completely new validation set without a costly
# complete refit.
#
# In the case that cross-validation is too costly, a potential alternative is to use a
# single train-test split by providing a floating number in range `[0, 1]` to the `cv`
# parameter. It splits the data into a training and testing set. Let's explore this
# option:
tuned_model.set_params(cv=0.75).fit(X_train, y_train)

# %%
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21, 6))

linestyles = ("dashed", "dotted")
markerstyles = ("o", ">")
colors = ("tab:blue", "tab:orange")
names = ("Vanilla GBDT", "Tuned GBDT")
for idx, (est, linestyle, marker, color, name) in enumerate(
    zip((model, tuned_model), linestyles, markerstyles, colors, names)
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
    tuned_model.decision_thresholds_, tuned_model.objective_scores_, color="tab:orange"
)
axs[2].plot(
    tuned_model.decision_threshold_,
    tuned_model.objective_score_,
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
# for the variability of the fit/predict process and thus we are unable to know if there
# is any variance in the cut-off point. The repeated cross-validation averages out
# this effect.
#
# Another observation concerns the ROC and Precision-Recall curves of the tuned model.
# As expected, these curves differ from those of the vanilla model, given that we
# trained the underlying classifier on a subset of the data provided during fitting and
# reserved a validation set for tuning the cut-off point.
#
# Cost-sensitive learning when gains and costs are not constant
# -------------------------------------------------------------
#
# As stated in [2]_, gains and costs are generally not constant in real-world problems.
# In this section, we use a similar example as in [2]_ by using credit cards
# records.
#
# The credit card dataset
# ^^^^^^^^^^^^^^^^^^^^^^^
credit_card = fetch_openml(data_id=1597, as_frame=True, parser="pandas")
credit_card.frame.info()

# %%
# The dataset contains information about credit card records from which some are
# fraudulent and others are legitimate. The goal is therefore to predict whether or
# not a credit card record is fraudulent.
#
# In addition, we have extra information regarding the amount of each card transaction.
# This information is used to define the business metric later.
columns_to_drop = ["Class", "Amount"]
data = credit_card.frame.drop(columns=columns_to_drop)
target = credit_card.frame["Class"].astype(int)
amount = credit_card.frame["Amount"].to_numpy()

# %%
# First, we check the class distribution of the datasets.
target.value_counts(normalize=True)

# %%
# The dataset is highly imbalanced with fraudulent transaction representing only 0.17%
# of the data. Additionally, we check the distribution of the amount of the fraudulent
# transactions.
fraud = target == 1
amount_fraud = amount[fraud]
_, ax = plt.subplots()
ax.hist(amount_fraud, bins=100)
ax.set_title("Amount of fraud transaction")
_ = ax.set_xlabel("Amount ($)")

# %%
# Now, we create the business metric that depends on the amount of each transaction. We
# define the cost matrix similarly to [2]_. Accepting a legitimate transaction provides
# a gain of 2% of the amount of the transaction. However, accepting a fraudulent
# transaction result in a loss of the amount of the transaction. As stated in [2]_, the
# gain and loss related to refusals (of fraudulent and legitimate transactions) are not
# trivial to define. Here, we define that a refusal of a legitimate transaction is
# estimated to a loss of $5 while the refusal of a fraudulent transaction is estimated
# to a gain of $50 dollars and the amount of the transaction. Therefore, we define the
# following function to compute the total benefit of a given decision:


def business_metric(y_true, y_pred, amount):
    mask_true_positive = (y_true == 1) & (y_pred == 1)
    mask_true_negative = (y_true == 0) & (y_pred == 0)
    mask_false_positive = (y_true == 0) & (y_pred == 1)
    mask_false_negative = (y_true == 1) & (y_pred == 0)
    fraudulent_refuse = (mask_true_positive.sum() * 50) + amount[
        mask_true_positive
    ].sum()
    fraudulent_accept = -amount[mask_false_negative].sum()
    legitimate_refuse = mask_false_positive.sum() * -5
    legitimate_accept = (amount[mask_true_negative] * 0.02).sum()
    return fraudulent_refuse + fraudulent_accept + legitimate_refuse + legitimate_accept


# %%
# From this business metric, we create a scikit-learn scorer that given a fitted
# classifier and a test set compute the business metric. In this regard, we use
# the :func:`~sklearn.metrics.make_scorer` factory. The variable `amount` is an
# additional metadata to be passed to the scorer and we need to use
# :ref:`metadata routing <metadata_routing>` to take into account this information.
sklearn.set_config(enable_metadata_routing=True)
business_scorer = make_scorer(business_metric).set_score_request(amount=True)

# %%
# We first start to train a dummy classifier to have some baseline results.
from sklearn.model_selection import train_test_split

data_train, data_test, target_train, target_test, amount_train, amount_test = (
    train_test_split(
        data, target, amount, stratify=target, test_size=0.5, random_state=42
    )
)

# %%
from sklearn.dummy import DummyClassifier

easy_going_classifier = DummyClassifier(strategy="constant", constant=0)
easy_going_classifier.fit(data_train, target_train)
benefit_cost = business_scorer(
    easy_going_classifier, data_test, target_test, amount=amount_test
)
print(f"Benefit/cost of our easy-going classifier: ${benefit_cost:,.2f}")

# %%
# A classifier that predict all transactions as legitimate would create a profit of
# around $220,000. We make the same evaluation for a classifier that predicts all
# transactions as fraudulent.
intolerant_classifier = DummyClassifier(strategy="constant", constant=1)
intolerant_classifier.fit(data_train, target_train)
benefit_cost = business_scorer(
    intolerant_classifier, data_test, target_test, amount=amount_test
)
print(f"Benefit/cost of our intolerant classifier: ${benefit_cost:,.2f}")

# %%
# Such a classifier create a loss of around $670,000. A predictive model should allow
# us to make a profit larger than $220,000. It is interesting to compare this business
# metric with another "standard" statistical metric such as the balanced accuracy.
from sklearn.metrics import get_scorer

balanced_accuracy_scorer = get_scorer("balanced_accuracy")
print(
    "Balanced accuracy of our easy-going classifier: "
    f"{balanced_accuracy_scorer(easy_going_classifier, data_test, target_test):.3f}"
)
print(
    "Balanced accuracy of our intolerant classifier: "
    f"{balanced_accuracy_scorer(intolerant_classifier, data_test, target_test):.3f}"
)

# %%
# This is not a surprise that the balanced accuracy is at 0.5 for both classifiers.
# However, we need to be careful in the rest of the evaluation: we potentially can
# obtain a model with a decent balanced accuracy but that does not make any profit.
# In this case, the model would be useless for our business.
#
# Let's now create a predictive model using a logistic regression without tuning the
# decision threshold.
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

model = make_pipeline(StandardScaler(), LogisticRegression(random_state=42)).fit(
    data_train, target_train
)

print(
    "Benefit/cost of our logistic regression: "
    f"${business_scorer(model, data_test, target_test, amount=amount_test):,.2f}"
)
print(
    "Balanced accuracy of our logistic regression: "
    f"{balanced_accuracy_scorer(model, data_test, target_test):.3f}"
)

# %%
# By observing the balanced accuracy, we see that our predictive model is learning
# some associations between the features and the target. The business metric also shows
# that our model is beating the baseline in terms of profit and it would be already
# beneficial to use it instead of ignoring the fraud detection problem.
#
# Now the question is: is our model optimum for the type of decision that we want to do?
# Up to now, we did not optimize the decision threshold. We use the
# :class:`~sklearn.model_selection.TunedThresholdClassifier` to optimize the decision
# given our business scorer.
tuned_model = TunedThresholdClassifier(
    estimator=model,
    objective_metric=business_scorer,
    n_thresholds=100,
    n_jobs=2,
)

# %%
# Since our business scorer requires the amount of each transaction, we need to pass
# this information in the `fit` method. The
# :class:`~sklearn.model_selection.TunedThresholdClassifier` is in charge of
# automatically dispatching this metadata to the underlying scorer.
tuned_model.fit(data_train, target_train, amount=amount_train)

# %%
print(
    "Benefit/cost of our logistic regression: "
    f"${business_scorer(tuned_model, data_test, target_test, amount=amount_test):,.2f}"
)
print(
    "Balanced accuracy of our logistic regression: "
    f"{balanced_accuracy_scorer(tuned_model, data_test, target_test):.3f}"
)

# %%
# We observe that tuning the decision threshold increases the profit of our model.
# Eventually, the balanced accuracy also increased. Note that it might not always be
# the case because the statistical metric is not necessarily a surrogate of the
# business metric. It is therefore important, whenever possible, optimize the decision
# threshold with respect to the business metric.
