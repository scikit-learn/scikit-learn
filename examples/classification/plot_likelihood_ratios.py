"""
=======================
Class Likelihood Ratios
=======================

This example shows a worked-out demonstration of the
:func:`class_likelihood_ratios` function, which computes the positive and
negative likelihood ratios and outputs them as a tuple (`LR+`, `LR-`).

In the first section we discuss the link between pre-test and post-test odds
given by the :ref:`class_likelihood_ratios`. We thereafter evaluate the behavior
of the likelihood ratios in some controlled scenarios. In the last section we
plot the likelihood ratios as a function of the prevalence of the positive
class.

A typical application of such functions in medicine is to assess the predictive
power of a binary diagnostic test. In such case, the pre-test probability of an
individual having the target condition can be chosen to be the prevalence, i.e.
the proportion of a particular population found to be affected by a medical
condition. The post-test probabilities represent then the probability that the
condition is truly present given a positive test result.

"""

# Authors:  Arturo Amor <david-arturo.amor-quiroz@inria.fr>
#           Olivier Grisel <olivier.grisel@ensta.org>
# %%
# Pre-test vs. post-test analysis
# ===============================
#
# Suppose we have a population of subjects with physiological measurements `X`
# that can hopefully serve as indirect bio-markers of the disease and actual
# disease indicators `y` (ground truth). Most of the people in the population do
# not carry the disease but a minority (in this case around 10%) does:

from sklearn.datasets import make_classification

X, y = make_classification(n_samples=10000, weights=[0.9], random_state=0)
print(f"Percentage of people carrying the disease: {100*y.mean():.2f}%")

# %%
# A machine learning model is built to diagnose if a person with some given
# physiological measurements is likely to carry the disease of interest. To
# evaluate the model, we need to assess its performance on a held-out test set:

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# %%
# Then we can fit our diagnosis model and compute the positive likelihood
# ratio to evaluate the usefulness of this classifier as a disease diagnosis
# tool:

from sklearn.metrics import class_likelihood_ratios
from sklearn.linear_model import LogisticRegression

clf = LogisticRegression().fit(X_train, y_train)
y_pred = clf.predict(X_test)
pos_LR, neg_LR = class_likelihood_ratios(y_test, y_pred)
print(f"LR+: {pos_LR:.3f}")

# %%
# Since the positive class likelihood ratio is much larger than 1.0, it means
# that the machine learning-based diagnosis tool is useful: the post-test odds
# are more than 12 times larger than the pre-test odds. We now choose the
# pre-test probability to be the prevalence of the disease in the held-out
# testing set.

pretest_odds = y_test.mean() / (1 - y_test.mean())
print(f"Observed pre-test odds: {pretest_odds:.3f}")
print(f"Estimated post-test odds using LR+: {pretest_odds * pos_LR:.3f}")

# %%
# The post-test probability is the probability that the condition is truly
# present given a positive test result, i.e. the number of true positives
# divided by the total number of samples. In real life applications this is
# unknown. Transforming post-test probabilities to odds, we can verify that the
# positive likelihood ratio relates pre-test and post-test odds.

posttest_prob = y_test[y_pred == 1].mean()
posttest_odds = posttest_prob / (1 - posttest_prob)

print(f"Observed post-test odds: {posttest_odds:.3f}")
# %%
# The estimated and the observed post-test odds do coincide.
#
# Cross-validation of likelihood ratios
# =====================================
#
# We assess the variability of the measurements for the class likelihood ratios
# in some particular cases. For that purpose we define a cross validating
# function that preserves the proportion of samples for both classes in each
# fold:

from sklearn.model_selection import StratifiedKFold
import pandas as pd


def cross_validate_LRs(clf, X, y, verbose=True):

    pos_LRs = []
    neg_LRs = []

    for train_index, test_index in StratifiedKFold(n_splits=10).split(X, y):

        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        pos_LR, neg_LR = class_likelihood_ratios(y_test, y_pred, raise_warning=False)
        pos_LRs.append(pos_LR)
        neg_LRs.append(neg_LR)

    class_LRs = pd.DataFrame({"LR+": pos_LRs, "LR-": neg_LRs})

    if verbose:
        print(
            "The mean cross-validated LR+ is: "
            f"{class_LRs['LR+'].mean():.2f} +/- "
            f"{class_LRs['LR+'].std():.2f}"
        )
        print(
            "The mean cross-validated LR- is: "
            f"{class_LRs['LR-'].mean():.2f} +/- "
            f"{class_LRs['LR-'].std():.2f}"
        )

    return class_LRs


# %%
# We first validate the `LogisticRegression` model with default hyperparameters
# as used in the previous section.

_ = cross_validate_LRs(LogisticRegression(), X, y)

# %%
# We confirm that the model is useful: the post-test odds are between 12 and 20
# times larger than the pre-test odds.
#
# On the contrary, let's consider a dummy model that will output random
# predictions with similar odds as the average disease prevalence in the
# training set:

from sklearn.dummy import DummyClassifier

_ = cross_validate_LRs(DummyClassifier(strategy="stratified", random_state=1234), X, y)

# %%
# Here both class likelihood ratios are compatible with 1.0 which makes this
# classifier useless as a diagnostic tool to improve disease detection.
#
# Another option for the dummy model is to always predict the most frequent
# class, which in this case is "no-disease".

_ = cross_validate_LRs(DummyClassifier(strategy="most_frequent"), X, y)

# %%
# The absence of positive predictions means there will be no true positives nor
# false positives, leading to an undefined `LR+` that by no means should be
# interpreted as an infinite `LR+` (the classifier perfectly identifying
# positive cases). In such situation the :func:`class_likelihood_ratios` function
# returns `nan` and raises a warning by default. Indeed, the value of `LR-` helps
# us discard this model.
#
# A similar scenario may arise when cross-validating highly imbalanced data with
# few samples: some folds will have no samples with the disease and therefore
# they will output no true positives nor false negatives when used for testing.
# Mathematically this leads to an infinite `LR+`, which should also not be
# interpreted as the model perfectly identifying positive cases. Such event
# leads to a higher variance of the estimated likelihood ratios, but can still
# be interpreted as an increment of the post-test odds of having the condition.

X, y = make_classification(n_samples=300, weights=[0.9], random_state=0)
_ = cross_validate_LRs(LogisticRegression(), X, y)

# %%
# Plot class likelihood ratios
# ============================
#
# The likelihood ratios are independent of the disease prevalence and can be
# extrapolated between populations regardless of any possible class imbalance,
# as long as the same model is applied to all of them. In the case where
# models/diagnostic tests are updated, for instance, to evolve along with the
# spread of a disease, one has to account for the prevalence seen during
# training before making a comparison analysis.
#
# Here we use the `make_classification` function to explore both scenarios:
# constant and evolving models. This function creates 2 clusters per class as
# shown in the plots below. We vary the prevalence of the hypothetical condition
# by passing different `weights` to the `make_classification` function, which
# ensures the data-generating process is always the same. The label `1`
# corresponds to the positive class "disease", whereas the label `0` stands for
# "no-disease". The decision boundaries are trained using a `LogisticRegression`
# classifier similarly to the previous sections.
#
# Re-trained models
# -----------------
# We first evaluate the case where the model is re-trained in each iteration
# with varying prevalence. Notice that the decision boundary changes accordingly.

import numpy as np
import matplotlib.pyplot as plt
from sklearn.inspection import DecisionBoundaryDisplay

pos_LRs = []
neg_LRs = []
pos_LRs_std = []
neg_LRs_std = []
prevalence = []
weights = np.linspace(0.1, 0.9, 6)
weights = weights[::-1]

plt.figure(figsize=(15, 12))
plt.subplots_adjust(hspace=0.25)


for n, weight in enumerate(weights):

    X, y = make_classification(
        n_samples=10_000,
        n_features=2,
        n_informative=2,
        n_redundant=0,
        weights=[weight],
        random_state=0,
    )

    clf = LogisticRegression()
    ax = plt.subplot(3, 2, n + 1)
    disp = DecisionBoundaryDisplay.from_estimator(
        clf.fit(X, y),
        X,
        response_method="predict",
        alpha=0.5,
        ax=ax,
    )
    scatter = disp.ax_.scatter(X[:, 0], X[:, 1], c=y, edgecolor="k")
    disp.ax_.set_title(f"prevalence = {y.mean():.2f}")
    disp.ax_.legend(*scatter.legend_elements())

    class_LRs = cross_validate_LRs(clf, X, y, verbose=False)

    pos_LR = class_LRs["LR+"].mean()
    neg_LR = class_LRs["LR-"].mean()
    pos_LR_std = class_LRs["LR+"].std()
    neg_LR_std = class_LRs["LR-"].std()

    pos_LRs.append(pos_LR)
    neg_LRs.append(neg_LR)
    pos_LRs_std.append(pos_LR_std)
    neg_LRs_std.append(neg_LR_std)
    prevalence.append(y.mean())

class_LRs = pd.DataFrame(
    {"LR+": pos_LRs, "LR-": neg_LRs, "LR+_std": pos_LRs_std, "LR-_std": neg_LRs_std}
)

# %%
# The following plots show that likelihood ratios vary with the prevalence when
# the diagnostic tools are re-calibrated with different prevalences even when the
# disease is truly dichotomous and all of the physiological measurements `X` are
# informative.

plt.figure(figsize=(15, 6))
ax1 = plt.subplot(1, 2, 1)
ax1.plot(prevalence, class_LRs["LR+"], color="r")
ax1.fill_between(
    prevalence,
    class_LRs["LR+"] - class_LRs["LR+_std"],
    class_LRs["LR+"] + class_LRs["LR+_std"],
    color="r",
    alpha=0.3,
)
ax1.set(
    title="Positive likelihood ratio",
    xlabel="prevalence",
    ylabel="LR+",
)

ax2 = plt.subplot(1, 2, 2)
ax2.plot(prevalence, class_LRs["LR-"], color="b")
ax2.fill_between(
    prevalence,
    class_LRs["LR-"] + class_LRs["LR-_std"],
    class_LRs["LR-"] - class_LRs["LR-_std"],
    color="b",
    alpha=0.3,
)
ax2.set(
    title="Negative likelihood ratio",
    xlabel="prevalence",
    ylabel="LR-",
)

plt.show()

# %%
# Both estimated class likelihood ratios decrease with increasing prevalence:
#
# For `LR+` this means that the post-test odds of the condition truly being
# present given a positive test result become closer to the pre-test odds as the
# proportion of the population affected by the condition increases. In the limit
# of prevalence=1.0 all the population is affected and the act of testing is no
# more informative than assuming the sample is a priori positive.
#
# For `LR-` it means that the probability of a sample of the positive class
# being misclassified as belonging to the negative class drops to zero when the
# fraction of population with the condition increases, which should be the case
# if the classifier is truly informative.
#
# Constant model
# --------------
# Notice that in this case, the model is trained with a prevalence of 50% and
# the resulting decision boundary is used to evaluate said model over
# populations with varying prevalence.

# %%
pos_LRs = []
neg_LRs = []
pos_LRs_std = []
neg_LRs_std = []
prevalence = []

weights = np.linspace(0.1, 0.8, 6)
weights = weights[::-1]

X, y = make_classification(
    n_samples=10_000,
    n_features=2,
    n_informative=2,
    n_redundant=0,
    weights=[0.5],
    random_state=0,
)

clf = LogisticRegression()
class_LRs = cross_validate_LRs(clf, X, y, verbose=False)

pos_LR_base = class_LRs["LR+"].mean()
neg_LR_base = class_LRs["LR-"].mean()
pos_LR_base_std = class_LRs["LR+"].std()
neg_LR_base_std = class_LRs["LR-"].std()

plt.figure(figsize=(15, 12))
plt.subplots_adjust(hspace=0.25)

for n, weight in enumerate(weights):

    X_test, y_test = make_classification(
        n_samples=10_000,
        n_features=2,
        n_informative=2,
        n_redundant=0,
        weights=[weight],
        random_state=0,
    )

    ax = plt.subplot(3, 2, n + 1)
    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X_test,
        response_method="predict",
        alpha=0.5,
        ax=ax,
    )
    scatter = disp.ax_.scatter(X_test[:, 0], X_test[:, 1], c=y_test, edgecolor="k")
    disp.ax_.set_title(f"prevalence = {y_test.mean():.2f}")
    disp.ax_.legend(*scatter.legend_elements())

    y_pred = clf.predict(X_test)
    pos_LR, neg_LR = class_likelihood_ratios(y_test, y_pred, raise_warning=False)
    pos_LRs.append(pos_LR)
    neg_LRs.append(neg_LR)
    prevalence.append(y_test.mean())

class_LRs = pd.DataFrame({"LR+": pos_LRs, "LR-": neg_LRs})
# %%
# In the plots below we can observe that the class likelihood ratios remain
# within one standard deviation of those computed with the constant model.

plt.figure(figsize=(15, 6))

ax1 = plt.subplot(1, 2, 1)
ax1.plot(prevalence, class_LRs["LR+"], "r+", label="extrapolation through populations")

ax1.axhline(y=pos_LR_base + pos_LR_base_std, color="r", linestyle="--")
ax1.axhline(
    y=pos_LR_base - pos_LR_base_std,
    color="r",
    linestyle="--",
    label="trained model confidence band",
)
ax1.set(
    title="Positive likelihood ratio",
    xlabel="prevalence",
    ylabel="LR+",
    ylim=[0, 5],
)
plt.legend(loc="lower right")

ax2 = plt.subplot(1, 2, 2)
ax2.plot(prevalence, class_LRs["LR-"], "b+", label="extrapolation through populations")
ax2.axhline(y=neg_LR_base + neg_LR_base_std, color="b", linestyle="--")
ax2.axhline(
    y=neg_LR_base - neg_LR_base_std,
    color="b",
    linestyle="--",
    label="trained model confidence band",
)
ax2.set(
    title="Negative likelihood ratio",
    xlabel="prevalence",
    ylabel="LR-",
    ylim=[0, 0.5],
)
plt.legend(loc="lower right")

plt.show()
# %%
