"""
==================================================
Statistical comparison of models using grid search
==================================================

This example illustrates how to statistically compare the performance of models
trained and evaluated using :class:`sklearn.model_selection.GridSearchCV`.

"""

print(__doc__)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from itertools import combinations
from math import factorial
from scipy.stats import t
from sklearn.datasets import make_moons
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold
from sklearn.svm import SVC

# set seaborn parameters
sns.set_style("darkgrid")

###############################################################################
# We will start by simulating moon shaped data (where the ideal separation
# between classes is non-linear), adding to it a moderate degree of noise:

X, y = make_moons(noise=0.352, random_state=1)

# plot simulated data
sns.scatterplot(
    X[:, 0], X[:, 1], hue=y,
    marker='o', s=25, edgecolor='k', legend=False
    ).set_title("Data")


###############################################################################
# We will compare the performance of SVC estimators that vary on their `kernel`
# parameter, to decide which choice of this hyper-parameter predicts our
# simulated data best. We will evaluate the performance of the models using
# :class:`sklearn.model_selection.RepeatedStratifiedKFold`, repeating 10 times
# a 10-fold stratified cross validation using a different randomization of the
# data in each repetition.

param_grid = [
    {'kernel': ['linear']},
    {'kernel': ['poly'], 'degree': [2, 3]},
    {'kernel': ['rbf']}
    ]

svc = SVC(random_state=0)

cv = RepeatedStratifiedKFold(
    n_splits=10, n_repeats=10, random_state=0
    )

search = GridSearchCV(
    estimator=svc, param_grid=param_grid,
    scoring='roc_auc', cv=cv
    )
search.fit(X, y)

###############################################################################
# We can now inspect the results of our search:

results_df = pd.DataFrame(search.cv_results_)
results_df[
    ['params', 'rank_test_score', 'mean_test_score', 'std_test_score']
    ]

###############################################################################
# We can see that the estimator using the `rbf` kernel performed better,
# closely followed by `linear`. Both estimators with a `poly` kernel performed
# worse, with the one using a two degree polynomial achieving a much lower
# perfomance than all other models.
#
# The output of GridSearchCV does not provide information on the certainty we
# have on the differences between the models. To evaluate this, we need to
# conduct a statistical test. Specifically, to contrast the performance of two
# models we should statistically compare their vectors of AUC scores (where
# each observation is the performance obtained per fold and run).
#
# However, the score vectors of the models are not independent: we iteratively
# used the same partitions of the data to evaluate them, increasing the
# correlation between the performance of the models. This means that some
# partitions of the data can make the distinction of the classes particularly
# easy or hard to find for all models, and thus their scores will co-vary.
#
# Let's inspect this partition effect by plotting the performance of all models
# in each fold, and calculating the correlation between models across folds:

# create df of model scores ordered by perfomance
model_scores = results_df.rename(
    index={0: 'linear', 1: 'poly2', 2: 'poly3', 3: 'rbf'}
    )
model_scores = (
    model_scores.sort_values(
        by=['rank_test_score']
        ).filter(like='split')
    )

# plot 30 examples of dependency between cv fold and AUC scores
fig, ax = plt.subplots()
sns.lineplot(
    data=model_scores.transpose().iloc[:30],
    dashes='', palette='Set1', marker='o', alpha=.5, ax=ax
    )
ax.set_xlabel("CV fold", size=12, labelpad=10)
ax.set_ylabel("Model AUC", size=12)
ax.tick_params(bottom=True, labelbottom=False)

# print correlation of AUC scores across folds
print(f"Correlation of models:\n {model_scores.transpose().corr()}")

###############################################################################
# We can observe that the performance of the models highly depends on the fold.
#
# As a consequence we will be underestimating the variance computed in our
# statistical tests, increasing the number of false positive errors (i.e.
# detecting a significant difference between models when such does no exist)
# [1]_.
#
# Several approaches have been developed to correct for the variance in these
# cases. In this example we will explore one of these possibilities implemented
# using two different statistical approaches: one Frequentist and one Bayesian.


###############################################################################
# Comparing two models: Frequentist approach
# ------------------------------------------
#
# We can start by asking: "Is the first model significantly better than the
# second model (when ranked by `mean_test_score`)?"
#
# If we wanted to answer this question using a Frequentist approach we could
# run a paired t-test. Many variants of the latter have been developed to
# account for the partition effect described in the previous section. We will
# use the one proven to obtain the highest replicability scores while
# mantaining a low rate of false postitives and false negatives: the Nadeau and
# Bengio's corrected t-test [2]_ that uses a 10 times repeated 10-fold cross
# validation [3]_.
#
# This corrected paired t-test is computed as:
#
# :math:`t=\frac{\frac{1}{k.r}\sum_{i=1}^{k}\sum_{j=1}^{r}x_{ij}}
# {\sqrt{(\frac{1}{k.r}+\frac{n_2}{n_1})\hat{\sigma}^2}}`
#
# where :math:`k` is the number of folds and :math:`r` the number of
# repetitions in the cross-validation, and :math:`n_2` is the number of
# observations used for testing while :math:`n_1` is the the number of
# observations used for training.
#
# Let's implement a corrected one tailed paired t-test to compare the first and
# sencond model, and let's compute the corresponding p-value:


def correct_std(differences, n, n_train, n_test):
    """
    Calculates the standard deviation of a set of observations implementing
    Nadeau and Bengio's variance correction
    """
    corrected_var = np.var(differences, ddof=1) * ((1/n) + (n_test/n_train))
    corrected_std = np.sqrt(corrected_var)
    return corrected_std


def compute_corrected_ttest(differences, n, df, n_train, n_test):
    """
    Computes right-tailed paired t-test with corrected variance
    """
    mean = np.mean(differences)
    corrected_std = correct_std(differences, n, n_train, n_test)
    t_stat = mean / corrected_std
    p_val = t.sf(np.abs(t_stat), df)  # right-tailed t-test
    return t_stat, p_val


model_1_scores = model_scores.iloc[0].values
model_2_scores = model_scores.iloc[1].values

differences = model_1_scores - model_2_scores

n = differences.shape[0]
df = n - 1
n_train = len(list(cv.split(X, y))[0][0])
n_test = len(list(cv.split(X, y))[0][1])

t_stat, p_val = compute_corrected_ttest(differences, n, df, n_train, n_test)
print(f"Corrected t-value: {np.round(t_stat, 3)}\n"
      f"Corrected p-value: {np.round(p_val, 3)}")


###############################################################################
# We can compare the corrected t- and p-values with the uncorrected ones:

t_stat_uncorrected = (
    np.mean(differences) / np.sqrt(np.var(differences, ddof=1) / n)
    )
p_val_uncorrected = t.sf(np.abs(t_stat_uncorrected), df)

print(f"Uncorrected t-value: {np.round(t_stat_uncorrected, 3)}\n"
      f"Uncorrected p-value: {np.round(p_val_uncorrected, 3)}")


###############################################################################
# Using the conventional significance alpha level at .05, we observe that the
# uncorrected t-test concludes that the first model is significantly better
# than the second.
#
# The corrected approach, in contrast, fails to detect this difference.
#
# In the latter case, however, the Frequentist approach does not let us
# conclude that the first and second model have an equivalent performance. If
# we wanted to make this assertion we need to use a Bayesian approach.


###############################################################################
# Comparing two models: Bayesian approach
# ---------------------------------------
# We can use Bayesian estimation to calculate the probability that the first
# model is better than the second. Bayesian estimation will output a
# distribution that specifies the probability of each parameter value, in this
# example being the mean of the differences in the performance of two models.
#
# To obtain the posterior distribution we need to define a prior that models
# our beliefs of how those means are distributed before looking at the data,
# and multiply it by a likelihood function that computes how likely are our
# observations given the values that the mean of differences could take.
#
# Bayesian estimation can be carried out in many forms to answer our question,
# but in this example we will implement the approach suggested by Benavoli and
# collegues [4]_.
#
# One way of defining our posterior using a `closed-form expression
# <https://en.wikipedia.org/wiki/Closed-form_expression>`_ is to select a
# prior `conjugate <https://en.wikipedia.org/wiki/Conjugate_prior>`_ to the
# likelihood function. Benavoli and collegues [4]_ show that when comparing the
# performance of two classifiers we can model the prior as a Normal-Gamma
# distribution (with both mean and variance unknown)
# `conjugate <https://en.wikipedia.org/wiki/Conjugate_prior>`_ to a normal
# likelihood, to thus express the posterior as a normal distribution.
# Marginalizing out the variance from this normal posterior, we can define the
# posterior of the mean parameter as a Student T distribution. Specifically:
#
# :math:`St(\mu;n-1,\overline{x},(\frac{1}{n}+\frac{n_2}{n_1})\hat{\sigma}^2)`
#
# where :math:`n` is the total number of observations.
#
# Notice that we are using Nadeau and Bengio's corrected variance in our
# Bayesian approach as well.
#
# Let's compute the posterior:

# intitialize random variable
t_post = t(
    df, loc=np.mean(differences),
    scale=correct_std(differences, n, n_train, n_test)
    )

###############################################################################
# Let's plot the posterior distribution:

x = np.linspace(t_post.ppf(0.001), t_post.ppf(0.999), 100)

plt.plot(x, t_post.pdf(x))
plt.xticks(np.arange(-0.04, 0.06, 0.01))
plt.fill_between(x, t_post.pdf(x), 0, facecolor='blue', alpha=.2)
plt.ylabel("Probability density")
plt.xlabel("Mean difference")
plt.title("Posterior distribution")


###############################################################################
# We can calculate the probability that the first model is better than the
# second by computing the area under the curve of the posterior distirbution
# from zero to infinity. And also the reverse: we can calculate the probability
# that the second model is better than the first by computing the area under
# the curve from minus infinity to zero.

better_prob = t_post.sf(0)

print(f"Probability of {model_scores.index[0]} being more accurate than "
      f"{model_scores.index[1]}: {np.round(better_prob, 3)}")
print(f"Probability of {model_scores.index[1]} being more accurate than "
      f"{model_scores.index[0]}: {np.round((1-better_prob), 3)}")

###############################################################################
# In contrast with the Frequentist approach, we can compute the probability
# that one model is better than the other.
#
# Note that we obtained similar results as those in the Frequentist approach.
# Given our choice of priors we are in essence performing the same
# computations, but we are allowed to make different assertions.


###############################################################################
# Sometimes we are interested in determining the probabilities that our models
# have an equivalent performance, where "equivalent" is defined in a practical
# way. A default approach [4]_ is to define estimators as practically
# equivalent when they differ by less than 1% in their accuracy. But we could
# also define this practical equivalence taking into account the problem we are
# trying to solve. E.g. a difference of 5% in accuracy would mean an increase
# of $1000 in sales, and we consider any quantity above that as relevant for
# our business.
#
# In this example we are going to follow the suggestion in [4]_ and define the
# region of practical equivalence (ROPE) to be [-0.01, 0.01]. That is, we will
# consider two models as practically equivelent if they differ by less than 1%
# in their performance.
#
# To compute the probabilities of the classifiers being practically equivalent,
# we calculate the area under the curve of the posterior over the ROPE
# interval:

rope_interval = [-0.01, 0.01]
rope_prob = t_post.cdf(rope_interval[1]) - t_post.cdf(rope_interval[0])

print(f"Probability of {model_scores.index[0]} and {model_scores.index[1]} "
      f"being practically equivalent: {np.round(rope_prob, 3)}")


###############################################################################
# We can plot how the posterior is distributed over the ROPE interval:

x_rope = np.linspace(rope_interval[0], rope_interval[1], 100)

plt.plot(x, t_post.pdf(x))
plt.xticks(np.arange(-0.04, 0.06, 0.01))
plt.vlines([-0.01, 0.01], ymin=0, ymax=(np.max(t_post.pdf(x)) + 1))
plt.fill_between(x_rope, t_post.pdf(x_rope), 0, facecolor='blue', alpha=.2)
plt.ylabel("Probability density")
plt.xlabel("Mean difference")
plt.title("Posterior distribution under the ROPE")


###############################################################################
# As suggested in [4]_, we can further interpret these probabilities using the
# same criteria as the Frequentist approach: Is the probability of falling
# inside the ROPE bigger than 95% (alpha value of 5%)?  In that case we can
# conclude that both models are practically equivalent.


###############################################################################
# The Bayesian estimation approach also allows us to compute how uncertain we
# are about our estimation of the difference. This can be calculated using
# credible intervals which compute the range of values that have a particular
# probability of containing the true mean of differences.
#
# Let's determine the credible intervals of our data using 50%, 75% and 95%:

cred_intervals = []
intervals = [0.5, 0.75, 0.95]

for interval in intervals:
    cred_interval = list(t_post.interval(interval))
    cred_intervals.append([interval, cred_interval[0], cred_interval[1]])

cred_int_df = pd.DataFrame(
    cred_intervals,
    columns=['interval', 'lower value', 'uper value']
    ).set_index('interval')
cred_int_df


###############################################################################
# Pairwise comparison of all models: Frequentist approach
# -------------------------------------------------------
#
# We could also be interested in comparing the performance of all our models
# evaluated with GridSearchCV. In these case we are running our statistical
# test multiple times, which leads us to the `multiple comparisons
# problem <https://en.wikipedia.org/wiki/Multiple_comparisons_problem>`_.
#
# There are many possible ways to tackle the latter problem, but a standard
# approach is to apply a `Bonferroni
# correction <https://en.wikipedia.org/wiki/Bonferroni_correction>`_. The
# Bonferroni can be computed by multiplying the p-value by the number of
# comparisons we are testing.
#
# Let's compare the performance of the models using the corrected t-test:

n_comparisons = (
    factorial(len(model_scores))
    / (factorial(2) * factorial(len(model_scores) - 2))
    )
pairwise_t_test = []

for model_i, model_k in combinations(range(len(model_scores)), 2):
    model_i_scores = model_scores.iloc[model_i].values
    model_k_scores = model_scores.iloc[model_k].values
    differences = model_i_scores - model_k_scores
    t_stat, p_val = compute_corrected_ttest(
        differences, n, df, n_train, n_test
        )
    p_val *= n_comparisons  # implement Bonferroni correction
    if p_val > 1:
        p_val = 1  # Bonferroni can output p-values higher than 1
    pairwise_t_test.append(
        [model_scores.index[model_i], model_scores.index[model_k],
            t_stat, p_val]
        )

pairwise_comp_df = pd.DataFrame(
    pairwise_t_test,
    columns=['model_1', 'model_2', 't_stat', 'p_val']
    ).round(3)
pairwise_comp_df


###############################################################################
# We observe that after correcting for multiple comparisons, the only model
# that significantly differs from the others is `poly2`.
#
# `rbf`, the model ranked first by our GridSearch, does not significantly
# differ from `linear` or `poly3`.


###############################################################################
# Pairwise comparison of all models: Bayesian approach
# ----------------------------------------------------
#
# When using Bayesian estimation to compare multiple models, we don't need to
# correct for multiple comparisons (for reasons why see [4]_).
#
# We can carry out our pairwise comparisons the same way as in the first
# section:

pairwise_bayesian = []

for model_i, model_k in combinations(range(len(model_scores)), 2):
    model_i_scores = model_scores.iloc[model_i].values
    model_k_scores = model_scores.iloc[model_k].values
    differences = model_i_scores - model_k_scores
    t_post = t(
        df, loc=np.mean(differences),
        scale=correct_std(differences, n, n_train, n_test)
        )
    left_prob = t_post.cdf(rope_interval[0])
    right_prob = t_post.sf(rope_interval[1])
    rope_prob = t_post.cdf(rope_interval[1]) - t_post.cdf(rope_interval[0])

    pairwise_bayesian.append([left_prob, right_prob, rope_prob])

pairwise_bayesian_df = (pd.DataFrame(
    pairwise_bayesian,
    columns=['worse_prob', 'better_prob', 'rope_prob']
    ).round(3))

pairwise_comp_df = pairwise_comp_df.join(pairwise_bayesian_df)
pairwise_comp_df


###############################################################################
# Using the Bayesian approach we can compute the probability that a model
# performs better, worse or practically equivalent to another.
#
# Results show that the model ranked first by our GridSearch, `rbf`, has
# approximately a 6.8% chance of being worse than `linear`, and a 1.8% chance
# of being worse than `poly3`.
#
# `rbf` and `linear` have a 43% probability of being practically equivalent,
# while `rbf` and `poly3` have a 10% chance of being so.
#
# Similarly to the conclusions obtained using the Frequentist approach, all
# models have a 100% probability of being better than `poly2`, and none have a
# practically equivalent performance with the latter.


###############################################################################
# Take-home messages
# ------------------
# - When statistically comparing the performance of two models evaluated in
#   GridSearchCV, it is necessary to correct the calculated variance that could
#   be underestimated given that the scores of the models are not independent
#   from each other.
# - A Frequentist approach that uses a (variance-corrected) paired t-test can
#   tell us if the performance of one model is better that another with a
#   degree of certainty above chance.
# - A Bayesian approach can provide the probabilities of one model being
#   better, worse or practically equivalent than another. It can also tell us
#   how confident we are of knowing that the true differences of our models
#   fall under a certain range of values.
# - If multiple models are statistically compared, a multiple comparisons
#   correction is needed when using the Frequentist approach.


###############################################################################
# References:
# ___________
# .. [1] Dietterich, T. G. (1998). Approximate statistical tests for comparing
#        supervised classification learning algorithms. Neural computation,
#        10(7).
# .. [2] Nadeau, C., & Bengio, Y. (2000). Inference for the generalization
#        error. In Advances in neural information processing systems.
# .. [3] Bouckaert, R. R., & Frank, E. (2004). Evaluating the replicability of
#        significance tests for comparing learning algorithms. In Pacific-Asia
#        Conference on Knowledge Discovery and Data Mining.
# .. [4] Benavoli, A., Corani, G., Demšar, J., & Zaffalon, M. (2017). Time for
#        a change: a tutorial for comparing multiple classifiers through
#        Bayesian analysis. The Journal of Machine Learning Research, 18(1).
#           - See the Python library that accompanies this paper
#             `here <https://github.com/janezd/baycomp>`_.
