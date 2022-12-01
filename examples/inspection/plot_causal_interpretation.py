"""
===================================================
Failure of Machine Learning to infer causal effects
===================================================

Machine Learning models are great for measuring statistical associations.
Unfortunately, unless we're willing to make strong assumptions about the data,
those models are unable to infer causal effects.

To illustrate this, we will simulate a situation in which we try to answer one
of the most important questions in economics of education: **what is the causal
effect of earning a college degree on hourly wages?** Although the answer to
this question is crucial to policy makers, `Omitted-Variable Biases
<https://en.wikipedia.org/wiki/Omitted-variable_bias>`_ (OVB) prevent us from
identifying that causal effect.


.. contents:: Table of Contents
   :local:
   :depth: 1
"""

print(__doc__)

# %%
# The dataset: simulated hourly wages
# -----------------------------------
#
# The data generating process is laid out in the code below. Work experience in
# years and a measure of ability are drawn from a Normal distribution; the
# hourly wage of one of the parents is drawn from Beta distribution. We then
# create an indicator of college degree which is positively impacted by ability
# and parental hourly wage. Finally, we model hourly wages as a linear function
# of all the previous variables and a random component. Note that all variables
# have a positive effect on hourly wages.

import numpy as np

n_samples = 10000
rng = np.random.RandomState(32)

experiences = rng.normal(20, 10, size=n_samples).astype(int)
experiences[experiences < 0] = 0
abilities = rng.normal(0, 0.15, size=n_samples)
parent_hourly_wages = 50 * rng.beta(2, 8, size=n_samples)
parent_hourly_wages[parent_hourly_wages < 0] = 0

college_degrees = (
    9 * abilities + 0.02 * parent_hourly_wages + rng.randn(n_samples) > 0.7
).astype(int)

hourly_wages = (
    0.2 * experiences
    + parent_hourly_wages
    + 2 * college_degrees
    + 5 * abilities
    + rng.normal(0, 1, size=n_samples)
)

hourly_wages[hourly_wages < 0] = 0

# %%
# Description of the simulated data
# ---------------------------------
#
# The following plot shows the distribution of each variable, and pairwise
# scatter plots. Key to our OVB story is the positive relationship between
# ability and college degree.

import pandas as pd
import seaborn as sns

df = pd.DataFrame(
    {
        "college degree": college_degrees,
        "ability": abilities,
        "hourly wage": hourly_wages,
        "experience": experiences,
        "parent hourly wage": parent_hourly_wages,
    }
)

grid = sns.pairplot(df, diag_kind="kde", corner=True)
# %%
# Predicting income with and without the ability feature
# ------------------------------------------------------
#
# Let's now train two :class:`~sklearn.linear_model.LinearRegression` models to
# predict our hourly wage. We include the ability feature in the first model
# and show that our estimate of the college degree coefficient is close to 2
# which is the true causal effect from our data generating process. In real
# life, intellectual ability is either never observed or only poorly measured
# (e.g. IQ score). Researchers are forced to "omit" the ability feature from
# their models, thereby inflating the estimate via a positive OVB. Despite an
# excellent R2 score, the model omitting the ability feature shows a
# coefficient that is far off the true value.

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

X = df[["experience", "parent hourly wage", "college degree", "ability"]]
y = df["hourly wage"]

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

regressor_with_ability = LinearRegression()
regressor_with_ability.fit(X_train, y_train)
y_pred_with_ability = regressor_with_ability.predict(X_test)
R2_with_ability = r2_score(y_test, y_pred_with_ability)

regressor_without_ability = LinearRegression()
regressor_without_ability.fit(X_train.drop(columns="ability"), y_train)
y_pred_without_ability = regressor_without_ability.predict(
    X_test.drop(columns="ability")
)
R2_without_ability = r2_score(y_test, y_pred_without_ability)

print(f"R2 score with ability: {R2_with_ability}")
print(f"College degree coefficient with ability: {regressor_with_ability.coef_[2]}")
print(f"R2 score without ability: {R2_without_ability}")
print(
    f"College degree coefficient without ability: {regressor_without_ability.coef_[2]}"
)


# %%
# Lessons learned
# ---------------
#
# Machine learning models are not designed for the estimation of causal
# effects. While we showed this with a linear model, OVB can affect any type of
# model.
#
# Whenever interpreting a coefficient or a change in predictions brought about
# by a change in one of the features, it is important to keep in mind
# potentially unobserved variables that could be correlated with both the
# feature in question and the target variable. Such variables are called
# `Confounding Variables <https://en.wikipedia.org/wiki/Confounding>`_. In
# order to still estimate causal effect in the presence of confounding,
# researchers usually conduct experiments in which the treatment variable (e.g.
# college degree) is randomized. When an experiment is prohibitively expensive
# or unethical, researchers can sometimes use other causal inference techniques
# such as `Instrumental Variables
# <https://en.wikipedia.org/wiki/Instrumental_variables_estimation>`_ (IV)
# estimations.
