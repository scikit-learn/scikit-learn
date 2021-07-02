"""
================================================
Causal Interpretation of Machine Learning Models
================================================

Machine Learning models are great for measuring statistical associations.
Unfortunately, unless we're willing to make strong assumptions about the data,
those models are unable to infer causal effects.

To illustrate this, we will simulate a situation in which we try to answer
one of the most important questions in economics of education:
what is the causal effect of earning a college degree on hourly wages?
Although the answer to this question is crucial to policy makers,
`Omitted-Variable Biases <https://en.wikipedia.org/wiki/Omitted-variable_bias>`_ (OVB)
prevent us from identifying that causal effect.


.. contents::
   :local:
   :depth: 1
"""

print(__doc__)

import numpy as np
from sklearn.linear_model import LinearRegression

# %%
# The dataset: simulated hourly wages
# -----------------------------------
#
# The data generating process is laid out in the code below.
# Work experience in years and a measure of ability are drawn
# from a Normal distribution;
# the hourly wage of one of the parents is drawn from Beta distribution.
# We then create an indicator of college degree which is positivily
# impacted by ability and parental hourly wage. Finally, we model
# hourly wages as a linear function of all the previous variables and a
# random component. Note that all variables have a positivie effect on
# hourly wages.

N = 10000
experiences = np.random.normal(20, 10, size=N).astype(int)
experiences[experiences < 0] = 0
abilities = np.random.normal(0, 0.15, size=N)
parent_hourly_wages = 50 * np.random.beta(2, 8, size=N)
parent_hourly_wages[parent_hourly_wages < 0] = 0


college_degrees = (9 * abilities + 0.02 * parent_hourly_wages
                   + np.random.randn() > 0.7).astype(int)

hourly_wages = 0.2 * experiences + parent_hourly_wages + 2 \
    * college_degrees + 5 * abilities + np.random.normal(0, 1, size=N)

hourly_wages[hourly_wages < 0] = 0

# %%
# Description of the simulated data
# ---------------------------------
#
# The following plot shows the distribution of each variable,
# and pairwise scatter plots. Key to our OVB story is the positive
# relationship between ability and college degree.

import pandas as pd
import seaborn as sns

df = pd.DataFrame({
    'hourly_wage': hourly_wages,
    'experience': experiences,
    'parent_hourly_wage': parent_hourly_wages,
    'college_degree': college_degrees,
    'ability': abilities
    })

_ = sns.pairplot(df, diag_kind='kde')

# %%
# Predicting income with and without the ability feature
# ------------------------------------------------------
#
# Let's now train two :class:`~sklearn.linear_model.LinearRegression`
# models to predict our hourly wage. We include the ability feature
# in the first model and show that our estimate of the college degree
# coefficient is close to 2 which is the true causal effect from our
# data generating process. In real life, intellectual ability is either
# never observed or poorly measured (e.g. IQ score).
# Researchers are forced to "omit" the ability feature from their models,
# thereby inflating the estimate via a positive OVB.

clf_with = LinearRegression()
clf_with.fit(np.stack([experiences, parent_hourly_wages, college_degrees,
             abilities], axis=1), hourly_wages)

clf_without = LinearRegression()
clf_without.fit(np.stack([experiences, parent_hourly_wages, college_degrees],
                axis=1), hourly_wages)


print('College degree coefficient with ability control: {}'.format(
    clf_with.coef_[2]))

print('College degree coefficient without ability control: {}'.format(
    clf_without.coef_[2]))

# %%
# Lessons learned
# ---------------
#
# Machine learning models are not reliable for inferring causation causal effect.
# While we showed this with a linear model, OVB can strike any type of models.
#
# Whenever interpreting a coefficient or a change in predictions
# brought about by a change in one of the features, it is important to think of
# potentially unobsered variables and whether they could be correlated with both
# the feature in question and the target variable.
# Such variables are called `Confounding Variables
# <https://en.wikipedia.org/wiki/Confounding>`_.
# To avoid these, researchers usually conduct experiments in which
# the treatment variable (e.g. college degree) is randomized. When an experiment
# is prohibitively expensive or unethical, researchers
# can somethimes use other causal inference techniques such
# as `Instrumental Variables
# <https://en.wikipedia.org/wiki/Instrumental_variables_estimation>`_ (IV)
# estimations.
