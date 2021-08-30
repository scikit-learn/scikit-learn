"""
==========================================================
Inspect and analyze a linear model within cross-validation
==========================================================

Evaluating a predictive linear model involves :ref:`cross-validation
<cross_validation>`. This example:

* details how to interpret the results from a cross-validation framework;
* highlights how to inspect the internals of a model when using the
  cross-validation framework.
"""

print(__doc__)
import sklearn

sklearn.set_config(display="diagram")

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

# %%
# Dataset
# -------
#
# We will use the :ref:`california_housing_dataset` where the goal is to
# predict the average house value in a neighborhood. From the start, we will
# split our data into two sets: a set that we will use for all our experiments
# and a set that we will leave out for further confirmation.

from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split

X, y = fetch_california_housing(as_frame=True, return_X_y=True)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)

# %%
# Predictive model
# ----------------
#
# In this example, we will use a linear model which should be a good baseline:
# a :ref:`ridge model <ridge_regression>`. A ridge model enforces a L2 penalty
# on the coefficients. Thus, the penalty parameter `alpha` has to be tuned.
# More importantly, this parameter needs to be tuned for our specific problem:
# tuning on another dataset does not ensure an optimal parameter value for the
# current dataset.
#
# Here, we use the class :class:`~sklearn.linear_model.RidgeCV` that can tune
# `alpha` by internal leave-one-out cross-validation when calling `fit`.
#
# We also add a preprocessing stage to :ref:`standardize
# <preprocessing_scaler>` the data such that the regularization strength is
# applied homogeneously on each coefficient.
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

alphas = np.logspace(-1, 2.5, num=50)
model = Pipeline(
    [
        ("scaler", StandardScaler()),
        ("ridge_regressor", RidgeCV(alphas=alphas)),
    ]
)

# %%
# `model` is a machine learning pipeline made of the preprocessing stage
# followed by a ridge regressor that will perform an internal cross-validation
# at `fit` to select the best parameter `alpha`. The candidates for the
# parameter `alpha` are given by the variable `alphas`.
#
# Cross-validation framework
# --------------------------
#
# Before putting such a predictive model into production, one needs to evaluate
# the performance of the model to have an idea of what to expect in production.
#
# Cross-validation should be used to make this analysis. First, it allows us to
# quantify the variance of the model performance. A large variance of the score
# metric will indicate that we cannot trust the reported performance nor try to
# interpret findings built on internal model's parameters. One possible cause
# of large variations are small sample sizes.
from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedKFold

cv = RepeatedKFold(n_splits=10, n_repeats=10, random_state=0)
cv_results = cross_validate(
    model, X_train, y_train, cv=cv, return_estimator=True, n_jobs=2
)

# %%
# Here, we used a repeated K-fold cross-validation. At each round of
# cross-validation, it should be noted that the parameter `alpha` of the ridge
# regressor is also optimized via another internal cross-validation. This
# process is called a nested cross-validation and should always be implemented
# whenever model's parameters need to be optimized.
#
# Analyze model performance
# .........................
#
# As previously mentioned, one should look at the variance of the model
# performance within the cross-validation framework.
import matplotlib.pyplot as plt

cv_score = cv_results["test_score"]
plt.hist(cv_score, bins=200, density=True)
plt.xlim([0, 1])
plt.ylabel("Density")
plt.xlabel("R2 score")
_ = plt.title(
    "Distribution of the scores on the test sets\n during the cross-validation"
)

# %%
# We start by plotting the empirical distribution of the test score computed
# during cross-validation.
#
# We observe little variation in R2 score and can thus safely use the results
# to interpret the model performance.
#
# Our baseline model has an R2 of around 0.6 which is better than a dummy
# regressor. We can therefore use this predictive model as a baseline against
# more advanced machine learning pipelines.
#
# To conclude, cross-validation allows us to answer to two questions: are the
# results reliable and, if it is, how good is the algorithm used to create
# the different models used during the cross-validation?
#
# Model inspection
# ................
#
# Once we are happy or at least aware of any model limitations, we can
# investigate the internals of our model. When we performed the
# cross-validation, we set the parameter `return_estimator` to `True`.
#
# It allows us to get the different predictive models trained and tested within
# cross-validation.
cv_pipelines = cv_results["estimator"]

# %%
# While the cross-validation allows us to know if our models are reliable, we
# have not checked if the model at each cross validation fold is similar.
#
# Recall that in our case, our baseline model required tuning of the parameter
# `alpha`. Therefore, up-to-now, we do not know if the optimal `alpha`
# parameter for all models are similar. In other words, we are asking the
# question: what is the variance of the `alpha` parameter value across
# iterations.
#
# Indeed, if the `alpha` parameters vary depending on the input data, it
# might be challenging to put our model in production because we would not be
# able to fix this hyperparameter.
#
# Let's check the `alpha` parameter variance.
alpha = [pipeline["ridge_regressor"].alpha_ for pipeline in cv_pipelines]
plt.hist(alpha, bins=30, density=True)
plt.xlabel("Alpha")
plt.xscale("log")
plt.ylabel("Density")
_ = plt.title("Distribution of alpha parameter \nduring cross-validation")

# %%
# We see that the regularization parameter, `alpha`, values are centered and
# condensed around 40. This is a good sign and means that most of the models
# tuned within the cross-validation had similar `alpha` values.
#
# However, not only hyperparameter such as `alpha` should be studied. The model
# parameter coming out of the fitting process should analyzed. In our case, we
# used a linear model. These models are parametrized and rely on two
# parameters: `coef_` and `intercept_`. Therefore, we should analyze the
# variance of these parameters as well.
#
# For the sake of simplicity, we are going to solely look at the `coef_`.
import pandas as pd

coefficients = pd.DataFrame(
    [pipeline["ridge_regressor"].coef_ for pipeline in cv_pipelines],
    columns=X.columns,
)
coefficients

# %%
import seaborn as sns

plt.figure(figsize=(9, 7))
sns.stripplot(data=coefficients, orient="h", color="k", alpha=0.5)
sns.boxplot(data=coefficients, orient="h")
plt.axvline(x=0, color=".5")
plt.title("Coefficient values our model")
_ = plt.subplots_adjust(left=0.3)

# %%
# We observe that the coefficients do not vary much taking the size of the
# dataset into account, meaning that all the trained models are similar. Each
# individual model is expected to more or less give the same predictions.
#
# This would not necessarily have been the case if the number of training
# samples had been smaller, if the model had not been properly regularized or
# if we had strongly dependent features in the data.
#
# PFinalizing the predictive model
# --------------------------------
#
# In the above analysis, we saw that the internally tuned value of `alpha`
# did not vary much across cross-validation splits.
# Indeed, we can retrain this model and optimize the value of `alpha` on the
# full training set.
from sklearn.base import clone

final_model = clone(model)
final_model.fit(X_train, y_train)

# %%
# At the beginning of the process, when we split our data, we had a set of left
# out data. Now, we can use it to further check if the model performs as we
# would expect from the analysis done within the cross-validation framework.
print(f"The performance of our final model: R2={final_model.score(X_test, y_test):.2f}")

# %%
# We see that the statistical performance is comparable to the cross-validation
# statistical performances which is not surprising. Similarly, we could look at
# the coefficients of the final model and compare them with the coefficients
# obtained within the cross-validation study.
#
# However, you should be aware that this latest step does not give any
# information about the variance of our final model. Thus, if we want to
# evaluate our final model, we should get a distribution of scores if we would
# like to get this information.
#
# The example
# :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_stats.py`
# gives more information regarding the comparison that can be made during a
# :class:`~sklearn.model_selection.GridSearchCV`.
