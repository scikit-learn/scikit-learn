"""
=====================================================
Inspect and analyze estimator within cross-validation
=====================================================

Evaluating a predictive model involves cross-validation. This example:

* recalls the question that a cross-validation framework allows to answer to;
* highlights the way to inspect the internals of a model when using the
  cross-validation framework.
"""

print(__doc__)

# %%
# Dataset
# -------
#
# We will use the california housing dataset where the goals is to predict the
# average house value in a neighborhood. From the start, we will split our data
# into two sets: a set that we will use to make all our experiments and a set
# that we will leave out for further confirmation.

# %%
from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split

X, y = fetch_california_housing(as_frame=True, return_X_y=True)
y -= y.mean()

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.1, random_state=0
)

# %%
# Predictive model
# ----------------
#
# In this example, we will use a linear model which should be a good baseline:
# a ridge model. A ridge model enforce a L2 penalty on the coefficients. Thus,
# the penalty parameter `alpha` has to be tuned. More importantly, this
# parameter needs to be tune for our specific problem: a fine tuning on another
# dataset does not ensure an optimal parameter value for the current dataset.
#
# Here, we use the class :class:`~sklearn.linear_model.RidgeCV` that allows to
# tune `alpha` by cross-validation.
#
# Besides, we add a preprocessing stage to standardize the data such that the
# optimization problem encountered by the ridge regressor is well-posed.

# %%
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

alphas = np.logspace(-1, 2.5, num=50)
model = make_pipeline(StandardScaler(), RidgeCV(alphas=alphas))

# %%
# `model` is a machine learning pipeline made of the preprocessing stage
# followed by a ridge regressor that will perform an internal cross-validation
# at `fit` to select the best parameter `alpha`. The candidates for the
# parameter `alpha` are given by the variable `alphas`.
#
# Cross-validation framework
# --------------------------
#
# Before to put such a predictive model into production, one needs to evaluate
# the performance of such model to have an idea of what to expect in
# production.
#
# Cross-validation should be used to make such analysis. First, it allows to
# quantify the variance of the model performance. A large variance of the
# metric will indicate that we cannot trust the reported performance nor try to
# interpret findings built on internal model's parameters. Usually large
# variations are linked to small sample size but not only.

# %%
from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedKFold

cv = RepeatedKFold(n_splits=10, n_repeats=10, random_state=0)
cv_results = cross_validate(
    model,
    X_train,
    y_train,
    cv=cv,
    return_estimator=True,
    n_jobs=2,
)

# %%
# Here, we used a repeated K-fold cross-validation. At each round of the
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

# %%
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
# We observe a little variation in terms of R2 score. Thus, we can now safely
# interpret if the results obtained are synonymous of having a good model.
#
# Our baseline perform around an R2 of 0.6 which is better than a dummy
# regressor. Therefore, such a predictive model can be used as a baseline if
# we would like to develop more advanced machine learning pipelines.
#
# To conclude, cross-validation allows to answer to two questions here: are the
# results reliable and, if it is the case, how good is my predictive model.
#
# Model inspection
# ................
#
# Once we are happy or at least aware of our model limitation, we can
# investigate the internals of our model. When we performed the
# cross-validation, we pass a parameter `return_estimator` to `True`.
#
# It allows to get the different predictive models trained and tested within
# cross-validation.

# %%
cv_estimators = cv_results["estimator"]

# %%
# While the cross-validation allows us to know if we got reliable models, we
# did not check if the model at each round were similar.
#
# We recall that in our case, our baseline was required to tune tha parameter
# `alpha`. Therefore, up-to-now, we do not know if the optimal `alpha`
# parameter for all models are similar. We can reformulate this as: what is the
# variance of the `alpha` parameter values across iterations.
#
# Indeed, if the parameters `alpha` are varying depending of the input data, it
# might be challenging to put our model in production because we will not be
# able to fix this hyperparameter.
#
# Let's check the `alpha` parameter variance.

# %%
alpha = [est[-1].alpha_ for est in cv_estimators]
plt.hist(alpha, bins=30, density=True)
plt.xlabel("Alpha")
plt.ylabel("Density")
_ = plt.title("Distribution of alpha parameter \nduring cross-validation")

# %%
# We see that the regularization parameter `alpha` values are centered and
# condensed around 40: it is a good sign. It means that most of the models
# tuned within the cross-validation had an equivalent parameter `alpha`.
#
# However, not only hyperparameter such as `alpha` should be studied. The model
# parameter coming out of the fitting process should analyzed. In our case, we
# used a linear model. These models are parametrized and rely on two
# parameters: `coef_` and `intercept_`. Therefore, we should analyze the
# variance of these parameters as well.
#
# For the sake of simplicity, we are going to solely look at the `coef_`.
# %%
import pandas as pd

coefficients = pd.DataFrame(
    [est[-1].coef_ for est in cv_estimators],
    columns=X.columns,
)
coefficients

# %%
import seaborn as sns

plt.figure(figsize=(9, 7))
sns.swarmplot(data=coefficients, orient="h", color="k", alpha=0.5)
sns.boxplot(data=coefficients, orient="h")
plt.axvline(x=0, color=".5")
plt.title("Coefficient variability")
_ = plt.subplots_adjust(left=0.3)

# %%
# We observe that the coefficients do not vary meaning that all models trained
# are similar. Each individual model is expected to more or less give the same
# predictions.
#
# Put a predictive model in production
# ------------------------------------
#
# With the above analysis, we can safely create a model by fixing the `alpha`
# hyperparameter. Subsequently, we can train the model on the full training
# set.

# %%
from sklearn.linear_model import Ridge

production_model = make_pipeline(StandardScaler(), Ridge(alpha=40))
production_model.fit(X_train, y_train)

# %%
# At the beginning of the process, we left out some data. Now, we can use it
# to further check if the model performance is similar to what we could expect
# from the analysis done within the cross-validation framework.

# %%
print(
    f"The performance of our production model: "
    f"R2={production_model.score(X_test, y_test):.2f}"
)

# %%
# We see that the performance are comparable to the above performance which is
# not surprising. Similarly, we could look at the coefficients of the
# production model and compare it with the coefficients obtained within the
# cross-validation study.
#
# However, you should be aware that this latest step does not give any
# information about the variance of the model. It should never be used to
# evaluate the model itself.
