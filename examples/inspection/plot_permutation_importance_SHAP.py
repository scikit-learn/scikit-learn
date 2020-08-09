"""
==============================
Permutation Importance vs SHAP
==============================

In this example, we will compare :ref:`permutation_importance` with SHAP
(**SH**apley **A**dditive ex**P**lanations). Both are model-agnostic
model interpretation methods used to evaluate feature importances for models
trained on tabular data.

Permutation importance uses the decrease in model score, after shuffling
the values of a feature, to measure how important each feature is. It can
be calculated using :func:`~sklearn.inspection.permutation_importance`.
Shapley values represent the contribution of each feature to the prediction
output by the model, using game theory. These values satisfy a number of
good properties and are thus deemed a fair way to 'split' the prediction output
between the features (for more on the properties see).
These values are computationally very expensive to calculate. SHAP are
a way to estimate these values using an additive model that is a linear
function of features.
SHAP allows for the calculation of the contribution of each feature for
individual samples. The average across samples is generally used to
represent the average contribution of a feature in a model. The authors of
SHAP implement a number of versions it in the Python library
`SHAP <https://github.com/slundberg/shap>`_. In this example we will discuss
the use of KernalSHAP and TreeSHAP, both of which support scikit-learn
estimators.

In practice, both methods will generally order features similarly in terms
of their importance in a model. However, there are some important
differences and practical considerations when using the methods, which are
outlined below.

.. topic:: References:

   [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32,
       2001. https://doi.org/10.1023/A:1010933404324
"""

# %%
# Data
# ----
#
# We will use the :ref:`california_housing_dataset` for this example. In this
# dataset the target is the median house price (in $100,000's) for a district
# and the features consist of various information about each district. We
# will only use a subset of the data to speed up computation.

import pandas as pd

from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split

cal_housing = fetch_california_housing()
y = cal_housing.target[::10]
X = pd.DataFrame(
    data=cal_housing.data[::10, :], columns=cal_housing.feature_names
)
X.head()

# %%
# Split the data into training and testing subsets.

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=7)

# %%
# Calculating feature importance
# ------------------------------
#
# We will used a :class:`~sklearn.ensemble.RandomForestRegressor` and fit the
# model with our training subset.

from sklearn.ensemble import RandomForestRegressor

reg = RandomForestRegressor(max_depth=10, random_state=7)
reg.fit(X_train, y_train)

# %%
# Permutation importance
# ^^^^^^^^^^^^^^^^^^^^^^
#
# First, we will calculate permutation importances on the test subset
# using the default score metric of
# :class:`~sklearn.ensemble.RandomForestRegressor`, :ref:`R^2 <r2_score>`.
# The values of each feature will be permuted ``n_repeats=10`` times and the
# decrease in R^2 value for each permutation is shown below with boxplots.

import matplotlib.pyplot as plt

from sklearn.inspection import permutation_importance

perm_import = permutation_importance(reg, X_test, y_test, n_repeats=10,
                                     random_state=42, n_jobs=2)
sorted_idx = perm_import.importances_mean.argsort()

fig, ax = plt.subplots()
ax.boxplot(
    perm_import.importances[sorted_idx].T, vert=False,
    labels=X_train.columns[sorted_idx]
)
ax.set_title("Permutation Importances")
fig.tight_layout()
plt.show()

# %%
# The plot shows that ``MedInc`` (median income) causes by far the biggest
# drop in R^2 score whereas ``AveBedrms`` and ``Population`` seem to have
# almost no effect.
#
# **Considerations**
#
# * Permutation importance is linked to the score metric used, thus selecting
#   an appropriate metric for your needs is essential.
# * Permutation importance can be calculated with either the training or
#   testing subset and each provides different information about the model.
#   Using the training subset shows what the model has learned to use for
#   making predictions. Using the testing subset shows what is actually useful
#   when making predictions on novel data. Features that the model is able to
#   use to overfit on the training data will appear more important when using
#   the training data than the testing data. This can be seen in
#   :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`
#   where random features are added into the dataset.
#
# **Advantages**
#
# * Permutation importances are easy to interpret and provide information on
#   the global importance of each feature.
# * :func:`~sklearn.inspection.permutation_importance` can accept both
#   an estimator and a :class:`sklearn.pipeline.Pipeline`. When given a
#   pipeline (that includes feature transformations), permutations are
#   performed on the original feature values. This means that
#   the importances are interpretable in the original feature space. If
#   only estimators were accepted, like with the SHAP package, the
#   permutations would need to be performed on transformed features meaning
#   the importance values would need to be interpreted in the transformed
#   feature space.
# * Permutation importances are error ratios, where 1 means no
#   difference between error with and without permutation. This means
#   these values are comparable between different models and problems.
#
# **Disadvantages**
#
# * Feature permutation destroys all dependency between features. This means
#   the permutation importance value represents the combination of
#   the main effect and any interaction effects with other features. Take
#   for example two features that interact with each other. The interaction
#   effects would be 'accounted for' twice, once in the importance value of
#   each feature. This means that the decrease in score from permuting each
#   feature would not add up to equal to overall decrease in score from
#   permuting all features, unless there were no interaction effects between
#   features.
# * Permutation importance assumes independence between features. The effect
#   is that the 'importance' is split between correlated features. See
#   :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`
#   for an example of this effect.
#
# KernalSHAP
# ^^^^^^^^^^^
#
# KernalSHAP is a kernal-based method to estimate Shapley values.
# It calculates the prediction output with different subsets of the
# features 'missing'. Missingness is simulated by using a 'background' (e.g.,
# average) value for that feature. These predictions are then used to fit a
# linear model whose predictions match that of the original model as closely
# as possible. The coefficients of the linear model are the Shapley values.
#
# This is much more computationally expensive than permutation importance
# as the number of possible combinations of features quickly becomes very
# large as the number of features increases. However, this does enable SHAP
# to account for interaction effects between features.
#
# First, we will instantiate ``KernalSHAP`` using our fitted regressor,
# ``reg`` and some 'background' data, used to simulate missingness.
# If your dataset is small (e.g., training data is <100 samples) you can use
# the whole training subset for the background data. Since our dataset is
# larger than this, we will use  to summary values-
#
#

from shap import KernalSHAP
from shap import

kernal_ex = KernalSHAP(reg, )


# %%
# Next we will calculate Shapley values using the testing subset.
#




# Explain background data. Show expected value. shap value show individual
# make plot-




# %%
# **Advantages**
#
# * This method allows you to compute feature importances for individual
#   samples.
# * Interaction effects are dealt with better than in permutation importances.
#   Contribution of each feature add up to the overall prediction.
#
# **Disadvantages**
#
# * The method is very computationally expensive, especially when there are
#   are a large number of features.
# * This method also assumes independence between features. Again, the result
#   is that 'contributions' will be split between correlated features.

# Tree SHAP
# ^^^^^^^^^


#
# Both permutation importance and SHAP assume independence between features.
#
# estimate the effect of
# withholding the feature of interest on the prediction output. It does this by
# taking a
# sample and calculating the effect of withholding a feature. Since this effect
# depends on the other features in the model, it calculates this difference
# (with and without the feature) for all possible combinations of the other
# features. This difference is averaged across all samples to
