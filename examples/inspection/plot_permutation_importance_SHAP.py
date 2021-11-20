r"""
==============================
Permutation Importance vs SHAP
==============================

In this example, we will compare :ref:`permutation_importance` with SHAP
(**SH**\ apley **A**\ dditive ex\ **P**\ lanations). Both are model-agnostic
model interpretation methods used to evaluate feature importances for models
trained on tabular data.

Permutation importance uses the decrease in model score, after shuffling
the values of a feature, to measure how important each feature is. It can
be calculated using :func:`~sklearn.inspection.permutation_importance`.
Shapley values represent the contribution of each feature to the prediction
output by the model, using game theory. It is calculated as the average
marginal contribution of each feature, across all possible feature subset
combinations. These values satisfy a number of good properties and are thus
deemed a fair way to 'split' the prediction output between the features (for
more on these properties see [1]_ and [5]_).
These values are computationally very expensive to calculate. SHAP offers
a way to estimate these values using an additive model that is a linear
function of features. SHAP calculates the contribution of each feature
for individual samples. The average of these contributions across samples
can then be used to indicate the 'importance' of a feature in a model.
Several versions of SHAP are implemented in the Python library
`SHAP <https://github.com/slundberg/shap>`_. In this example we will discuss
the use of KernelSHAP [1]_ and TreeSHAP [2]_, both of which support
scikit-learn estimators.

In practice, both methods will generally order features similarly in terms
of their importance in a model. However, there are some important
differences and practical considerations when using the methods, which are
outlined below.
"""

# %%
# Data
# ----
#
# We will use the :ref:`california_housing_dataset` for this example. In this
# dataset the target is the median house price (in $100,000's) for a district
# and the 8 features consist of various information about each district. We
# will only use a subset of the data to speed up computation.

import pandas as pd

from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split

cal_housing = fetch_california_housing()
y = cal_housing.target[::10]
X = pd.DataFrame(data=cal_housing.data[::10, :], columns=cal_housing.feature_names)
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
# :class:`~sklearn.ensemble.RandomForestRegressor`: :ref:`R² <r2_score>`.
# The values of each feature will be permuted ``n_repeats=10`` times and the
# decrease in R^2 value for each permutation is shown below with boxplots.

import matplotlib.pyplot as plt

from sklearn.inspection import permutation_importance

perm_import = permutation_importance(
    reg, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)
sorted_idx = perm_import.importances_mean.argsort()

fig, ax = plt.subplots()
ax.boxplot(
    perm_import.importances[sorted_idx].T,
    vert=False,
    labels=X_train.columns[sorted_idx],
)
ax.set_title("Permutation Importances")
fig.tight_layout()
plt.show()

# %%
# The plot shows that ``MedInc`` (median income) causes by far the biggest
# drop in R² score whereas ``AveBedrms`` and ``Population`` seem to have
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
#   an estimator and a :class:`~sklearn.pipeline.Pipeline`. When given a
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
#   for an example of this.
#
# KernelSHAP
# ^^^^^^^^^^^
#
# KernelSHAP is a kernel-based method to estimate Shapley values for
# individual samples.
# First, it calculates the predictions for a sample when different subsets of
# the features are 'missing'. Missingness is simulated by using a 'background'
# (e.g., average) value for that feature. These predictions are then used to
# fit a linear model whose predictions match that of the original model as
# closely as possible. The coefficients of the linear 'explanation' model are
# the Shapley values. The linear explanation model equation is:
#
# .. math::
#   g(z')=\phi_0 + \sum_{i=1}^{M} \phi_i z'
#
# where :math:`g(z')` is the explanatory model, :math:`\phi_0` is the model
# prediction when all features are 'missing', :math:`M` is the number of
# possible subset sizes (n_features - 1), :math:`z'\in\{0,1\}^M`
# (denotes the presence or absence of each feature) and :math:`\phi_i` is the
# estimated Shapley value.
#
# This is much more computationally expensive than permutation importance
# particularly as the number of possible combinations of features for all
# feature subsets sizes quickly becomes very large as the number of features
# increases. However, this does enable SHAP to account for interaction effects
# between features.
#
# First, we will instantiate ``KernelSHAP`` using our fitted regressor,
# ``reg`` and some background data, which will be used to simulate missingness.
# If your dataset is small (e.g., training data is <100 samples) you can use
# the whole training subset for the background data. Our data is larger than
# this so we must summarize it in some way, otherwise computation will be
# too slow. For simplicity, we will use the
# median values of our features but SHAP offers a ``kmeans`` function that
# can summarize each feature with ``k`` means.

import shap

med = X_train.median().values.reshape((1, X_train.shape[1]))
explainer = shap.KernelExplainer(reg.predict, med)

# %%
# ``explainer`` stores various information about the data as attributes. Of
# interest is ``expected_value``. This represents :math:`\phi_0` in our
# equation above.

explainer.expected_value

# %%
# Can you work out how this value was calculated?
#
# It is the prediction output by our model when using the background
# values we gave it.

reg.predict(med)

# %%
# Next we will calculate Shapley values using the testing subset. This is the
# computationally expensive step.

shap_values = explainer.shap_values(X_test, l1_reg="aic")

# %%
# Let's look at the Shapley values of one sample. There are 8 Shapley values,
# one for each feature. '0' means that the feature did not contribute to the
# prediction output. A negative value 'pushes' the prediction lower and a
# positive value pushes the prediction higher.

shap_values[0, :]

# %%
# The Shapley values should sum to the difference between the
# prediction output by our our model ``reg`` and  the ``expected_value`
# (depending on how well the linear model was able to be fit).

print(f"The sum of Shapley values: {shap_values[0, :].sum()}")
prediction = reg.predict(X_test.iloc[0, :].values.reshape(1, -1))[0]
print(
    "Difference between prediction and expected value: "
    f"{prediction - explainer.expected_value}"
)

# %%
# We can also plot the Shapley values:

import matplotlib.pyplot as plt

shap.summary_plot(shap_values, X_test, show=False)
plt.tight_layout()

# %%
# In the plot above, each dot represents the Shapley value of one sample,
# for that feature. The features are ordered from most important at the top
# to least important at the bottom. Note that the dots cluster around 0
# (no contribution) more and more as you go down.
#
# Additionally, if you compare with the permutation importance plot, you will
# notice that the order of the features is roughly the same.
#
# **Advantages**
#
# * This method allows you to compute feature importances for individual
#   samples.
# * Interaction effects are accounted for, unlike in permutation importances.
#   For each sample, the contribution of each feature add up to the overall
#   prediction, as shown above.
#
# **Disadvantages**
#
# * The method is very computationally expensive, especially when there
#   are a large number of features.
# * This method also assumes independence between features [4]_. Similar to
#   permutation importance, the result is that 'contributions' will be split
#   between correlated features.
#
# TreeSHAP
# ^^^^^^^^
#
# TreeSHAP is another variant of SHAP designed for tree-based models. It is
# much faster than KernelSHAP and uses conditional expectation instead of
# marginal expectation. For a single tree, the expectation conditioned
# on a subset of features is the average value of all 'reachable' nodes,
# weighted by the number of samples in each node. 'Reachable' means
# nodes that are not contradicted by the values of the features present.
# For example, if a node splits on a feature that is missing, both child
# nodes are reachable. Conversely, if a node splits on a feature
# that is present, only the child that satisfies the split condition is
# 'reachable'. The difference between the conditional expectation of feature
# subsets with and without the feature of interest is used to estimate
# Shapley values.
#
# The major advantage of TreeSHAP is its time complexity.
# Compared to KernelSHAP, which computes Shapley values
# in exponential time, TreeSHAP does this in polynomial time [2]_.
#
# Below we calculate Shapley values using ``TreeExplainer``. We do not have
# to provide a 'background' dataset as the model can use the number of
# training samples at each node/leaf of the tree. This information is
# stored in the tree model. Note that the background dataset is used
# differently here. Here it provides the distribution of the features used
# to calculated expected values.

explainer = shap.TreeExplainer(reg)
shap_values = explainer.shap_values(X_test)

shap.summary_plot(shap_values, X_test, show=False)
plt.tight_layout()

# %%
# Note that the order of features is exactly the same as that calculated
# with KernelShap. TreeSHAP is able to calculate Shapley values much faster
# though.
#
# TreeSHAP Shapley values do have some problems:
#
# * They are sensitive to the degree of sparsity, which often arises when
#   features are continuous. This means that the calucalted Shapley values
#   "are very sensitive to noise in the data" [3]_.
# * A 'dummy' feature, which is not used by the model but is correlated
#   to a 'useful' feature, can have a non-zero Shapley value [3]_.
#
# References
# ----------
#
# .. topic:: References:
#
#   .. [1] Lundberg, Scott M., and Su-In Lee. :arxiv:`"A unified approach to
#          interpreting model predictions."
#          <1903.10464>` Advances in Neural
#          Information Processing Systems. 2017.
#   .. [2] Lundberg, Scott M., Gabriel G. Erion, and Su-In Lee.
#          :arxiv:`"Consistent individualized feature attribution for tree
#          ensembles." <1802.03888>` arXiv
#          preprint arXiv:1802.03888 (2018).
#   .. [3] Sundararajan, Mukund, and Amir Najmi. :arxiv:`"The many Shapley
#          values for model explanation."
#          <1908.08474>` arXiv:1908.08474 (2019).
#   .. [4] Aas, Kjersti, Martin Jullum, and Anders Løland. :arxiv:`"Explaining
#          individual predictions when features are dependent: More accurate
#          approximations to Shapley values."
#          <1903.10464>` arXiv:1903.10464 (2019).
#   .. [5] Molnar, Christoph. "Interpretable machine learning. A Guide for
#          Making Black Box Models Explainable", 2019.
#          https://christophm.github.io/interpretable-ml-book/
