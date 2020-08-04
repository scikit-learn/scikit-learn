"""
==============================
Permutation Importance vs SHAP
==============================

In this example, we will compare :ref:`permutation_importance` with SHAP
(**SH**apley **A**dditive ex**P**lanations). Both are model-agnostic
model interpretation methods used to evaluate feature importances for models
trained on tabular data.

Permutation importance uses the decrease in model score after shuffling
the values of a feature, to measure how important each feature is. It can
be calculated using :func:`~sklearn.inspection.permutation_importance`. SHAP
calculates the contribution of each feature to the prediction output
by the model. It is based on Shapley values, which attribute a portion of
the prediction value to each feature in a fair way, as determined by game
theory. SHAP calculates the contribution of each feature for a particular
sample and uses the average across samples represent the average contribution
of a feature in a model. The authors of SHAP implement in the
Python library `SHAP <https://github.com/slundberg/shap>`_ which supports
scikit-learn models.

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
# 

# %%
Usage


pipeline
time

# %%
Interpretation

shap add up, perm doesnt add up
individual sample

perm training or test diff
error values
error ratio - comparable across diff problems

# %%
# Multicollinearity

# estimate the effect of
# withholding the feature of interest on the prediction output. It does this by
# taking a
# sample and calculating the effect of withholding a feature. Since this effect
# depends on the other features in the model, it calculates this difference
# (with and without the feature) for all possible combinations of the other
# features. This difference is averaged across all samples to 