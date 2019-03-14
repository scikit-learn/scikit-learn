"""
==========================================================
Permutation Importance vs Random Forest Feature Importance
==========================================================

The permutation importance of a feature is how much a model's score decreases
when the feature's values are permuted. [1]_ In this context, the more
important features will cause the model's score to decrease the most. If a
feature's importance is negative, the model improved when the feature is
permutated. This suggests that the model would benefit from removing the
feature.

.. topic:: References:

   .. [1] Breiman, L. Machine Learning (2001) 45: 5.
        https://doi.org/10.1023/A:1010933404324
   .. [2] Strobl, C., Boulesteix, AL., Zeileis, A. et al. BMC Bioinformatics
        (2007) 8: 25. https://doi.org/10.1186/1471-2105-8-25
"""
print(__doc__)
##############################################################################
# In this example, the :class:`sklearn.ensemble.RandomForestRegressor`'s feature
# importance is compared with the permutation importance using the titantic
# dataset. A combination of numerical features and categorical features were
# used to train the initial model, ``rfr_init``:
from sklearn.inspect import permutation_importance
print("what is going on")

##############################################################################
# The tree based feature importance ranks the numerical features, age and BLAH,
# to be the most important features:

##############################################################################
# The corresponding test score for, ``rfr_init``, is:

##############################################################################
# Next, a model, ``rfr_new``, is trained with the numerical features removed
# and its test score is calcuated.

##############################################################################
# When the "important" numerical features are removed, the score increases! A
# decision tree's feature importance is known to be inflated for continuous
# variables. [2]_ As an alternatie, the permutation permutation of
# ``rfr_init`` is computed, which shows that the categorical feature, ``sex```
# is the most important feature:
