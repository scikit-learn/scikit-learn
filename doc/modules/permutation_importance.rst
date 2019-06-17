
.. _permutation_importance:

Permutation feature importance
==============================

.. currentmodule:: sklearn.inspection

Permutation feature importance is a model inspection technique that can be used
for any `fitted` `estimator` when the data is rectangular. This is especially 
useful for non-linear or opaque `estimators`. The permutation feature 
importance is defined to be the decrease in a model score when the feature 
value is randomly shuffled [1]_. This procedure breaks the relationship between 
the feature and the target, thus the drop in the model score is analogous to 
how much the model depends on the feature. This technique benefits from being 
model agnostic and can be calculated many times with different permutations of 
the feature.

The :func:`permutation_importance` function calculates the feature importance
of `estimators` for a given dataset. The ``n_rounds`` parameter sets the number
of times a feature is randomly shuffled and returns a sample of feature
importances. Note that :func:`permutation_importance` can be computed on any
dataset, though it is common practice to use a hold-out set since estimates
on the training set might not be accurate. Feature importance based on
the training set can inflate the importances, when in fact, the model is
overfitting and the feature is not important. Feature importance based on the
test set gives the importance of a feature on unseen data, which gives a sense
of the features that are *actually* important.

Relation to feature importance in trees
---------------------------------------

Tree based models provides a different measure of own feature importances based
on the mean decrease in the criterion. This gives importances to features
that may not be predictive on unseen data. The permutation feature importance
avoids this issue, since it can be applied to unseen data. Additionally, the
permtuation feature importance uses an arbitrary metric on the tree's 
predictions. These two methods of obtaining feature importance is explored in:
:ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`.

Strongly correlated features
----------------------------

When two features are correlated and one of the features is permuted, the model
will still have access to the feature through its correlated feature. This will 
result in a lower importance for both features, where they might *actually* be
important. One way  to handle this is to cluster features that are correlated
and only keep one feature from each cluster. This use case is explored in: 
:ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`.

.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`
  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`

.. topic:: References:

   .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32,
       2001. https://doi.org/10.1023/A:1010933404324
