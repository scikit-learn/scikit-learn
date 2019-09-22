
.. _permutation_importance:

Permutation feature importance
==============================

.. currentmodule:: sklearn.inspection

Permutation feature importance is a model inspection technique that can be used
for any :term:`fitted` :term:`estimator` when the data is rectangular. This is
especially useful for non-linear or opaque :term:`estimators`. The permutation
feature importance is defined to be the decrease in a model score when a single
feature value is randomly shuffled [1]_. This procedure breaks the relationship
between the feature and the target, thus the drop in the model score is
indicative of how much the model depends on the feature. This technique
benefits from being model agnostic and can be calculated many times with
different permutations of the feature.

The :func:`permutation_importance` function calculates the feature importance
of :term:`estimators` for a given dataset. The ``n_repeats`` parameter sets the
number of times a feature is randomly shuffled and returns a sample of feature
importances. Permutation importances can either be computed on the training set
or an held-out testing or validation set. Using a held-out set makes it
possible to highlight which features contribute the most to the generalization
power of the inspected model. Features that are important on the training set
but not on the held-out set might cause the model to overfit.

Note that features that are deemed non-important for some model with a
low predictive performance could be highly predictive for a model that
generalizes better. The conclusions should always be drawn in the context of
the specific model under inspection and cannot be automatically generalized to
the intrinsic predictive value of the features by them-selves. Therefore it is
always important to evaluate the predictive power of a model using a held-out
set (or better with cross-validation) prior to computing importances.

Relation to impurity-based importance in trees
----------------------------------------------

Tree based models provides a different measure of feature importances based
on the mean decrease in impurity (MDI, the splitting criterion). This gives
importance to features that may not be predictive on unseen data. The
permutation feature importance avoids this issue, since it can be applied to
unseen data. Furthermore, impurity-based feature importance for trees
are strongly biased and favor high cardinality features
(typically numerical features). Permutation-based feature importances do not
exhibit such a bias. Additionally, the permutation feature importance may use
an arbitrary metric on the tree's predictions. These two methods of obtaining
feature importance are explored in:
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
