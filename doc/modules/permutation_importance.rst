
.. _permutation_importance:

Permutation feature importance
==============================

Feature importance is an model interpretation technique that determines the features that are important to a model. The *permutation* feature importance is defined to be the decrease in a model's score when the feature's value is randomly shuffled [1]_. This procedure breaks the relationship between the feature and the target, thus the drop in the model's score is analogous to how much the model depends on the feature. This technique benefits from being model agnostic, only requiring the model to be `fitted` once, and can be calcuated many times with different permutations.

The :func:`permutation_importance` function calculates the feature importance of `estimators` for a given dataset. The ``n_rounds`` paramter sets the number of times a feature is randomly shuffled and returns a distribution of feature importances. Note that :func:`permutation_importance` can accept a training dataset or test dataset for computing importances. Feature importance based on the training set can inflate the importances, when the in fact, the model is overfitting and the feature is not important. Feature importance based on the test set gives the importance of a feature on unseen data, which gives a sense of the features that are *actually* important.

The :class:`~sklearn.ensemble.RandomForestClassifer` provides its own tree based feature importance that computes statistcs dervived from the training dataset. This can give importances to features that are not predictive of the target on unseen data. This use case is explored in: ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`.

When a correlated feature is present in a dataset it can decrease the importance of that feature. For example, in the extreme case where two features are exactly the same, when one feature is permuted, the model still has access to the feature through its correlated feature. This will result in a lower importance for both features, where they might *actually* important. One way to handle this is to cluster features that are correlated and only keep one feature from each cluster. This use case is explored in: :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`.

.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`: 
     Permutation Importance vs Random Forest Feature Importance
  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`: 
     Permutation Importance with Multicollinear Features

.. topic:: References:

   .. [1] Breiman, L. Machine Learning (2001) 45: 5.
     https://doi.org/10.1023/A:1010933404324
