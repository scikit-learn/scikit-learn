
.. _permutation_importance:

Permutation feature importance
==============================

Feature importance is an interpretation technique that determines the features 
that are important to a model. The permutation importance of a feature is 
defined to be the decrease in a score when the feature's value is randomly 
shuffled [1]_. When the permutation importance is negative, the model improves 
when the feature is permuted. This suggests that the model would beneift from 
the removing the feature. This technique benefits from being model agnostic and 
only needing to fit the model once. 

The :func:`permutation_importance` function calcuates the feature importance
of estimators based on permutation. since the features are randomly shuffled,
this procedure can be repeated to obtain a distribution of importances. The 
``n_rounds`` paramter controls how many times a feature is permuted. Here
the permutation importance distributions for features in the
titantic dataset are computed:

.. figure:: ../auto_examples/inspection/images/sphx_glr_plot_permutation_importance_002.png
   :target: ../auto_examples/inspection/plot_permutation_importance.html
   :align: center
   :scale: 100

.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`: 
     Permutation Importance vs Random Forest Feature Importance

.. topic:: References:

   .. [1] Breiman, L. Machine Learning (2001) 45: 5.
     https://doi.org/10.1023/A:1010933404324
