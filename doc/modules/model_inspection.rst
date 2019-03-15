.. _model_inspection:

================
Model Inspection
================

.. currentmodule:: sklearn.model_inspection

A trained machine learning model can achieved good predictive performance in 
varies domains. In some cases, it is important for a user to interprete
the representations leared and decisions made by these models. The 
:mod:`sklearn.model_inspection` module provides tools to inspect an estimator 
after it it is trained. 

.. _permutation_importance:

Permutation feature importance
==============================

Feature importance is an interpretation technique that determines the features 
that are important to a model. The permutation importance of a feature is 
defined to be the decrease in a score when the feature's value is randomly 
shuffled [1]_. When the permutation importance is negative, the model improved when
the feature is permuted. This suggests that the model would beneift from the
removing the feature. This technique benefits from being model agnostic and 
only needing to fit the model once. 

The :func:`permutation_importance` function calcuates the feature importance
of estimators based on permutation. since the features are randomly shuffled,
this technique naturally extends itself to being bootstrapped. The 
``n_bootstrap`` paramter controls how many times a feature is permuted. Here
the permutation importance distributions for features in the
titantic dataset are computed:

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_permutation_importance_002.png
   :target: ../auto_examples/model_inspection/plot_permutation_importance.html
   :align: center
   :scale: 50%

.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_neighbors_plot_permutation_importance.py`: 
     Permutation Importance vs Random Forest Feature Importance

.. topic:: References:

   .. [1] Breiman, L. Machine Learning (2001) 45: 5.
        https://doi.org/10.1023/A:1010933404324
