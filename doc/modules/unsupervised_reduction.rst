
.. _data_reduction:

=====================================
Unsupervised dimensionality reduction
=====================================

If your number of features is high, it may be useful to reduce it with an
unsupervised step prior to supervised steps. Many of the
:ref:`unsupervised-learning` methods implement a ``transform`` method that
can be used to reduce the dimensionality. Below we discuss two specific
example of this pattern that are heavily used.

.. topic:: **Pipelining**

    The unsupervised data reduction and the supervised estimator can be
    chained in one step. See :ref:`pipeline`.

.. currentmodule:: sklearn

PCA: principal component analysis
----------------------------------

:class:`decomposition.PCA` looks for a combination of features that
capture well the variance of the original features. See :ref:`decompositions`.

.. topic:: **Examples**

   * :ref:`example_applications_face_recognition.py`

Random projections
-------------------

The module: :mod:`random_projection` provides several tools for data
reduction by random projections. See the relevant section of the
documentation: :ref:`random_projection`.

.. topic:: **Examples**

   * :ref:`example_plot_johnson_lindenstrauss_bound.py`

Feature agglomeration
------------------------

:class:`cluster.FeatureAgglomeration` applies
:ref:`hierarchical_clustering` to group together features that behave
similarly.

.. topic:: **Examples**

   * :ref:`example_cluster_plot_feature_agglomeration_vs_univariate_selection.py`
   * :ref:`example_cluster_plot_digits_agglomeration.py`

.. topic:: **Feature scaling**

   Note that if features have very different scaling or statistical
   properties, :class:`cluster.FeatureAgglomeration` may not be able to
   capture the links between related features. Using a 
   :class:`preprocessing.StandardScaler` can be useful in these settings.

