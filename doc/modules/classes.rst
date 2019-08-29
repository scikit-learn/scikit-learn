.. _api_ref:

=============
API Reference
=============

This is the class and function reference of scikit-learn. Please refer to
the :ref:`full user guide <user_guide>` for further details, as the class and
function raw specifications may not be enough to give full guidelines on their
uses.
For reference on concepts repeated across the API, see :ref:`glossary`.


:mod:`sklearn.base`: Base classes and utility functions
=======================================================

.. automodule:: sklearn.base
    :no-members:
    :no-inherited-members:

Base classes
------------
.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   base.BaseEstimator
   base.BiclusterMixin
   base.ClassifierMixin
   base.ClusterMixin
   base.DensityMixin
   base.RegressorMixin
   base.TransformerMixin

Functions
---------
.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   base.clone
   base.is_classifier
   base.is_regressor
   config_context
   get_config
   set_config
   show_versions

.. _calibration_ref:

:mod:`sklearn.calibration`: Probability Calibration
===================================================

.. automodule:: sklearn.calibration
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`calibration` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   calibration.CalibratedClassifierCV


.. autosummary::
   :toctree: generated/
   :template: function.rst

   calibration.calibration_curve

.. _cluster_ref:

:mod:`sklearn.cluster`: Clustering
==================================

.. automodule:: sklearn.cluster
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`clustering` section for further details.

Classes
-------
.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   cluster.AffinityPropagation
   cluster.AgglomerativeClustering
   cluster.Birch
   cluster.DBSCAN
   cluster.FeatureAgglomeration
   cluster.KMeans
   cluster.MiniBatchKMeans
   cluster.MeanShift
   cluster.OPTICS
   cluster.SpectralClustering

Functions
---------
.. autosummary::
   :toctree: generated/
   :template: function.rst

   cluster.affinity_propagation
   cluster.cluster_optics_dbscan
   cluster.cluster_optics_xi
   cluster.compute_optics_graph
   cluster.dbscan
   cluster.estimate_bandwidth
   cluster.k_means
   cluster.mean_shift
   cluster.spectral_clustering
   cluster.ward_tree

.. _bicluster_ref:

:mod:`sklearn.cluster.bicluster`: Biclustering
==============================================

.. automodule:: sklearn.cluster.bicluster
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`biclustering` section for further details.

Classes
-------
.. currentmodule:: sklearn.cluster.bicluster

.. autosummary::
   :toctree: generated/
   :template: class.rst

   SpectralBiclustering
   SpectralCoclustering

.. _compose_ref:

:mod:`sklearn.compose`: Composite Estimators
============================================

.. automodule:: sklearn.compose
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`combining_estimators` section for further
details.

.. currentmodule:: sklearn

.. autosummary::
    :toctree: generated
    :template: class.rst

    compose.ColumnTransformer
    compose.TransformedTargetRegressor

.. autosummary::
   :toctree: generated/
   :template: function.rst

   compose.make_column_transformer
