.. _clustering:

===================================================
Clustering
===================================================

`Clustering <http://en.wikipedia.org/wiki/Cluster_analysis>`_ of
unlabeled data can be performed with the module `scikits.learn.cluster`.

.. contents:: This page contents:

Affinity propagation
====================

.. figure:: ../auto_examples/cluster/images/plot_affinity_propagation.png
   :target: ../auto_examples/cluster/plot_affinity_propagation.html
   :align: center
   :scale: 50

.. autoclass:: scikits.learn.cluster.AffinityPropagation
    :members:

.. topic:: Examples:

 * :ref:`example_plot_affinity_propagation.py`: Affinity
   Propagation on a synthetic 2D datasets with 3 classes.

 * :ref:`example_stock_market.py` Affinity Propagation on Financial 
   time series to find groups of companies

Mean Shift
====================

.. figure:: ../auto_examples/cluster/images/plot_mean_shift.png
   :target: ../auto_examples/cluster/plot_mean_shift.html
   :align: center
   :scale: 50


.. autoclass:: scikits.learn.cluster.MeanShift
    :members:


.. topic:: Examples:

 * :ref:`example_plot_mean_shift.py`: Mean Shift clustering
   on a synthetic 2D datasets with 3 classes.

K-means
====================

.. autoclass:: scikits.learn.cluster.KMeans
    :members:

Spectral clustering
====================

Spectral clustering is especially efficient if the affinity matrix is
sparse.

.. autoclass:: scikits.learn.cluster.SpectralClustering
    :members:


.. topic:: Examples:

 * :ref:`example_cluster_plot_lena_segmentation.py`: Spectral clustering 
   to split the image of lena in regions.

 * :ref:`example_cluster_plot_segmentation_toy.py`: Segmenting objects
   from a noisy background using spectral clustering.

