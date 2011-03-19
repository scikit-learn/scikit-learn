.. _clustering:

===================================================
Clustering
===================================================

`Clustering <http://en.wikipedia.org/wiki/Cluster_analysis>`__ of
unlabeled data can be performed with the module `scikits.learn.cluster`.

Each clustering algorithm comes in two variants: a class, that implements
the `fit` method to learn the clusters on train data, and a function,
that, given train data, returns an array of integer labels corresponding
to the different clusters. For the class, the labels over the training
data can be found in the `labels_` attribute.

.. currentmodule:: scikits.learn.cluster

One important thing to note is that the algorithms implemented in this module
take different kinds of matrix as input.  On one hand, :class:`MeanShift` and
:class:`KMeans` take data matrices of shape [n_samples, n_features]. These can
be obtained from the classes in the `scikits.learn.feature_extraction` module.
On the other hand, :class:`AffinityPropagation` and :class:`SpectralClustering`
take similarity matrices of shape [n_samples, n_samples].  These can be
obtained from the functions in the `scikits.learn.metrics.pairwise` module.
In other words, :class:`MeanShift` and :class:`KMeans` work with points in a
vector space, whereas :class:`AffinityPropagation` and
:class:`SpectralClustering` can work with arbitrary objects, as long as a
similarity measure exists for such objects.


Affinity propagation
====================


:class:`AffinityPropagation` clusters data by diffusion in the similarity
matrix. This algorithm automatically sets its numbers of cluster. It
will have difficulties scaling to thousands of samples.

.. figure:: ../auto_examples/cluster/images/plot_affinity_propagation.png
   :target: ../auto_examples/cluster/plot_affinity_propagation.html
   :align: center
   :scale: 50

.. topic:: Examples:

 * :ref:`example_cluster_plot_affinity_propagation.py`: Affinity
   Propagation on a synthetic 2D datasets with 3 classes.

 * :ref:`example_applications_stock_market.py` Affinity Propagation on Financial 
   time series to find groups of companies

Mean Shift
====================

:class:`MeanShift` clusters data by estimating *blobs* in a smooth
density of points matrix. This algorithm automatically sets its numbers
of cluster. It will have difficulties scaling to thousands of samples.


.. figure:: ../auto_examples/cluster/images/plot_mean_shift.png
   :target: ../auto_examples/cluster/plot_mean_shift.html
   :align: center
   :scale: 50


.. topic:: Examples:

 * :ref:`example_cluster_plot_mean_shift.py`: Mean Shift clustering
   on a synthetic 2D datasets with 3 classes.

K-means
====================

The :class:`KMeans` algorithm clusters data by trying to separate samples
in n groups of equal variance, minimizing a criterion known as the
'inertia' of the groups. This algorithm requires the number of cluster to
be specified. It scales well to large number of samples, however its
results may be dependent on an initialisation.


Spectral clustering
====================

:class:`SpectralClustering` does an low-dimension embedding of the
affinity matrix between samples, followed by a KMeans in the low
dimensional space. It is especially efficient if the affinity matrix is
sparse and the `pyamg <http://code.google.com/p/pyamg/>`_ module is
installed. SpectralClustering requires the number of clusters to be
specified. It works well for a small number of clusters but is not
advised when using many clusters. 

For two clusters, it solves a convex relaxation of the `normalised
cuts <http://www.cs.berkeley.edu/~malik/papers/SM-ncut.pdf>`_ problem on
the similarity graph: cutting the graph in two so that the weight of the
edges cut is small compared to the weights in of edges inside each
cluster. This criteria is especially interesting when working on images:
graph vertices are pixels, and edges of the similarity graph are a
function of the gradient of the image.

.. figure:: ../auto_examples/cluster/images/plot_segmentation_toy.png
   :target: ../auto_examples/cluster/plot_segmentation_toy.html
   :align: center
   :scale: 50



.. topic:: Examples:

 * :ref:`example_cluster_plot_lena_segmentation.py`: Spectral clustering 
   to split the image of lena in regions.

 * :ref:`example_cluster_plot_segmentation_toy.py`: Segmenting objects
   from a noisy background using spectral clustering.

