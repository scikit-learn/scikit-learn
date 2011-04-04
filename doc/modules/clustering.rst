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

.. figure:: ../auto_examples/cluster/images/plot_affinity_propagation_1.png
   :target: ../auto_examples/cluster/plot_affinity_propagation.html
   :align: center
   :scale: 50

.. topic:: Examples:

 * :ref:`example_cluster_plot_affinity_propagation.py`: Affinity
   Propagation on a synthetic 2D datasets with 3 classes.

 * :ref:`example_applications_stock_market.py` Affinity Propagation on Financial
   time series to find groups of companies


Mean Shift
==========

:class:`MeanShift` clusters data by estimating *blobs* in a smooth
density of points matrix. This algorithm automatically sets its numbers
of cluster. It will have difficulties scaling to thousands of samples.


.. figure:: ../auto_examples/cluster/images/plot_mean_shift_1.png
   :target: ../auto_examples/cluster/plot_mean_shift.html
   :align: center
   :scale: 50


.. topic:: Examples:

 * :ref:`example_cluster_plot_mean_shift.py`: Mean Shift clustering
   on a synthetic 2D datasets with 3 classes.


K-means
=======

The :class:`KMeans` algorithm clusters data by trying to separate samples
in n groups of equal variance, minimizing a criterion known as the
'inertia' of the groups. This algorithm requires the number of cluster to
be specified. It scales well to large number of samples, however its
results may be dependent on an initialisation.


Spectral clustering
====================

:class:`SpectralClustering` does a low-dimension embedding of the
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


.. |noisy_img| image:: ../auto_examples/cluster/images/plot_segmentation_toy_1.png
    :target: ../auto_examples/cluster/plot_segmentation_toy.html
    :scale: 50

.. |segmented_img| image:: ../auto_examples/cluster/images/plot_segmentation_toy_2.png
    :target: ../auto_examples/cluster/plot_segmentation_toy.html
    :scale: 50

.. centered:: |noisy_img| |segmented_img|

.. topic:: Examples:

 * :ref:`example_cluster_plot_segmentation_toy.py`: Segmenting objects
   from a noisy background using spectral clustering.

 * :ref:`example_cluster_plot_lena_segmentation.py`: Spectral clustering
   to split the image of lena in regions.


.. _hierarchical_clustering:

Hierarchical clustering
=======================

Hierarchical clustering is a general family of clustering algorithms that
build nested clusters by merging them successively. This hierarchy of
clusters represented as a tree (or dendrogram). The root of the tree is
the unique cluster that gathers all the samples, the leaves being the
clusters with only one sample. See the `Wikipedia page
<http://en.wikipedia.org/wiki/Hierarchical_clustering for more
details>`_.


The :class:`Ward` object performs a hierarchical clustering based on Ward
algorithm, that is a variance-minimizing approach. At each step, it
minimizes the sum of squared differences within all clusters (inertia
criterion). 

This algorithm can scale to large number of samples when it is used jointly
with an connectivity matrix, but can be computationally expensive when no
connectivity constraints are added between samples: it considers at each step
all the possible merges.

Adding connectivity constraints 
----------------------------------

An interesting aspect of the :class:`Ward` object is that connectivity
constraints can be added to this algorithm (only adjacent clusters can be
merged together), through an connectivity matrix that defines for each
sample the neighboring samples following a given structure of the data. For 
instance, in the swiss-roll example below, the connectivity constraints
forbid the merging of points that are not adjacent on the swiss roll, and
thus avoid forming clusters that extend across overlapping folds of the
roll.

.. |unstructured| image:: ../auto_examples/cluster/images/plot_ward_structured_vs_unstructured_1.png
        :target: ../auto_examples/cluster/plot_ward_structured_vs_unstructured.html
        :scale: 50

.. |structured| image:: ../auto_examples/cluster/images/plot_ward_structured_vs_unstructured_2.png
        :target: ../auto_examples/cluster/plot_ward_structured_vs_unstructured.html
        :scale: 50

.. centered:: |unstructured| |structured|


The connectivity constraints are imposed via an connectivity matrix: a
scipy sparse matrix that has elements only at the intersection of a row
and a column with indices of the dataset that should be connected. This
matrix can be constructed from apriori information, for instance if you
whish to cluster web pages, but only merging pages with a link pointing
from one to another. It can also be learned from the data, for instance
using :func:`scikits.learn.neighbors.kneighbors_graph` to restrict
merging to nearest neighbors as in the :ref:`swiss roll
<example_cluster_plot_ward_structured_vs_unstructured.py>` example, or
using :func:`scikits.learn.feature_extraction.image.grid_to_graph` to
enable only merging of neighboring pixels on an image, as in the
:ref:`Lena <example_cluster_plot_lena_ward_segmentation.py>` example.

.. topic:: Examples:

 * :ref:`example_cluster_plot_lena_ward_segmentation.py`: Ward clustering 
   to split the image of lena in regions.

 * :ref:`example_cluster_plot_ward_structured_vs_unstructured.py`: Example of 
    Ward algorithm on a swiss-roll, comparison of structured approaches
    versus unstructured approaches.

 * :ref:`example_cluster_plot_feature_agglomeration_vs_univariate_selection.py`:
    Example of dimensionality reduction with feature agglomeration based on
    Ward hierarchical clustering.


