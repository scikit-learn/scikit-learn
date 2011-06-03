.. _clustering:

==========
Clustering
==========

`Clustering <http://en.wikipedia.org/wiki/Cluster_analysis>`__ of
unlabeled data can be performed with the module :mod:`scikits.learn.cluster`.

Each clustering algorithm comes in two variants: a class, that implements
the `fit` method to learn the clusters on train data, and a function,
that, given train data, returns an array of integer labels corresponding
to the different clusters. For the class, the labels over the training
data can be found in the `labels_` attribute.

.. currentmodule:: scikits.learn.cluster

.. topic:: Input data

    One important thing to note is that the algorithms implemented in
    this module take different kinds of matrix as input.  On one hand,
    :class:`MeanShift` and :class:`KMeans` take data matrices of shape
    [n_samples, n_features]. These can be obtained from the classes in
    the :mod:`scikits.learn.feature_extraction` module. On the other hand,
    :class:`AffinityPropagation` and :class:`SpectralClustering` take
    similarity matrices of shape [n_samples, n_samples].  These can be
    obtained from the functions in the :mod:`scikits.learn.metrics.pairwise`
    module. In other words, :class:`MeanShift` and :class:`KMeans` work
    with points in a vector space, whereas :class:`AffinityPropagation`
    and :class:`SpectralClustering` can work with arbitrary objects, as
    long as a similarity measure exists for such objects.


.. _k_means:

K-means
=======

The :class:`KMeans` algorithm clusters data by trying to separate samples
in n groups of equal variance, minimizing a criterion known as the
'inertia' of the groups. This algorithm requires the number of cluster to
be specified. It scales well to large number of samples, however its
results may be dependent on an initialisation. As a result, the computation is
often done several times, with different initialisation of the centroids.

.. topic:: Examples:

 * :ref:`example_cluster_kmeans_digits.py`: Clustering handwritten digits


.. _mini_batch_kmeans:

Mini Batch K-Means
-------------------

The :class:`MiniBatchKMeans` is a variant of the :class:`K-Means` algorithm
using mini-batches, random subset of the dataset, to compute the centroids.

Althought the :class:`MiniBatchKMeans` converge faster than the KMeans 
version, the quality of the results, measured by the inertia, the sum of 
the distance of each points to the nearest centroid, is not as good as 
the :class:`KMeans` algorithm.

.. figure:: ../auto_examples/cluster/images/plot_mini_batch_kmeans_1.png
   :target: ../auto_examples/cluster/plot_mini_batch_kmeans.html
   :align: center
   :scale: 100


.. topic:: Examples:

 * :ref:`example_cluster_plot_mini_batch_kmeans.py`: Comparison of the KMeans and
   BatchKMeans


.. topic:: References:

 * `"Web Scale K-Means clustering"
   <http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf>`_
   D. Sculley, *Proceedings of the 19th international conference on World 
   wide web* (2010)


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

 * :ref:`example_applications_stock_market.py` Affinity Propagation on
   Financial time series to find groups of companies


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


Spectral clustering
===================

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

.. topic:: References:

 * `"A Tutorial on Spectral Clustering"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.165.9323>`_
   Ulrike von Luxburg, 2007

 * `"Normalized cuts and image segmentation"
   <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.160.2324>`_
   Jianbo Shi, Jitendra Malik, 2000

 * `"A Random Walks View of Spectral Segmentation"
   <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.33.1501>`_
   Marina Meila, Jianbo Shi, 2001

 * `"On Spectral Clustering: Analysis and an algorithm"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.19.8100>`_
   Andrew Y. Ng, Michael I. Jordan, Yair Weiss, 2001


.. _hierarchical_clustering:

Hierarchical clustering
=======================

Hierarchical clustering is a general family of clustering algorithms that
build nested clusters by merging them successively. This hierarchy of
clusters represented as a tree (or dendrogram). The root of the tree is
the unique cluster that gathers all the samples, the leaves being the
clusters with only one sample. See the `Wikipedia page
<http://en.wikipedia.org/wiki/Hierarchical_clustering>`_ for more
details.

The :class:`Ward` object performs a hierarchical clustering based on
the Ward algorithm, that is a variance-minimizing approach. At each
step, it minimizes the sum of squared differences within all clusters
(inertia criterion).

This algorithm can scale to large number of samples when it is used jointly
with an connectivity matrix, but can be computationally expensive when no
connectivity constraints are added between samples: it considers at each step
all the possible merges.


Adding connectivity constraints
-------------------------------

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


Clustering performance evaluation
=================================

Evaluating the performance of a clustering algorithm is not as trivial as
counting the number of errors or the precision and recall of a supervised
classification algorithm. In particular any evaluation metric should not
take the absolute values of the cluster labels into account but rather
if this clustering define separations of the data similar to some ground
truth set of classes or satisfying some assumption such that members
belong to the same class are more similar that members of different
classes according to some similarity metric.

.. currentmodule:: scikits.learn.metrics

Inertia
-------

Presentation and usage
~~~~~~~~~~~~~~~~~~~~~~

TODO: factorize inertia computation out of kmeans and then write me!


Advantages
~~~~~~~~~~

- No need for the ground truth knowledge of the "real" classes.

Drawbacks
~~~~~~~~~

- Inertia makes the assumption that clusters are convex and isotropic
  which is not always the case especially of the clusters are manifolds
  with weird shapes: for instance inertia is a useless metrics to evaluate
  clustering algorithm that tries to identify nested circles on a 2D plane.

- Inertia is not a normalized metrics: we just know that lower values are
  better and bounded by zero. One potential solution would be to adjust
  inertia for random clustering (assuming the number of ground truth classes
  is known).


Homogeneity, completeness and V-measure
---------------------------------------

Presentation and usage
~~~~~~~~~~~~~~~~~~~~~~

Given the knowledge of the ground truth class assignments of the samples,
it is possible to define some intuitive metric using conditional entropy
analysis.

In particular Rosenberg and Hirschberg (2007) define the following two
desirable objectives for any cluster assignment:

- **homogeneity**: each cluster contains only members of a single class.

- **completeness**: all members of a given class are assigned to the same
  cluster.

We can turn those concept as scores :func:`homogeneity_score` and
:func:`completeness_score`. Both are bounded below by 0.0 and above by
1.0 (higher is better)::

  >>> from scikits.learn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]

  >>> metrics.homogeneity_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.66...

  >>> metrics.completeness_score(labels_true, labels_pred) # doctest: +ELLIPSIS
  0.42...

Their harmonic mean called **V-measure** is computed by
:func:`v_measure_score`::

  >>> metrics.v_measure_score(labels_true, labels_pred)    # doctest: +ELLIPSIS
  0.51...

All three metrics can be computed at once using
:func:`homogeneity_completeness_v_measure` as follows::

  >>> metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
  ...                                                      # doctest: +ELLIPSIS
  (0.66..., 0.42..., 0.51...)

The following clustering assignment is slighlty better, since it is
homogeneous but not complete::

  >>> labels_pred = [0, 0, 0, 1, 2, 2]
  >>> metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
  ...                                                      # doctest: +ELLIPSIS
  (1.0, 0.68..., 0.81...)

.. note::

  :func:`v_measure_score` is **symmetric**: it can be used to evaluate
  the **agreement** of two independent assignements on the same dataset.

  This is not the case for :func:`completeness_score` and
  :func:`homogeneity_score`: both are bound by the relationship::

    homogeneity_score(a, b) == completeness_score(b, a)


Mathematical formulation
------------------------

Homogeneity and completeness scores are formally given by:

.. math:: h = 1 - \frac{H(C|K)}{H(C)}

.. math:: c = 1 - \frac{H(K|C)}{H(K)}

where :math:`H(C|K)` is the **conditional entropy of the classes given
the cluster assignments** and is given by:

.. math:: H(C|K) = - \sum_{c=1}^{|C|} \sum_{k=1}^{|K|} \frac{n_{c,k}}{n}
          \cdot log(\frac{n_{c,k}}{n_k})

and :math:`H(C)` is the **entropy of the classes** and is given by:

.. math:: H(C) = - \sum_{c=1}^{|C|} \frac{n_c}{n} \cdot log(\frac{n_c}{n})

with :math:`n` the total number of samples, :math:`n_c` and :math:`n_k`
the number of samples respectively belonging to class :math:`c` and
cluster :math:`k`, and finally :math:`n_{c,k}` the number of samples
from class :math:`c` assigned to cluster :math:`k`.

The **conditional entropy of clusters given class** :math:`H(K|C)` and the
**entropy of clusters** :math:`H(K)` are defined in a symmetric manner.

Rosenberg and Hirschberg further define **V-measure** as the **harmonic
mean of homogeneity and completeness**:

.. math:: v = 2 \cdot \frac{h \cdot c}{h + c}

.. topic:: References

 * `"V-Measure: A conditional entropy-based external cluster evaluation
   measure" <http://acl.ldc.upenn.edu/D/D07/D07-1043.pdf>`_
   Andrew Rosenberg and Julia Hirschberg, 2007


Advantages
~~~~~~~~~~

- Bounded scores: 0.0 is as bad as it can be, 1.0 is a perfect score

- Intuitive interpretation: clustering with bad V-measure can be
  qualitatively analyzed in terms of homogeneity and completeness to
  better feel what 'kind' of mistakes is done by the assigmenent.

- No assumption is made on the similarity metric and the cluster
  structure.


Drawbacks
~~~~~~~~~

- These metrics require the knowlege of the ground truth classes while
  almost never available in practice or requires manual assignment by
  human annotators (as in the supervised learning setting).

- The previously introduced metrics are not normalized w.r.t. random
  labeling: this means that depending on the number of samples,
  clusters and ground truth classes, a completely random labeling will
  not always yield the same values for homogeneity, completeness and
  hence v-measure. In particular random labeling won't yield zero scores.

  TODO: check the values we get for random labeling on various problem
  sizes to know whether this is a real problem in practice.
