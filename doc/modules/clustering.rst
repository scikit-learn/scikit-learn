.. _clustering:

==========
Clustering
==========

`Clustering <https://en.wikipedia.org/wiki/Cluster_analysis>`__ of
unlabeled data can be performed with the module :mod:`sklearn.cluster`.

Each clustering algorithm comes in two variants: a class, that implements
the ``fit`` method to learn the clusters on train data, and a function,
that, given train data, returns an array of integer labels corresponding
to the different clusters. For the class, the labels over the training
data can be found in the ``labels_`` attribute.

.. currentmodule:: sklearn.cluster

.. topic:: Input data

    One important thing to note is that the algorithms implemented in
    this module can take different kinds of matrix as input. All the
    methods accept standard data matrices of shape ``(n_samples, n_features)``.
    These can be obtained from the classes in the :mod:`sklearn.feature_extraction`
    module. For :class:`AffinityPropagation`, :class:`SpectralClustering`
    and :class:`DBSCAN` one can also input similarity matrices of shape
    ``(n_samples, n_samples)``. These can be obtained from the functions
    in the :mod:`sklearn.metrics.pairwise` module.

Overview of clustering methods
===============================

.. figure:: ../auto_examples/cluster/images/sphx_glr_plot_cluster_comparison_001.png
   :target: ../auto_examples/cluster/plot_cluster_comparison.html
   :align: center
   :scale: 50

   A comparison of the clustering algorithms in scikit-learn


.. list-table::
   :header-rows: 1
   :widths: 14 15 19 25 20

   * - Method name
     - Parameters
     - Scalability
     - Usecase
     - Geometry (metric used)

   * - :ref:`K-Means <k_means>`
     - number of clusters
     - Very large ``n_samples``, medium ``n_clusters`` with
       :ref:`MiniBatch code <mini_batch_kmeans>`
     - General-purpose, even cluster size, flat geometry,
       not too many clusters, inductive
     - Distances between points

   * - :ref:`Affinity propagation <affinity_propagation>`
     - damping, sample preference
     - Not scalable with n_samples
     - Many clusters, uneven cluster size, non-flat geometry, inductive
     - Graph distance (e.g. nearest-neighbor graph)

   * - :ref:`Mean-shift <mean_shift>`
     - bandwidth
     - Not scalable with ``n_samples``
     - Many clusters, uneven cluster size, non-flat geometry, inductive
     - Distances between points

   * - :ref:`Spectral clustering <spectral_clustering>`
     - number of clusters
     - Medium ``n_samples``, small ``n_clusters``
     - Few clusters, even cluster size, non-flat geometry, transductive
     - Graph distance (e.g. nearest-neighbor graph)

   * - :ref:`Ward hierarchical clustering <hierarchical_clustering>`
     - number of clusters or distance threshold
     - Large ``n_samples`` and ``n_clusters``
     - Many clusters, possibly connectivity constraints, transductive
     - Distances between points

   * - :ref:`Agglomerative clustering <hierarchical_clustering>`
     - number of clusters or distance threshold, linkage type, distance
     - Large ``n_samples`` and ``n_clusters``
     - Many clusters, possibly connectivity constraints, non Euclidean
       distances, transductive
     - Any pairwise distance

   * - :ref:`DBSCAN <dbscan>`
     - neighborhood size
     - Very large ``n_samples``, medium ``n_clusters``
     - Non-flat geometry, uneven cluster sizes, outlier removal,
       transductive
     - Distances between nearest points

   * - :ref:`HDBSCAN <hdbscan>`
     - minimum cluster membership, minimum point neighbors
     - large ``n_samples``, medium ``n_clusters``
     - Non-flat geometry, uneven cluster sizes, outlier removal,
       transductive, hierarchical, variable cluster density
     - Distances between nearest points

   * - :ref:`OPTICS <optics>`
     - minimum cluster membership
     - Very large ``n_samples``, large ``n_clusters``
     - Non-flat geometry, uneven cluster sizes, variable cluster density,
       outlier removal, transductive
     - Distances between points

   * - :ref:`Gaussian mixtures <mixture>`
     - many
     - Not scalable
     - Flat geometry, good for density estimation, inductive
     - Mahalanobis distances to  centers

   * - :ref:`BIRCH <birch>`
     - branching factor, threshold, optional global clusterer.
     - Large ``n_clusters`` and ``n_samples``
     - Large dataset, outlier removal, data reduction, inductive
     - Euclidean distance between points

   * - :ref:`Bisecting K-Means <bisect_k_means>`
     - number of clusters
     - Very large ``n_samples``, medium ``n_clusters``
     - General-purpose, even cluster size, flat geometry,
       no empty clusters, inductive, hierarchical
     - Distances between points

Non-flat geometry clustering is useful when the clusters have a specific
shape, i.e. a non-flat manifold, and the standard euclidean distance is
not the right metric. This case arises in the two top rows of the figure
above.

Gaussian mixture models, useful for clustering, are described in
:ref:`another chapter of the documentation <mixture>` dedicated to
mixture models. KMeans can be seen as a special case of Gaussian mixture
model with equal covariance per component.

:term:`Transductive <transductive>` clustering methods (in contrast to
:term:`inductive` clustering methods) are not designed to be applied to new,
unseen data.

.. _k_means:

K-means
=======

The :class:`KMeans` algorithm clusters data by trying to separate samples in n
groups of equal variance, minimizing a criterion known as the *inertia* or
within-cluster sum-of-squares (see below). This algorithm requires the number
of clusters to be specified. It scales well to large numbers of samples and has
been used across a large range of application areas in many different fields.

The k-means algorithm divides a set of :math:`N` samples :math:`X` into
:math:`K` disjoint clusters :math:`C`, each described by the mean :math:`\mu_j`
of the samples in the cluster. The means are commonly called the cluster
"centroids"; note that they are not, in general, points from :math:`X`,
although they live in the same space.

The K-means algorithm aims to choose centroids that minimise the **inertia**,
or **within-cluster sum-of-squares criterion**:

.. math:: \sum_{i=0}^{n}\min_{\mu_j \in C}(||x_i - \mu_j||^2)

Inertia can be recognized as a measure of how internally coherent clusters are.
It suffers from various drawbacks:

- Inertia makes the assumption that clusters are convex and isotropic,
  which is not always the case. It responds poorly to elongated clusters,
  or manifolds with irregular shapes.

- Inertia is not a normalized metric: we just know that lower values are
  better and zero is optimal. But in very high-dimensional spaces, Euclidean
  distances tend to become inflated
  (this is an instance of the so-called "curse of dimensionality").
  Running a dimensionality reduction algorithm such as :ref:`PCA` prior to
  k-means clustering can alleviate this problem and speed up the
  computations.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_kmeans_assumptions_002.png
   :target: ../auto_examples/cluster/plot_kmeans_assumptions.html
   :align: center
   :scale: 50

K-means is often referred to as Lloyd's algorithm. In basic terms, the
algorithm has three steps. The first step chooses the initial centroids, with
the most basic method being to choose :math:`k` samples from the dataset
:math:`X`. After initialization, K-means consists of looping between the
two other steps. The first step assigns each sample to its nearest centroid.
The second step creates new centroids by taking the mean value of all of the
samples assigned to each previous centroid. The difference between the old
and the new centroids are computed and the algorithm repeats these last two
steps until this value is less than a threshold. In other words, it repeats
until the centroids do not move significantly.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_kmeans_digits_001.png
   :target: ../auto_examples/cluster/plot_kmeans_digits.html
   :align: right
   :scale: 35

K-means is equivalent to the expectation-maximization algorithm
with a small, all-equal, diagonal covariance matrix.

The algorithm can also be understood through the concept of `Voronoi diagrams
<https://en.wikipedia.org/wiki/Voronoi_diagram>`_. First the Voronoi diagram of
the points is calculated using the current centroids. Each segment in the
Voronoi diagram becomes a separate cluster. Secondly, the centroids are updated
to the mean of each segment. The algorithm then repeats this until a stopping
criterion is fulfilled. Usually, the algorithm stops when the relative decrease
in the objective function between iterations is less than the given tolerance
value. This is not the case in this implementation: iteration stops when
centroids move less than the tolerance.

Given enough time, K-means will always converge, however this may be to a local
minimum. This is highly dependent on the initialization of the centroids.
As a result, the computation is often done several times, with different
initializations of the centroids. One method to help address this issue is the
k-means++ initialization scheme, which has been implemented in scikit-learn
(use the ``init='k-means++'`` parameter). This initializes the centroids to be
(generally) distant from each other, leading to probably better results than
random initialization, as shown in the reference.

K-means++ can also be called independently to select seeds for other
clustering algorithms, see :func:`sklearn.cluster.kmeans_plusplus` for details
and example usage.

The algorithm supports sample weights, which can be given by a parameter
``sample_weight``. This allows to assign more weight to some samples when
computing cluster centers and values of inertia. For example, assigning a
weight of 2 to a sample is equivalent to adding a duplicate of that sample
to the dataset :math:`X`.

K-means can be used for vector quantization. This is achieved using the
transform method of a trained model of :class:`KMeans`.

Low-level parallelism
---------------------

:class:`KMeans` benefits from OpenMP based parallelism through Cython. Small
chunks of data (256 samples) are processed in parallel, which in addition
yields a low memory footprint. For more details on how to control the number of
threads, please refer to our :ref:`parallelism` notes.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_kmeans_assumptions.py`: Demonstrating when
   k-means performs intuitively and when it does not
 * :ref:`sphx_glr_auto_examples_cluster_plot_kmeans_digits.py`: Clustering handwritten digits

.. topic:: References:

 * `"k-means++: The advantages of careful seeding"
   <http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf>`_
   Arthur, David, and Sergei Vassilvitskii,
   *Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete
   algorithms*, Society for Industrial and Applied Mathematics (2007)

.. _mini_batch_kmeans:

Mini Batch K-Means
------------------

The :class:`MiniBatchKMeans` is a variant of the :class:`KMeans` algorithm
which uses mini-batches to reduce the computation time, while still attempting
to optimise the same objective function. Mini-batches are subsets of the input
data, randomly sampled in each training iteration. These mini-batches
drastically reduce the amount of computation required to converge to a local
solution. In contrast to other algorithms that reduce the convergence time of
k-means, mini-batch k-means produces results that are generally only slightly
worse than the standard algorithm.

The algorithm iterates between two major steps, similar to vanilla k-means.
In the first step, :math:`b` samples are drawn randomly from the dataset, to form
a mini-batch. These are then assigned to the nearest centroid. In the second
step, the centroids are updated. In contrast to k-means, this is done on a
per-sample basis. For each sample in the mini-batch, the assigned centroid
is updated by taking the streaming average of the sample and all previous
samples assigned to that centroid. This has the effect of decreasing the
rate of change for a centroid over time. These steps are performed until
convergence or a predetermined number of iterations is reached.

:class:`MiniBatchKMeans` converges faster than :class:`KMeans`, but the quality
of the results is reduced. In practice this difference in quality can be quite
small, as shown in the example and cited reference.

.. figure:: ../auto_examples/cluster/images/sphx_glr_plot_mini_batch_kmeans_001.png
   :target: ../auto_examples/cluster/plot_mini_batch_kmeans.html
   :align: center
   :scale: 100


.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_mini_batch_kmeans.py`: Comparison of KMeans and
   MiniBatchKMeans

 * :ref:`sphx_glr_auto_examples_text_plot_document_clustering.py`: Document clustering using sparse
   MiniBatchKMeans

 * :ref:`sphx_glr_auto_examples_cluster_plot_dict_face_patches.py`


.. topic:: References:

 * `"Web Scale K-Means clustering"
   <https://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf>`_
   D. Sculley, *Proceedings of the 19th international conference on World
   wide web* (2010)

.. _affinity_propagation:

Affinity Propagation
====================

:class:`AffinityPropagation` creates clusters by sending messages between
pairs of samples until convergence. A dataset is then described using a small
number of exemplars, which are identified as those most representative of other
samples. The messages sent between pairs represent the suitability for one
sample to be the exemplar of the other, which is updated in response to the
values from other pairs. This updating happens iteratively until convergence,
at which point the final exemplars are chosen, and hence the final clustering
is given.

.. figure:: ../auto_examples/cluster/images/sphx_glr_plot_affinity_propagation_001.png
   :target: ../auto_examples/cluster/plot_affinity_propagation.html
   :align: center
   :scale: 50


Affinity Propagation can be interesting as it chooses the number of
clusters based on the data provided. For this purpose, the two important
parameters are the *preference*, which controls how many exemplars are
used, and the *damping factor* which damps the responsibility and
availability messages to avoid numerical oscillations when updating these
messages.

The main drawback of Affinity Propagation is its complexity. The
algorithm has a time complexity of the order :math:`O(N^2 T)`, where :math:`N`
is the number of samples and :math:`T` is the number of iterations until
convergence. Further, the memory complexity is of the order
:math:`O(N^2)` if a dense similarity matrix is used, but reducible if a
sparse similarity matrix is used. This makes Affinity Propagation most
appropriate for small to medium sized datasets.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_affinity_propagation.py`: Affinity
   Propagation on a synthetic 2D datasets with 3 classes.

 * :ref:`sphx_glr_auto_examples_applications_plot_stock_market.py` Affinity Propagation on
   Financial time series to find groups of companies


**Algorithm description:**
The messages sent between points belong to one of two categories. The first is
the responsibility :math:`r(i, k)`,
which is the accumulated evidence that sample :math:`k`
should be the exemplar for sample :math:`i`.
The second is the availability :math:`a(i, k)`
which is the accumulated evidence that sample :math:`i`
should choose sample :math:`k` to be its exemplar,
and considers the values for all other samples that :math:`k` should
be an exemplar. In this way, exemplars are chosen by samples if they are (1)
similar enough to many samples and (2) chosen by many samples to be
representative of themselves.

More formally, the responsibility of a sample :math:`k`
to be the exemplar of sample :math:`i` is given by:

.. math::

    r(i, k) \leftarrow s(i, k) - max [ a(i, k') + s(i, k') \forall k' \neq k ]

Where :math:`s(i, k)` is the similarity between samples :math:`i` and :math:`k`.
The availability of sample :math:`k`
to be the exemplar of sample :math:`i` is given by:

.. math::

    a(i, k) \leftarrow min [0, r(k, k) + \sum_{i'~s.t.~i' \notin \{i, k\}}{r(i', k)}]

To begin with, all values for :math:`r` and :math:`a` are set to zero,
and the calculation of each iterates until convergence.
As discussed above, in order to avoid numerical oscillations when updating the
messages, the damping factor :math:`\lambda` is introduced to iteration process:

.. math:: r_{t+1}(i, k) = \lambda\cdot r_{t}(i, k) + (1-\lambda)\cdot r_{t+1}(i, k)
.. math:: a_{t+1}(i, k) = \lambda\cdot a_{t}(i, k) + (1-\lambda)\cdot a_{t+1}(i, k)

where :math:`t` indicates the iteration times.

.. _mean_shift:

Mean Shift
==========
:class:`MeanShift` clustering aims to discover *blobs* in a smooth density of
samples. It is a centroid based algorithm, which works by updating candidates
for centroids to be the mean of the points within a given region. These
candidates are then filtered in a post-processing stage to eliminate
near-duplicates to form the final set of centroids.

The position of centroid candidates is iteratively adjusted using a technique called hill
climbing, which finds local maxima of the estimated probability density.
Given a candidate centroid :math:`x` for iteration :math:`t`, the candidate
is updated according to the following equation:

.. math::

    x^{t+1} = x^t + m(x^t)

Where :math:`m` is the *mean shift* vector that is computed for each
centroid that points towards a region of the maximum increase in the density of points.
To compute :math:`m` we define :math:`N(x)` as the neighborhood of samples within
a given distance around :math:`x`. Then :math:`m` is computed using the following
equation, effectively updating a centroid to be the mean of the samples within
its neighborhood:

.. math::

    m(x) = \frac{1}{|N(x)|} \sum_{x_j \in N(x)}x_j - x

In general, the equation for :math:`m` depends on a kernel used for density estimation.
The generic formula is:

.. math::

    m(x) = \frac{\sum_{x_j \in N(x)}K(x_j - x)x_j}{\sum_{x_j \in N(x)}K(x_j - x)} - x

In our implementation, :math:`K(x)` is equal to 1 if :math:`x` is small enough and is
equal to 0 otherwise. Effectively :math:`K(y - x)` indicates whether :math:`y` is in
the neighborhood of :math:`x`.

The algorithm automatically sets the number of clusters, instead of relying on a
parameter ``bandwidth``, which dictates the size of the region to search through.
This parameter can be set manually, but can be estimated using the provided
``estimate_bandwidth`` function, which is called if the bandwidth is not set.

The algorithm is not highly scalable, as it requires multiple nearest neighbor
searches during the execution of the algorithm. The algorithm is guaranteed to
converge, however the algorithm will stop iterating when the change in centroids
is small.

Labelling a new sample is performed by finding the nearest centroid for a
given sample.


.. figure:: ../auto_examples/cluster/images/sphx_glr_plot_mean_shift_001.png
   :target: ../auto_examples/cluster/plot_mean_shift.html
   :align: center
   :scale: 50


.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_mean_shift.py`: Mean Shift clustering
   on a synthetic 2D datasets with 3 classes.

.. topic:: References:

 * :doi:`"Mean shift: A robust approach toward feature space analysis"
   <10.1109/34.1000236>`
   D. Comaniciu and P. Meer, *IEEE Transactions on Pattern Analysis and Machine Intelligence* (2002)


.. _spectral_clustering:

Spectral clustering
===================

:class:`SpectralClustering` performs a low-dimension embedding of the
affinity matrix between samples, followed by clustering, e.g., by KMeans,
of the components of the eigenvectors in the low dimensional space.
It is especially computationally efficient if the affinity matrix is sparse
and the `amg` solver is used for the eigenvalue problem (Note, the `amg` solver
requires that the `pyamg <https://github.com/pyamg/pyamg>`_ module is installed.)

The present version of SpectralClustering requires the number of clusters
to be specified in advance. It works well for a small number of clusters,
but is not advised for many clusters.

For two clusters, SpectralClustering solves a convex relaxation of the
`normalized cuts <https://people.eecs.berkeley.edu/~malik/papers/SM-ncut.pdf>`_
problem on the similarity graph: cutting the graph in two so that the weight of
the edges cut is small compared to the weights of the edges inside each
cluster. This criteria is especially interesting when working on images, where
graph vertices are pixels, and weights of the edges of the similarity graph are
computed using a function of a gradient of the image.


.. |noisy_img| image:: ../auto_examples/cluster/images/sphx_glr_plot_segmentation_toy_001.png
    :target: ../auto_examples/cluster/plot_segmentation_toy.html
    :scale: 50

.. |segmented_img| image:: ../auto_examples/cluster/images/sphx_glr_plot_segmentation_toy_002.png
    :target: ../auto_examples/cluster/plot_segmentation_toy.html
    :scale: 50

.. centered:: |noisy_img| |segmented_img|

.. warning:: Transforming distance to well-behaved similarities

    Note that if the values of your similarity matrix are not well
    distributed, e.g. with negative values or with a distance matrix
    rather than a similarity, the spectral problem will be singular and
    the problem not solvable. In which case it is advised to apply a
    transformation to the entries of the matrix. For instance, in the
    case of a signed distance matrix, is common to apply a heat kernel::

        similarity = np.exp(-beta * distance / distance.std())

    See the examples for such an application.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_segmentation_toy.py`: Segmenting objects
   from a noisy background using spectral clustering.

 * :ref:`sphx_glr_auto_examples_cluster_plot_coin_segmentation.py`: Spectral clustering
   to split the image of coins in regions.

.. |coin_kmeans| image:: ../auto_examples/cluster/images/sphx_glr_plot_coin_segmentation_001.png
    :target: ../auto_examples/cluster/plot_coin_segmentation.html
    :scale: 35

.. |coin_discretize| image:: ../auto_examples/cluster/images/sphx_glr_plot_coin_segmentation_002.png
    :target: ../auto_examples/cluster/plot_coin_segmentation.html
    :scale: 35

.. |coin_cluster_qr| image:: ../auto_examples/cluster/images/sphx_glr_plot_coin_segmentation_003.png
    :target: ../auto_examples/cluster/plot_coin_segmentation.html
    :scale: 35

Different label assignment strategies
-------------------------------------

Different label assignment strategies can be used, corresponding to the
``assign_labels`` parameter of :class:`SpectralClustering`.
``"kmeans"`` strategy can match finer details, but can be unstable.
In particular, unless you control the ``random_state``, it may not be
reproducible from run-to-run, as it depends on random initialization.
The alternative ``"discretize"`` strategy is 100% reproducible, but tends
to create parcels of fairly even and geometrical shape.
The recently added ``"cluster_qr"`` option is a deterministic alternative that
tends to create the visually best partitioning on the example application
below.

================================  ================================  ================================
 ``assign_labels="kmeans"``        ``assign_labels="discretize"``    ``assign_labels="cluster_qr"``
================================  ================================  ================================
|coin_kmeans|                          |coin_discretize|                  |coin_cluster_qr|
================================  ================================  ================================

.. topic:: References:

 * `"Multiclass spectral clustering"
   <https://people.eecs.berkeley.edu/~jordan/courses/281B-spring04/readings/yu-shi.pdf>`_
   Stella X. Yu, Jianbo Shi, 2003

 * :doi:`"Simple, direct, and efficient multi-way spectral clustering"<10.1093/imaiai/iay008>`
   Anil Damle, Victor Minden, Lexing Ying, 2019

.. _spectral_clustering_graph:

Spectral Clustering Graphs
--------------------------

Spectral Clustering can also be used to partition graphs via their spectral
embeddings.  In this case, the affinity matrix is the adjacency matrix of the
graph, and SpectralClustering is initialized with `affinity='precomputed'`::

    >>> from sklearn.cluster import SpectralClustering
    >>> sc = SpectralClustering(3, affinity='precomputed', n_init=100,
    ...                         assign_labels='discretize')
    >>> sc.fit_predict(adjacency_matrix)  # doctest: +SKIP

.. topic:: References:

 * :doi:`"A Tutorial on Spectral Clustering"
   <10.1007/s11222-007-9033-z>`
   Ulrike von Luxburg, 2007

 * :doi:`"Normalized cuts and image segmentation"
   <10.1109/34.868688>`
   Jianbo Shi, Jitendra Malik, 2000

 * `"A Random Walks View of Spectral Segmentation"
   <https://citeseerx.ist.psu.edu/doc_view/pid/84a86a69315e994cfd1e0c7debb86d62d7bd1f44>`_
   Marina Meila, Jianbo Shi, 2001

 * `"On Spectral Clustering: Analysis and an algorithm"
   <https://citeseerx.ist.psu.edu/doc_view/pid/796c5d6336fc52aa84db575fb821c78918b65f58>`_
   Andrew Y. Ng, Michael I. Jordan, Yair Weiss, 2001

 * :arxiv:`"Preconditioned Spectral Clustering for Stochastic
   Block Partition Streaming Graph Challenge"
   <1708.07481>`
   David Zhuzhunashvili, Andrew Knyazev

.. _hierarchical_clustering:

Hierarchical clustering
=======================

Hierarchical clustering is a general family of clustering algorithms that
build nested clusters by merging or splitting them successively. This
hierarchy of clusters is represented as a tree (or dendrogram). The root of the
tree is the unique cluster that gathers all the samples, the leaves being the
clusters with only one sample. See the `Wikipedia page
<https://en.wikipedia.org/wiki/Hierarchical_clustering>`_ for more details.

The :class:`AgglomerativeClustering` object performs a hierarchical clustering
using a bottom up approach: each observation starts in its own cluster, and
clusters are successively merged together. The linkage criteria determines the
metric used for the merge strategy:

- **Ward** minimizes the sum of squared differences within all clusters. It is a
  variance-minimizing approach and in this sense is similar to the k-means
  objective function but tackled with an agglomerative hierarchical
  approach.
- **Maximum** or **complete linkage** minimizes the maximum distance between
  observations of pairs of clusters.
- **Average linkage** minimizes the average of the distances between all
  observations of pairs of clusters.
- **Single linkage** minimizes the distance between the closest
  observations of pairs of clusters.

:class:`AgglomerativeClustering` can also scale to large number of samples
when it is used jointly with a connectivity matrix, but is computationally
expensive when no connectivity constraints are added between samples: it
considers at each step all the possible merges.

.. topic:: :class:`FeatureAgglomeration`

   The :class:`FeatureAgglomeration` uses agglomerative clustering to
   group together features that look very similar, thus decreasing the
   number of features. It is a dimensionality reduction tool, see
   :ref:`data_reduction`.

Different linkage type: Ward, complete, average, and single linkage
-------------------------------------------------------------------

:class:`AgglomerativeClustering` supports Ward, single, average, and complete
linkage strategies.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_linkage_comparison_001.png
    :target: ../auto_examples/cluster/plot_linkage_comparison.html
    :scale: 43

Agglomerative cluster has a "rich get richer" behavior that leads to
uneven cluster sizes. In this regard, single linkage is the worst
strategy, and Ward gives the most regular sizes. However, the affinity
(or distance used in clustering) cannot be varied with Ward, thus for non
Euclidean metrics, average linkage is a good alternative. Single linkage,
while not robust to noisy data, can be computed very efficiently and can
therefore be useful to provide hierarchical clustering of larger datasets.
Single linkage can also perform well on non-globular data.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_digits_linkage.py`: exploration of the
   different linkage strategies in a real dataset.

Visualization of cluster hierarchy
----------------------------------

It's possible to visualize the tree representing the hierarchical merging of clusters
as a dendrogram. Visual inspection can often be useful for understanding the structure
of the data, though more so in the case of small sample sizes.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_dendrogram_001.png
    :target: ../auto_examples/cluster/plot_agglomerative_dendrogram.html
    :scale: 42



Adding connectivity constraints
-------------------------------

An interesting aspect of :class:`AgglomerativeClustering` is that
connectivity constraints can be added to this algorithm (only adjacent
clusters can be merged together), through a connectivity matrix that defines
for each sample the neighboring samples following a given structure of the
data. For instance, in the swiss-roll example below, the connectivity
constraints forbid the merging of points that are not adjacent on the swiss
roll, and thus avoid forming clusters that extend across overlapping folds of
the roll.

.. |unstructured| image:: ../auto_examples/cluster/images/sphx_glr_plot_ward_structured_vs_unstructured_001.png
        :target: ../auto_examples/cluster/plot_ward_structured_vs_unstructured.html
        :scale: 49

.. |structured| image:: ../auto_examples/cluster/images/sphx_glr_plot_ward_structured_vs_unstructured_002.png
        :target: ../auto_examples/cluster/plot_ward_structured_vs_unstructured.html
        :scale: 49

.. centered:: |unstructured| |structured|

These constraint are useful to impose a certain local structure, but they
also make the algorithm faster, especially when the number of the samples
is high.

The connectivity constraints are imposed via an connectivity matrix: a
scipy sparse matrix that has elements only at the intersection of a row
and a column with indices of the dataset that should be connected. This
matrix can be constructed from a-priori information: for instance, you
may wish to cluster web pages by only merging pages with a link pointing
from one to another. It can also be learned from the data, for instance
using :func:`sklearn.neighbors.kneighbors_graph` to restrict
merging to nearest neighbors as in :ref:`this example
<sphx_glr_auto_examples_cluster_plot_agglomerative_clustering.py>`, or
using :func:`sklearn.feature_extraction.image.grid_to_graph` to
enable only merging of neighboring pixels on an image, as in the
:ref:`coin <sphx_glr_auto_examples_cluster_plot_coin_ward_segmentation.py>` example.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_coin_ward_segmentation.py`: Ward clustering
   to split the image of coins in regions.

 * :ref:`sphx_glr_auto_examples_cluster_plot_ward_structured_vs_unstructured.py`: Example of
   Ward algorithm on a swiss-roll, comparison of structured approaches
   versus unstructured approaches.

 * :ref:`sphx_glr_auto_examples_cluster_plot_feature_agglomeration_vs_univariate_selection.py`:
   Example of dimensionality reduction with feature agglomeration based on
   Ward hierarchical clustering.

 * :ref:`sphx_glr_auto_examples_cluster_plot_agglomerative_clustering.py`

.. warning:: **Connectivity constraints with single, average and complete linkage**

    Connectivity constraints and single, complete or average linkage can enhance
    the 'rich getting richer' aspect of agglomerative clustering,
    particularly so if they are built with
    :func:`sklearn.neighbors.kneighbors_graph`. In the limit of a small
    number of clusters, they tend to give a few macroscopically occupied
    clusters and almost empty ones. (see the discussion in
    :ref:`sphx_glr_auto_examples_cluster_plot_agglomerative_clustering.py`).
    Single linkage is the most brittle linkage option with regard to this issue.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_001.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering.html
    :scale: 38

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_002.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering.html
    :scale: 38

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_003.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering.html
    :scale: 38

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_004.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering.html
    :scale: 38


Varying the metric
-------------------

Single, average and complete linkage can be used with a variety of distances (or
affinities), in particular Euclidean distance (*l2*), Manhattan distance
(or Cityblock, or *l1*), cosine distance, or any precomputed affinity
matrix.

* *l1* distance is often good for sparse features, or sparse noise: i.e.
  many of the features are zero, as in text mining using occurrences of
  rare words.

* *cosine* distance is interesting because it is invariant to global
  scalings of the signal.

The guidelines for choosing a metric is to use one that maximizes the
distance between samples in different classes, and minimizes that within
each class.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_metrics_005.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering_metrics.html
    :scale: 32

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_metrics_006.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering_metrics.html
    :scale: 32

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_agglomerative_clustering_metrics_007.png
    :target: ../auto_examples/cluster/plot_agglomerative_clustering_metrics.html
    :scale: 32

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_agglomerative_clustering_metrics.py`

Bisecting K-Means
-----------------

.. _bisect_k_means:

The :class:`BisectingKMeans` is an iterative variant of :class:`KMeans`, using
divisive hierarchical clustering. Instead of creating all centroids at once, centroids
are picked progressively based on a previous clustering: a cluster is split into two
new clusters repeatedly until the target number of clusters is reached.

:class:`BisectingKMeans` is more efficient than :class:`KMeans` when the number of
clusters is large since it only works on a subset of the data at each bisection
while :class:`KMeans` always works on the entire dataset.

Although :class:`BisectingKMeans` can't benefit from the advantages of the `"k-means++"`
initialization by design, it will still produce comparable results than
`KMeans(init="k-means++")` in terms of inertia at cheaper computational costs, and will
likely produce better results than `KMeans` with a random initialization.

This variant is more efficient to agglomerative clustering if the number of clusters is
small compared to the number of data points.

This variant also does not produce empty clusters.

There exist two strategies for selecting the cluster to split:
 - ``bisecting_strategy="largest_cluster"`` selects the cluster having the most points
 - ``bisecting_strategy="biggest_inertia"`` selects the cluster with biggest inertia
   (cluster with biggest Sum of Squared Errors within)

Picking by largest amount of data points in most cases produces result as
accurate as picking by inertia and is faster (especially for larger amount of data
points, where calculating error may be costly).

Picking by largest amount of data points will also likely produce clusters of similar
sizes while `KMeans` is known to produce clusters of different sizes.

Difference between Bisecting K-Means and regular K-Means can be seen on example
:ref:`sphx_glr_auto_examples_cluster_plot_bisect_kmeans.py`.
While the regular K-Means algorithm tends to create non-related clusters,
clusters from Bisecting K-Means are well ordered and create quite a visible hierarchy.

.. topic:: References:

 * `"A Comparison of Document Clustering Techniques"
   <http://www.philippe-fournier-viger.com/spmf/bisectingkmeans.pdf>`_
   Michael Steinbach, George Karypis and Vipin Kumar,
   Department of Computer Science and Egineering, University of Minnesota
   (June 2000)
 * `"Performance Analysis of K-Means and Bisecting K-Means Algorithms in Weblog Data"
   <https://ijeter.everscience.org/Manuscripts/Volume-4/Issue-8/Vol-4-issue-8-M-23.pdf>`_
   K.Abirami and Dr.P.Mayilvahanan,
   International Journal of Emerging Technologies in Engineering Research (IJETER)
   Volume 4, Issue 8, (August 2016)
 * `"Bisecting K-means Algorithm Based on K-valued Self-determining
   and Clustering Center Optimization"
   <http://www.jcomputers.us/vol13/jcp1306-01.pdf>`_
   Jian Di, Xinyue Gou
   School of Control and Computer Engineering,North China Electric Power University,
   Baoding, Hebei, China (August 2017)

.. _dbscan:

DBSCAN
======

The :class:`DBSCAN` algorithm views clusters as areas of high density
separated by areas of low density. Due to this rather generic view, clusters
found by DBSCAN can be any shape, as opposed to k-means which assumes that
clusters are convex shaped. The central component to the DBSCAN is the concept
of *core samples*, which are samples that are in areas of high density. A
cluster is therefore a set of core samples, each close to each other
(measured by some distance measure)
and a set of non-core samples that are close to a core sample (but are not
themselves core samples). There are two parameters to the algorithm,
``min_samples`` and ``eps``,
which define formally what we mean when we say *dense*.
Higher ``min_samples`` or lower ``eps``
indicate higher density necessary to form a cluster.

More formally, we define a core sample as being a sample in the dataset such
that there exist ``min_samples`` other samples within a distance of
``eps``, which are defined as *neighbors* of the core sample. This tells
us that the core sample is in a dense area of the vector space. A cluster
is a set of core samples that can be built by recursively taking a core
sample, finding all of its neighbors that are core samples, finding all of
*their* neighbors that are core samples, and so on. A cluster also has a
set of non-core samples, which are samples that are neighbors of a core sample
in the cluster but are not themselves core samples. Intuitively, these samples
are on the fringes of a cluster.

Any core sample is part of a cluster, by definition. Any sample that is not a
core sample, and is at least ``eps`` in distance from any core sample, is
considered an outlier by the algorithm.

While the parameter ``min_samples`` primarily controls how tolerant the
algorithm is towards noise (on noisy and large data sets it may be desirable
to increase this parameter), the parameter ``eps`` is *crucial to choose
appropriately* for the data set and distance function and usually cannot be
left at the default value. It controls the local neighborhood of the points.
When chosen too small, most data will not be clustered at all (and labeled
as ``-1`` for "noise"). When chosen too large, it causes close clusters to
be merged into one cluster, and eventually the entire data set to be returned
as a single cluster. Some heuristics for choosing this parameter have been
discussed in the literature, for example based on a knee in the nearest neighbor
distances plot (as discussed in the references below).

In the figure below, the color indicates cluster membership, with large circles
indicating core samples found by the algorithm. Smaller circles are non-core
samples that are still part of a cluster. Moreover, the outliers are indicated
by black points below.

.. |dbscan_results| image:: ../auto_examples/cluster/images/sphx_glr_plot_dbscan_002.png
        :target: ../auto_examples/cluster/plot_dbscan.html
        :scale: 50

.. centered:: |dbscan_results|

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_cluster_plot_dbscan.py`

.. topic:: Implementation

    The DBSCAN algorithm is deterministic, always generating the same clusters
    when given the same data in the same order.  However, the results can differ when
    data is provided in a different order. First, even though the core samples
    will always be assigned to the same clusters, the labels of those clusters
    will depend on the order in which those samples are encountered in the data.
    Second and more importantly, the clusters to which non-core samples are assigned
    can differ depending on the data order.  This would happen when a non-core sample
    has a distance lower than ``eps`` to two core samples in different clusters. By the
    triangular inequality, those two core samples must be more distant than
    ``eps`` from each other, or they would be in the same cluster. The non-core
    sample is assigned to whichever cluster is generated first in a pass
    through the data, and so the results will depend on the data ordering.

    The current implementation uses ball trees and kd-trees
    to determine the neighborhood of points,
    which avoids calculating the full distance matrix
    (as was done in scikit-learn versions before 0.14).
    The possibility to use custom metrics is retained;
    for details, see :class:`~sklearn.neighbors.NearestNeighbors`.

.. topic:: Memory consumption for large sample sizes

    This implementation is by default not memory efficient because it constructs
    a full pairwise similarity matrix in the case where kd-trees or ball-trees cannot
    be used (e.g., with sparse matrices). This matrix will consume :math:`n^2` floats.
    A couple of mechanisms for getting around this are:

    - Use :ref:`OPTICS <optics>` clustering in conjunction with the
      `extract_dbscan` method. OPTICS clustering also calculates the full
      pairwise matrix, but only keeps one row in memory at a time (memory
      complexity n).

    - A sparse radius neighborhood graph (where missing entries are presumed to
      be out of eps) can be precomputed in a memory-efficient way and dbscan
      can be run over this with ``metric='precomputed'``.  See
      :meth:`sklearn.neighbors.NearestNeighbors.radius_neighbors_graph`.

    - The dataset can be compressed, either by removing exact duplicates if
      these occur in your data, or by using BIRCH. Then you only have a
      relatively small number of representatives for a large number of points.
      You can then provide a ``sample_weight`` when fitting DBSCAN.

.. topic:: References:

 * `"A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases
   with Noise" <https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf>`_
   Ester, M., H. P. Kriegel, J. Sander, and X. Xu,
   In Proceedings of the 2nd International Conference on Knowledge Discovery
   and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996

 * :doi:`"DBSCAN revisited, revisited: why and how you should (still) use DBSCAN."
   <10.1145/3068335>`
   Schubert, E., Sander, J., Ester, M., Kriegel, H. P., & Xu, X. (2017).
   In ACM Transactions on Database Systems (TODS), 42(3), 19.

.. _hdbscan:

HDBSCAN
=======

The :class:`HDBSCAN` algorithm can be seen as an extension of :class:`DBSCAN`
and :class:`OPTICS`. Specifically, :class:`DBSCAN` assumes that the clustering
criterion (i.e. density requirement) is *globally homogeneous*.
In other words, :class:`DBSCAN` may struggle to successfully capture clusters
with different densities.
:class:`HDBSCAN` alleviates this assumption and explores all possible density
scales by building an alternative representation of the clustering problem.

.. note::

  This implementation is adapted from the original implementation of HDBSCAN,
  `scikit-learn-contrib/hdbscan <https://github.com/scikit-learn-contrib/hdbscan>`_ based on [LJ2017]_.

Mutual Reachability Graph
-------------------------

HDBSCAN first defines :math:`d_c(x_p)`, the *core distance* of a sample :math:`x_p`, as the
distance to its `min_samples` th-nearest neighbor, counting itself. For example,
if `min_samples=5` and :math:`x_*` is the 5th-nearest neighbor of :math:`x_p`
then the core distance is:

.. math:: d_c(x_p)=d(x_p, x_*).

Next it defines :math:`d_m(x_p, x_q)`, the *mutual reachability distance* of two points
:math:`x_p, x_q`, as:

.. math:: d_m(x_p, x_q) = \max\{d_c(x_p), d_c(x_q), d(x_p, x_q)\}

These two notions allow us to construct the *mutual reachability graph*
:math:`G_{ms}` defined for a fixed choice of `min_samples` by associating each
sample :math:`x_p` with a vertex of the graph, and thus edges between points
:math:`x_p, x_q` are the mutual reachability distance :math:`d_m(x_p, x_q)`
between them. We may build subsets of this graph, denoted as
:math:`G_{ms,\varepsilon}`, by removing any edges with value greater than :math:`\varepsilon`:
from the original graph. Any points whose core distance is less than :math:`\varepsilon`:
are at this staged marked as noise. The remaining points are then clustered by
finding the connected components of this trimmed graph.

.. note::

  Taking the connected components of a trimmed graph :math:`G_{ms,\varepsilon}` is
  equivalent to running DBSCAN* with `min_samples` and :math:`\varepsilon`. DBSCAN* is a
  slightly modified version of DBSCAN mentioned in [CM2013]_.

Hierarchical Clustering
-----------------------
HDBSCAN can be seen as an algorithm which performs DBSCAN* clustering across all
values of :math:`\varepsilon`. As mentioned prior, this is equivalent to finding the connected
components of the mutual reachability graphs for all values of :math:`\varepsilon`. To do this
efficiently, HDBSCAN first extracts a minimum spanning tree (MST) from the fully
-connected mutual reachability graph, then greedily cuts the edges with highest
weight. An outline of the HDBSCAN algorithm is as follows:

  1. Extract the MST of :math:`G_{ms}`
  2. Extend the MST by adding a "self edge" for each vertex, with weight equal
     to the core distance of the underlying sample.
  3. Initialize a single cluster and label for the MST.
  4. Remove the edge with the greatest weight from the MST (ties are
     removed simultaneously).
  5. Assign cluster labels to the connected components which contain the
     end points of the now-removed edge. If the component does not have at least
     one edge it is instead assigned a "null" label marking it as noise.
  6. Repeat 4-5 until there are no more connected components.

HDBSCAN is therefore able to obtain all possible partitions achievable by
DBSCAN* for a fixed choice of `min_samples` in a hierarchical fashion.
Indeed, this allows HDBSCAN to perform clustering across multiple densities
and as such it no longer needs :math:`\varepsilon` to be given as a hyperparameter. Instead
it relies solely on the choice of `min_samples`, which tends to be a more robust
hyperparameter.

.. |hdbscan_ground_truth| image:: ../auto_examples/cluster/images/sphx_glr_plot_hdbscan_005.png
        :target: ../auto_examples/cluster/plot_hdbscan.html
        :scale: 75
.. |hdbscan_results| image:: ../auto_examples/cluster/images/sphx_glr_plot_hdbscan_007.png
        :target: ../auto_examples/cluster/plot_hdbscan.html
        :scale: 75

.. centered:: |hdbscan_ground_truth|
.. centered:: |hdbscan_results|

HDBSCAN can be smoothed with an additional hyperparameter `min_cluster_size`
which specifies that during the hierarchical clustering, components with fewer
than `minimum_cluster_size` many samples are considered noise. In practice, one
can set `minimum_cluster_size = min_samples` to couple the parameters and
simplify the hyperparameter space.

.. topic:: References:

 .. [CM2013] Campello, R.J.G.B., Moulavi, D., Sander, J. (2013). Density-Based Clustering
   Based on Hierarchical Density Estimates. In: Pei, J., Tseng, V.S., Cao, L.,
   Motoda, H., Xu, G. (eds) Advances in Knowledge Discovery and Data Mining.
   PAKDD 2013. Lecture Notes in Computer Science(), vol 7819. Springer, Berlin,
   Heidelberg.
   :doi:`Density-Based Clustering Based on Hierarchical Density Estimates <10.1007/978-3-642-37456-2_14>`

 .. [LJ2017] L. McInnes and J. Healy, (2017). Accelerated Hierarchical Density Based
   Clustering. In: IEEE International Conference on Data Mining Workshops (ICDMW),
   2017, pp. 33-42.
   :doi:`Accelerated Hierarchical Density Based Clustering <10.1109/ICDMW.2017.12>`

.. _optics:

OPTICS
======

The :class:`OPTICS` algorithm shares many similarities with the :class:`DBSCAN`
algorithm, and can be considered a generalization of DBSCAN that relaxes the
``eps`` requirement from a single value to a value range. The key difference
between DBSCAN and OPTICS is that the OPTICS algorithm builds a *reachability*
graph, which assigns each sample both a ``reachability_`` distance, and a spot
within the cluster ``ordering_`` attribute; these two attributes are assigned
when the model is fitted, and are used to determine cluster membership. If
OPTICS is run with the default value of *inf* set for ``max_eps``, then DBSCAN
style cluster extraction can be performed repeatedly in linear time for any
given ``eps`` value using the ``cluster_optics_dbscan`` method. Setting
``max_eps`` to a lower value will result in shorter run times, and can be
thought of as the maximum neighborhood radius from each point to find other
potential reachable points.

.. |optics_results| image:: ../auto_examples/cluster/images/sphx_glr_plot_optics_001.png
        :target: ../auto_examples/cluster/plot_optics.html
        :scale: 50

.. centered:: |optics_results|

The *reachability* distances generated by OPTICS allow for variable density
extraction of clusters within a single data set. As shown in the above plot,
combining *reachability* distances and data set ``ordering_`` produces a
*reachability plot*, where point density is represented on the Y-axis, and
points are ordered such that nearby points are adjacent. 'Cutting' the
reachability plot at a single value produces DBSCAN like results; all points
above the 'cut' are classified as noise, and each time that there is a break
when reading from left to right signifies a new cluster. The default cluster
extraction with OPTICS looks at the steep slopes within the graph to find
clusters, and the user can define what counts as a steep slope using the
parameter ``xi``. There are also other possibilities for analysis on the graph
itself, such as generating hierarchical representations of the data through
reachability-plot dendrograms, and the hierarchy of clusters detected by the
algorithm can be accessed through the ``cluster_hierarchy_`` parameter. The
plot above has been color-coded so that cluster colors in planar space match
the linear segment clusters of the reachability plot. Note that the blue and
red clusters are adjacent in the reachability plot, and can be hierarchically
represented as children of a larger parent cluster.

.. topic:: Examples:

     * :ref:`sphx_glr_auto_examples_cluster_plot_optics.py`


.. topic:: Comparison with DBSCAN

    The results from OPTICS ``cluster_optics_dbscan`` method and DBSCAN are
    very similar, but not always identical; specifically, labeling of periphery
    and noise points. This is in part because the first samples of each dense
    area processed by OPTICS have a large reachability value while being close
    to other points in their area, and will thus sometimes be marked as noise
    rather than periphery. This affects adjacent points when they are
    considered as candidates for being marked as either periphery or noise.

    Note that for any single value of ``eps``, DBSCAN will tend to have a
    shorter run time than OPTICS; however, for repeated runs at varying ``eps``
    values, a single run of OPTICS may require less cumulative runtime than
    DBSCAN. It is also important to note that OPTICS' output is close to
    DBSCAN's only if ``eps`` and ``max_eps`` are close.

.. topic:: Computational Complexity

    Spatial indexing trees are used to avoid calculating the full distance
    matrix, and allow for efficient memory usage on large sets of samples.
    Different distance metrics can be supplied via the ``metric`` keyword.

    For large datasets, similar (but not identical) results can be obtained via
    :class:`HDBSCAN`. The HDBSCAN implementation is
    multithreaded, and has better algorithmic runtime complexity than OPTICS,
    at the cost of worse memory scaling. For extremely large datasets that
    exhaust system memory using HDBSCAN, OPTICS will maintain :math:`n` (as opposed
    to :math:`n^2`) memory scaling; however, tuning of the ``max_eps`` parameter
    will likely need to be used to give a solution in a reasonable amount of
    wall time.

.. topic:: References:

 *  "OPTICS: ordering points to identify the clustering structure."
    Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander.
    In ACM Sigmod Record, vol. 28, no. 2, pp. 49-60. ACM, 1999.

.. _birch:

BIRCH
=====

The :class:`Birch` builds a tree called the Clustering Feature Tree (CFT)
for the given data. The data is essentially lossy compressed to a set of
Clustering Feature nodes (CF Nodes). The CF Nodes have a number of
subclusters called Clustering Feature subclusters (CF Subclusters)
and these CF Subclusters located in the non-terminal CF Nodes
can have CF Nodes as children.

The CF Subclusters hold the necessary information for clustering which prevents
the need to hold the entire input data in memory. This information includes:

- Number of samples in a subcluster.
- Linear Sum - An n-dimensional vector holding the sum of all samples
- Squared Sum - Sum of the squared L2 norm of all samples.
- Centroids - To avoid recalculation linear sum / n_samples.
- Squared norm of the centroids.

The BIRCH algorithm has two parameters, the threshold and the branching factor.
The branching factor limits the number of subclusters in a node and the
threshold limits the distance between the entering sample and the existing
subclusters.

This algorithm can be viewed as an instance or data reduction method,
since it reduces the input data to a set of subclusters which are obtained directly
from the leaves of the CFT. This reduced data can be further processed by feeding
it into a global clusterer. This global clusterer can be set by ``n_clusters``.
If ``n_clusters`` is set to None, the subclusters from the leaves are directly
read off, otherwise a global clustering step labels these subclusters into global
clusters (labels) and the samples are mapped to the global label of the nearest subcluster.

**Algorithm description:**

- A new sample is inserted into the root of the CF Tree which is a CF Node.
  It is then merged with the subcluster of the root, that has the smallest
  radius after merging, constrained by the threshold and branching factor conditions.
  If the subcluster has any child node, then this is done repeatedly till it reaches
  a leaf. After finding the nearest subcluster in the leaf, the properties of this
  subcluster and the parent subclusters are recursively updated.

- If the radius of the subcluster obtained by merging the new sample and the
  nearest subcluster is greater than the square of the threshold and if the
  number of subclusters is greater than the branching factor, then a space is temporarily
  allocated to this new sample. The two farthest subclusters are taken and
  the subclusters are divided into two groups on the basis of the distance
  between these subclusters.

- If this split node has a parent subcluster and there is room
  for a new subcluster, then the parent is split into two. If there is no room,
  then this node is again split into two and the process is continued
  recursively, till it reaches the root.

**BIRCH or MiniBatchKMeans?**

 - BIRCH does not scale very well to high dimensional data. As a rule of thumb if
   ``n_features`` is greater than twenty, it is generally better to use MiniBatchKMeans.
 - If the number of instances of data needs to be reduced, or if one wants a
   large number of subclusters either as a preprocessing step or otherwise,
   BIRCH is more useful than MiniBatchKMeans.


**How to use partial_fit?**

To avoid the computation of global clustering, for every call of ``partial_fit``
the user is advised

 1. To set ``n_clusters=None`` initially
 2. Train all data by multiple calls to partial_fit.
 3. Set ``n_clusters`` to a required value using
    ``brc.set_params(n_clusters=n_clusters)``.
 4. Call ``partial_fit`` finally with no arguments, i.e. ``brc.partial_fit()``
    which performs the global clustering.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_birch_vs_minibatchkmeans_001.png
    :target: ../auto_examples/cluster/plot_birch_vs_minibatchkmeans.html

.. topic:: References:

 * Tian Zhang, Raghu Ramakrishnan, Maron Livny
   BIRCH: An efficient data clustering method for large databases.
   https://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf

 * Roberto Perdisci
   JBirch - Java implementation of BIRCH clustering algorithm
   https://code.google.com/archive/p/jbirch


.. _clustering_evaluation:

Clustering performance evaluation
=================================

Evaluating the performance of a clustering algorithm is not as trivial as
counting the number of errors or the precision and recall of a supervised
classification algorithm. In particular any evaluation metric should not
take the absolute values of the cluster labels into account but rather
if this clustering define separations of the data similar to some ground
truth set of classes or satisfying some assumption such that members
belong to the same class are more similar than members of different
classes according to some similarity metric.

.. currentmodule:: sklearn.metrics

.. _rand_score:
.. _adjusted_rand_score:

Rand index
----------

Given the knowledge of the ground truth class assignments
``labels_true`` and our clustering algorithm assignments of the same
samples ``labels_pred``, the **(adjusted or unadjusted) Rand index**
is a function that measures the **similarity** of the two assignments,
ignoring permutations::

  >>> from sklearn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]
  >>> metrics.rand_score(labels_true, labels_pred)
  0.66...

The Rand index does not ensure to obtain a value close to 0.0 for a
random labelling. The adjusted Rand index **corrects for chance** and
will give such a baseline.

  >>> metrics.adjusted_rand_score(labels_true, labels_pred)
  0.24...

As with all clustering metrics, one can permute 0 and 1 in the predicted
labels, rename 2 to 3, and get the same score::

  >>> labels_pred = [1, 1, 0, 0, 3, 3]
  >>> metrics.rand_score(labels_true, labels_pred)
  0.66...
  >>> metrics.adjusted_rand_score(labels_true, labels_pred)
  0.24...

Furthermore, both :func:`rand_score` :func:`adjusted_rand_score` are
**symmetric**: swapping the argument does not change the scores. They can
thus be used as **consensus measures**::

  >>> metrics.rand_score(labels_pred, labels_true)
  0.66...
  >>> metrics.adjusted_rand_score(labels_pred, labels_true)
  0.24...

Perfect labeling is scored 1.0::

  >>> labels_pred = labels_true[:]
  >>> metrics.rand_score(labels_true, labels_pred)
  1.0
  >>> metrics.adjusted_rand_score(labels_true, labels_pred)
  1.0

Poorly agreeing labels (e.g. independent labelings) have lower scores,
and for the adjusted Rand index the score will be negative or close to
zero. However, for the unadjusted Rand index the score, while lower,
will not necessarily be close to zero.::

  >>> labels_true = [0, 0, 0, 0, 0, 0, 1, 1]
  >>> labels_pred = [0, 1, 2, 3, 4, 5, 5, 6]
  >>> metrics.rand_score(labels_true, labels_pred)
  0.39...
  >>> metrics.adjusted_rand_score(labels_true, labels_pred)
  -0.07...


Advantages
~~~~~~~~~~

- **Interpretability**: The unadjusted Rand index is proportional
  to the number of sample pairs whose labels are the same in both
  `labels_pred` and `labels_true`, or are different in both.

- **Random (uniform) label assignments have an adjusted Rand index
  score close to 0.0** for any value of ``n_clusters`` and
  ``n_samples`` (which is not the case for the unadjusted Rand index
  or the V-measure for instance).

- **Bounded range**: Lower values indicate different labelings,
  similar clusterings have a high (adjusted or unadjusted) Rand index,
  1.0 is the perfect match score. The score range is [0, 1] for the
  unadjusted Rand index and [-1, 1] for the adjusted Rand index.

- **No assumption is made on the cluster structure**: The (adjusted or
  unadjusted) Rand index can be used to compare all kinds of
  clustering algorithms, and can be used to compare clustering
  algorithms such as k-means which assumes isotropic blob shapes with
  results of spectral clustering algorithms which can find cluster
  with "folded" shapes.


Drawbacks
~~~~~~~~~

- Contrary to inertia, the **(adjusted or unadjusted) Rand index
  requires knowledge of the ground truth classes** which is almost
  never available in practice or requires manual assignment by human
  annotators (as in the supervised learning setting).

  However (adjusted or unadjusted) Rand index can also be useful in a
  purely unsupervised setting as a building block for a Consensus
  Index that can be used for clustering model selection (TODO).

- The **unadjusted Rand index is often close to 1.0** even if the
  clusterings themselves differ significantly. This can be understood
  when interpreting the Rand index as the accuracy of element pair
  labeling resulting from the clusterings: In practice there often is
  a majority of element pairs that are assigned the ``different`` pair
  label under both the predicted and the ground truth clustering
  resulting in a high proportion of pair labels that agree, which
  leads subsequently to a high score.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_adjusted_for_chance_measures.py`:
   Analysis of the impact of the dataset size on the value of
   clustering measures for random assignments.


Mathematical formulation
~~~~~~~~~~~~~~~~~~~~~~~~

If C is a ground truth class assignment and K the clustering, let us
define :math:`a` and :math:`b` as:

- :math:`a`, the number of pairs of elements that are in the same set
  in C and in the same set in K

- :math:`b`, the number of pairs of elements that are in different sets
  in C and in different sets in K

The unadjusted Rand index is then given by:

.. math:: \text{RI} = \frac{a + b}{C_2^{n_{samples}}}

where :math:`C_2^{n_{samples}}` is the total number of possible pairs
in the dataset. It does not matter if the calculation is performed on
ordered pairs or unordered pairs as long as the calculation is
performed consistently.

However, the Rand index does not guarantee that random label assignments
will get a value close to zero (esp. if the number of clusters is in
the same order of magnitude as the number of samples).

To counter this effect we can discount the expected RI :math:`E[\text{RI}]` of
random labelings by defining the adjusted Rand index as follows:

.. math:: \text{ARI} = \frac{\text{RI} - E[\text{RI}]}{\max(\text{RI}) - E[\text{RI}]}

.. topic:: References

 * `Comparing Partitions
   <https://link.springer.com/article/10.1007%2FBF01908075>`_
   L. Hubert and P. Arabie, Journal of Classification 1985

 * `Properties of the Hubert-Arabie adjusted Rand index
   <https://psycnet.apa.org/record/2004-17801-007>`_
   D. Steinley, Psychological Methods 2004

 * `Wikipedia entry for the Rand index
   <https://en.wikipedia.org/wiki/Rand_index>`_

 * `Wikipedia entry for the adjusted Rand index
   <https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index>`_


.. _mutual_info_score:

Mutual Information based scores
-------------------------------

Given the knowledge of the ground truth class assignments ``labels_true`` and
our clustering algorithm assignments of the same samples ``labels_pred``, the
**Mutual Information** is a function that measures the **agreement** of the two
assignments, ignoring permutations.  Two different normalized versions of this
measure are available, **Normalized Mutual Information (NMI)** and **Adjusted
Mutual Information (AMI)**. NMI is often used in the literature, while AMI was
proposed more recently and is **normalized against chance**::

  >>> from sklearn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]

  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +SKIP
  0.22504...

One can permute 0 and 1 in the predicted labels, rename 2 to 3 and get
the same score::

  >>> labels_pred = [1, 1, 0, 0, 3, 3]
  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +SKIP
  0.22504...

All, :func:`mutual_info_score`, :func:`adjusted_mutual_info_score` and
:func:`normalized_mutual_info_score` are symmetric: swapping the argument does
not change the score. Thus they can be used as a **consensus measure**::

  >>> metrics.adjusted_mutual_info_score(labels_pred, labels_true)  # doctest: +SKIP
  0.22504...

Perfect labeling is scored 1.0::

  >>> labels_pred = labels_true[:]
  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +SKIP
  1.0

  >>> metrics.normalized_mutual_info_score(labels_true, labels_pred)  # doctest: +SKIP
  1.0

This is not true for ``mutual_info_score``, which is therefore harder to judge::

  >>> metrics.mutual_info_score(labels_true, labels_pred)  # doctest: +SKIP
  0.69...

Bad (e.g. independent labelings) have non-positive scores::

  >>> labels_true = [0, 1, 2, 0, 3, 4, 5, 1]
  >>> labels_pred = [1, 1, 0, 0, 2, 2, 2, 2]
  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +SKIP
  -0.10526...


Advantages
~~~~~~~~~~

- **Random (uniform) label assignments have a AMI score close to 0.0**
  for any value of ``n_clusters`` and ``n_samples`` (which is not the
  case for raw Mutual Information or the V-measure for instance).

- **Upper bound  of 1**:  Values close to zero indicate two label
  assignments that are largely independent, while values close to one
  indicate significant agreement. Further, an AMI of exactly 1 indicates
  that the two label assignments are equal (with or without permutation).


Drawbacks
~~~~~~~~~

- Contrary to inertia, **MI-based measures require the knowledge
  of the ground truth classes** while almost never available in practice or
  requires manual assignment by human annotators (as in the supervised learning
  setting).

  However MI-based measures can also be useful in purely unsupervised setting as a
  building block for a Consensus Index that can be used for clustering
  model selection.

- NMI and MI are not adjusted against chance.


.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_adjusted_for_chance_measures.py`: Analysis of
   the impact of the dataset size on the value of clustering measures
   for random assignments. This example also includes the Adjusted Rand
   Index.


Mathematical formulation
~~~~~~~~~~~~~~~~~~~~~~~~

Assume two label assignments (of the same N objects), :math:`U` and :math:`V`.
Their entropy is the amount of uncertainty for a partition set, defined by:

.. math:: H(U) = - \sum_{i=1}^{|U|}P(i)\log(P(i))

where :math:`P(i) = |U_i| / N` is the probability that an object picked at
random from :math:`U` falls into class :math:`U_i`. Likewise for :math:`V`:

.. math:: H(V) = - \sum_{j=1}^{|V|}P'(j)\log(P'(j))

With :math:`P'(j) = |V_j| / N`. The mutual information (MI) between :math:`U`
and :math:`V` is calculated by:

.. math:: \text{MI}(U, V) = \sum_{i=1}^{|U|}\sum_{j=1}^{|V|}P(i, j)\log\left(\frac{P(i,j)}{P(i)P'(j)}\right)

where :math:`P(i, j) = |U_i \cap V_j| / N` is the probability that an object
picked at random falls into both classes :math:`U_i` and :math:`V_j`.

It also can be expressed in set cardinality formulation:

.. math:: \text{MI}(U, V) = \sum_{i=1}^{|U|} \sum_{j=1}^{|V|} \frac{|U_i \cap V_j|}{N}\log\left(\frac{N|U_i \cap V_j|}{|U_i||V_j|}\right)

The normalized mutual information is defined as

.. math:: \text{NMI}(U, V) = \frac{\text{MI}(U, V)}{\text{mean}(H(U), H(V))}

This value of the mutual information and also the normalized variant is not
adjusted for chance and will tend to increase as the number of different labels
(clusters) increases, regardless of the actual amount of "mutual information"
between the label assignments.

The expected value for the mutual information can be calculated using the
following equation [VEB2009]_. In this equation,
:math:`a_i = |U_i|` (the number of elements in :math:`U_i`) and
:math:`b_j = |V_j|` (the number of elements in :math:`V_j`).


.. math:: E[\text{MI}(U,V)]=\sum_{i=1}^{|U|} \sum_{j=1}^{|V|} \sum_{n_{ij}=(a_i+b_j-N)^+
   }^{\min(a_i, b_j)} \frac{n_{ij}}{N}\log \left( \frac{ N.n_{ij}}{a_i b_j}\right)
   \frac{a_i!b_j!(N-a_i)!(N-b_j)!}{N!n_{ij}!(a_i-n_{ij})!(b_j-n_{ij})!
   (N-a_i-b_j+n_{ij})!}

Using the expected value, the adjusted mutual information can then be
calculated using a similar form to that of the adjusted Rand index:

.. math:: \text{AMI} = \frac{\text{MI} - E[\text{MI}]}{\text{mean}(H(U), H(V)) - E[\text{MI}]}

For normalized mutual information and adjusted mutual information, the normalizing
value is typically some *generalized* mean of the entropies of each clustering.
Various generalized means exist, and no firm rules exist for preferring one over the
others.  The decision is largely a field-by-field basis; for instance, in community
detection, the arithmetic mean is most common. Each
normalizing method provides "qualitatively similar behaviours" [YAT2016]_. In our
implementation, this is controlled by the ``average_method`` parameter.

Vinh et al. (2010) named variants of NMI and AMI by their averaging method [VEB2010]_. Their
'sqrt' and 'sum' averages are the geometric and arithmetic means; we use these
more broadly common names.

.. topic:: References

 * Strehl, Alexander, and Joydeep Ghosh (2002). "Cluster ensembles – a
   knowledge reuse framework for combining multiple partitions". Journal of
   Machine Learning Research 3: 583–617.
   `doi:10.1162/153244303321897735 <http://strehl.com/download/strehl-jmlr02.pdf>`_.

 * `Wikipedia entry for the (normalized) Mutual Information
   <https://en.wikipedia.org/wiki/Mutual_Information>`_

 * `Wikipedia entry for the Adjusted Mutual Information
   <https://en.wikipedia.org/wiki/Adjusted_Mutual_Information>`_

 .. [VEB2009] Vinh, Epps, and Bailey, (2009). "Information theoretic measures
   for clusterings comparison". Proceedings of the 26th Annual International
   Conference on Machine Learning - ICML '09.
   `doi:10.1145/1553374.1553511 <https://dl.acm.org/citation.cfm?doid=1553374.1553511>`_.
   ISBN 9781605585161.

 .. [VEB2010] Vinh, Epps, and Bailey, (2010). "Information Theoretic Measures for
   Clusterings Comparison: Variants, Properties, Normalization and
   Correction for Chance". JMLR
   <https://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf>

 .. [YAT2016] Yang, Algesheimer, and Tessone, (2016). "A comparative analysis of
   community
   detection algorithms on artificial networks". Scientific Reports 6: 30750.
   `doi:10.1038/srep30750 <https://www.nature.com/articles/srep30750>`_.



.. _homogeneity_completeness:

Homogeneity, completeness and V-measure
---------------------------------------

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

  >>> from sklearn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]

  >>> metrics.homogeneity_score(labels_true, labels_pred)
  0.66...

  >>> metrics.completeness_score(labels_true, labels_pred)
  0.42...

Their harmonic mean called **V-measure** is computed by
:func:`v_measure_score`::

  >>> metrics.v_measure_score(labels_true, labels_pred)
  0.51...

This function's formula is as follows:

.. math:: v = \frac{(1 + \beta) \times \text{homogeneity} \times \text{completeness}}{(\beta \times \text{homogeneity} + \text{completeness})}

`beta` defaults to a value of 1.0, but for using a value less than 1 for beta::

  >>> metrics.v_measure_score(labels_true, labels_pred, beta=0.6)
  0.54...

more weight will be attributed to homogeneity, and using a value greater than 1::

  >>> metrics.v_measure_score(labels_true, labels_pred, beta=1.8)
  0.48...

more weight will be attributed to completeness.

The V-measure is actually equivalent to the mutual information (NMI)
discussed above, with the aggregation function being the arithmetic mean [B2011]_.

Homogeneity, completeness and V-measure can be computed at once using
:func:`homogeneity_completeness_v_measure` as follows::

  >>> metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
  (0.66..., 0.42..., 0.51...)

The following clustering assignment is slightly better, since it is
homogeneous but not complete::

  >>> labels_pred = [0, 0, 0, 1, 2, 2]
  >>> metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
  (1.0, 0.68..., 0.81...)

.. note::

  :func:`v_measure_score` is **symmetric**: it can be used to evaluate
  the **agreement** of two independent assignments on the same dataset.

  This is not the case for :func:`completeness_score` and
  :func:`homogeneity_score`: both are bound by the relationship::

    homogeneity_score(a, b) == completeness_score(b, a)


Advantages
~~~~~~~~~~

- **Bounded scores**: 0.0 is as bad as it can be, 1.0 is a perfect score.

- Intuitive interpretation: clustering with bad V-measure can be
  **qualitatively analyzed in terms of homogeneity and completeness**
  to better feel what 'kind' of mistakes is done by the assignment.

- **No assumption is made on the cluster structure**: can be used
  to compare clustering algorithms such as k-means which assumes isotropic
  blob shapes with results of spectral clustering algorithms which can
  find cluster with "folded" shapes.


Drawbacks
~~~~~~~~~

- The previously introduced metrics are **not normalized with regards to
  random labeling**: this means that depending on the number of samples,
  clusters and ground truth classes, a completely random labeling will
  not always yield the same values for homogeneity, completeness and
  hence v-measure. In particular **random labeling won't yield zero
  scores especially when the number of clusters is large**.

  This problem can safely be ignored when the number of samples is more
  than a thousand and the number of clusters is less than 10. **For
  smaller sample sizes or larger number of clusters it is safer to use
  an adjusted index such as the Adjusted Rand Index (ARI)**.

.. figure:: ../auto_examples/cluster/images/sphx_glr_plot_adjusted_for_chance_measures_001.png
   :target: ../auto_examples/cluster/plot_adjusted_for_chance_measures.html
   :align: center
   :scale: 100

- These metrics **require the knowledge of the ground truth classes** while
  almost never available in practice or requires manual assignment by
  human annotators (as in the supervised learning setting).


.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_adjusted_for_chance_measures.py`: Analysis of
   the impact of the dataset size on the value of clustering measures
   for random assignments.


Mathematical formulation
~~~~~~~~~~~~~~~~~~~~~~~~

Homogeneity and completeness scores are formally given by:

.. math:: h = 1 - \frac{H(C|K)}{H(C)}

.. math:: c = 1 - \frac{H(K|C)}{H(K)}

where :math:`H(C|K)` is the **conditional entropy of the classes given
the cluster assignments** and is given by:

.. math:: H(C|K) = - \sum_{c=1}^{|C|} \sum_{k=1}^{|K|} \frac{n_{c,k}}{n}
          \cdot \log\left(\frac{n_{c,k}}{n_k}\right)

and :math:`H(C)` is the **entropy of the classes** and is given by:

.. math:: H(C) = - \sum_{c=1}^{|C|} \frac{n_c}{n} \cdot \log\left(\frac{n_c}{n}\right)

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

 * `V-Measure: A conditional entropy-based external cluster evaluation
   measure <https://aclweb.org/anthology/D/D07/D07-1043.pdf>`_
   Andrew Rosenberg and Julia Hirschberg, 2007

 .. [B2011] `Identication and Characterization of Events in Social Media
   <http://www.cs.columbia.edu/~hila/hila-thesis-distributed.pdf>`_, Hila
   Becker, PhD Thesis.

.. _fowlkes_mallows_scores:

Fowlkes-Mallows scores
----------------------

The Fowlkes-Mallows index (:func:`sklearn.metrics.fowlkes_mallows_score`) can be
used when the ground truth class assignments of the samples is known. The
Fowlkes-Mallows score FMI is defined as the geometric mean of the
pairwise precision and recall:

.. math:: \text{FMI} = \frac{\text{TP}}{\sqrt{(\text{TP} + \text{FP}) (\text{TP} + \text{FN})}}

Where ``TP`` is the number of **True Positive** (i.e. the number of pair
of points that belong to the same clusters in both the true labels and the
predicted labels), ``FP`` is the number of **False Positive** (i.e. the number
of pair of points that belong to the same clusters in the true labels and not
in the predicted labels) and ``FN`` is the number of **False Negative** (i.e the
number of pair of points that belongs in the same clusters in the predicted
labels and not in the true labels).

The score ranges from 0 to 1. A high value indicates a good similarity
between two clusters.

  >>> from sklearn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]

  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)
  0.47140...

One can permute 0 and 1 in the predicted labels, rename 2 to 3 and get
the same score::

  >>> labels_pred = [1, 1, 0, 0, 3, 3]

  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)
  0.47140...

Perfect labeling is scored 1.0::

  >>> labels_pred = labels_true[:]
  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)
  1.0

Bad (e.g. independent labelings) have zero scores::

  >>> labels_true = [0, 1, 2, 0, 3, 4, 5, 1]
  >>> labels_pred = [1, 1, 0, 0, 2, 2, 2, 2]
  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)
  0.0

Advantages
~~~~~~~~~~

- **Random (uniform) label assignments have a FMI score close to 0.0**
  for any value of ``n_clusters`` and ``n_samples`` (which is not the
  case for raw Mutual Information or the V-measure for instance).

- **Upper-bounded at 1**:  Values close to zero indicate two label
  assignments that are largely independent, while values close to one
  indicate significant agreement. Further, values of exactly 0 indicate
  **purely** independent label assignments and a FMI of exactly 1 indicates
  that the two label assignments are equal (with or without permutation).

- **No assumption is made on the cluster structure**: can be used
  to compare clustering algorithms such as k-means which assumes isotropic
  blob shapes with results of spectral clustering algorithms which can
  find cluster with "folded" shapes.


Drawbacks
~~~~~~~~~

- Contrary to inertia, **FMI-based measures require the knowledge
  of the ground truth classes** while almost never available in practice or
  requires manual assignment by human annotators (as in the supervised learning
  setting).

.. topic:: References

  * E. B. Fowkles and C. L. Mallows, 1983. "A method for comparing two
    hierarchical clusterings". Journal of the American Statistical Association.
    https://www.tandfonline.com/doi/abs/10.1080/01621459.1983.10478008

  * `Wikipedia entry for the Fowlkes-Mallows Index
    <https://en.wikipedia.org/wiki/Fowlkes-Mallows_index>`_

.. _silhouette_coefficient:

Silhouette Coefficient
----------------------

If the ground truth labels are not known, evaluation must be performed using
the model itself. The Silhouette Coefficient
(:func:`sklearn.metrics.silhouette_score`)
is an example of such an evaluation, where a
higher Silhouette Coefficient score relates to a model with better defined
clusters. The Silhouette Coefficient is defined for each sample and is composed
of two scores:

- **a**: The mean distance between a sample and all other points in the same
  class.

- **b**: The mean distance between a sample and all other points in the *next
  nearest cluster*.

The Silhouette Coefficient *s* for a single sample is then given as:

.. math:: s = \frac{b - a}{max(a, b)}

The Silhouette Coefficient for a set of samples is given as the mean of the
Silhouette Coefficient for each sample.


  >>> from sklearn import metrics
  >>> from sklearn.metrics import pairwise_distances
  >>> from sklearn import datasets
  >>> X, y = datasets.load_iris(return_X_y=True)

In normal usage, the Silhouette Coefficient is applied to the results of a
cluster analysis.

  >>> import numpy as np
  >>> from sklearn.cluster import KMeans
  >>> kmeans_model = KMeans(n_clusters=3, random_state=1).fit(X)
  >>> labels = kmeans_model.labels_
  >>> metrics.silhouette_score(X, labels, metric='euclidean')
  0.55...

.. topic:: References

 * Peter J. Rousseeuw (1987). :doi:`"Silhouettes: a Graphical Aid to the
   Interpretation and Validation of Cluster Analysis"<10.1016/0377-0427(87)90125-7>`
   . Computational and Applied Mathematics 20: 53–65.


Advantages
~~~~~~~~~~

- The score is bounded between -1 for incorrect clustering and +1 for highly
  dense clustering. Scores around zero indicate overlapping clusters.

- The score is higher when clusters are dense and well separated, which relates
  to a standard concept of a cluster.


Drawbacks
~~~~~~~~~

- The Silhouette Coefficient is generally higher for convex clusters than other
  concepts of clusters, such as density based clusters like those obtained
  through DBSCAN.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_kmeans_silhouette_analysis.py` : In this example
   the silhouette analysis is used to choose an optimal value for n_clusters.


.. _calinski_harabasz_index:

Calinski-Harabasz Index
-----------------------


If the ground truth labels are not known, the Calinski-Harabasz index
(:func:`sklearn.metrics.calinski_harabasz_score`) - also known as the Variance
Ratio Criterion - can be used to evaluate the model, where a higher
Calinski-Harabasz score relates to a model with better defined clusters.

The index is the ratio of the sum of between-clusters dispersion and of
within-cluster dispersion for all clusters (where dispersion is defined as the
sum of distances squared):

  >>> from sklearn import metrics
  >>> from sklearn.metrics import pairwise_distances
  >>> from sklearn import datasets
  >>> X, y = datasets.load_iris(return_X_y=True)

In normal usage, the Calinski-Harabasz index is applied to the results of a
cluster analysis:

  >>> import numpy as np
  >>> from sklearn.cluster import KMeans
  >>> kmeans_model = KMeans(n_clusters=3, random_state=1).fit(X)
  >>> labels = kmeans_model.labels_
  >>> metrics.calinski_harabasz_score(X, labels)
  561.62...

Advantages
~~~~~~~~~~

- The score is higher when clusters are dense and well separated, which relates
  to a standard concept of a cluster.

- The score is fast to compute.


Drawbacks
~~~~~~~~~

- The Calinski-Harabasz index is generally higher for convex clusters than other
  concepts of clusters, such as density based clusters like those obtained
  through DBSCAN.

Mathematical formulation
~~~~~~~~~~~~~~~~~~~~~~~~

For a set of data :math:`E` of size :math:`n_E` which has been clustered into
:math:`k` clusters, the Calinski-Harabasz score :math:`s` is defined as the
ratio of the between-clusters dispersion mean and the within-cluster dispersion:

.. math::
  s = \frac{\mathrm{tr}(B_k)}{\mathrm{tr}(W_k)} \times \frac{n_E - k}{k - 1}

where :math:`\mathrm{tr}(B_k)` is trace of the between group dispersion matrix
and :math:`\mathrm{tr}(W_k)` is the trace of the within-cluster dispersion
matrix defined by:

.. math:: W_k = \sum_{q=1}^k \sum_{x \in C_q} (x - c_q) (x - c_q)^T

.. math:: B_k = \sum_{q=1}^k n_q (c_q - c_E) (c_q - c_E)^T

with :math:`C_q` the set of points in cluster :math:`q`, :math:`c_q` the center
of cluster :math:`q`, :math:`c_E` the center of :math:`E`, and :math:`n_q` the
number of points in cluster :math:`q`.

.. topic:: References

 * Caliński, T., & Harabasz, J. (1974).
   `"A Dendrite Method for Cluster Analysis"
   <https://www.researchgate.net/publication/233096619_A_Dendrite_Method_for_Cluster_Analysis>`_.
   :doi:`Communications in Statistics-theory and Methods 3: 1-27 <10.1080/03610927408827101>`.


.. _davies-bouldin_index:

Davies-Bouldin Index
--------------------

If the ground truth labels are not known, the Davies-Bouldin index
(:func:`sklearn.metrics.davies_bouldin_score`) can be used to evaluate the
model, where a lower Davies-Bouldin index relates to a model with better
separation between the clusters.

This index signifies the average 'similarity' between clusters, where the
similarity is a measure that compares the distance between clusters with the
size of the clusters themselves.

Zero is the lowest possible score. Values closer to zero indicate a better
partition.

In normal usage, the Davies-Bouldin index is applied to the results of a
cluster analysis as follows:

  >>> from sklearn import datasets
  >>> iris = datasets.load_iris()
  >>> X = iris.data
  >>> from sklearn.cluster import KMeans
  >>> from sklearn.metrics import davies_bouldin_score
  >>> kmeans = KMeans(n_clusters=3, random_state=1).fit(X)
  >>> labels = kmeans.labels_
  >>> davies_bouldin_score(X, labels)
  0.6619...


Advantages
~~~~~~~~~~

- The computation of Davies-Bouldin is simpler than that of Silhouette scores.
- The index is solely based on quantities and features inherent to the dataset
  as its computation only uses point-wise distances.

Drawbacks
~~~~~~~~~

- The Davies-Boulding index is generally higher for convex clusters than other
  concepts of clusters, such as density based clusters like those obtained from
  DBSCAN.
- The usage of centroid distance limits the distance metric to Euclidean space.

Mathematical formulation
~~~~~~~~~~~~~~~~~~~~~~~~

The index is defined as the average similarity between each cluster :math:`C_i`
for :math:`i=1, ..., k` and its most similar one :math:`C_j`. In the context of
this index, similarity is defined as a measure :math:`R_{ij}` that trades off:

- :math:`s_i`, the average distance between each point of cluster :math:`i` and
  the centroid of that cluster -- also know as cluster diameter.
- :math:`d_{ij}`, the distance between cluster centroids :math:`i` and :math:`j`.

A simple choice to construct :math:`R_{ij}` so that it is nonnegative and
symmetric is:

.. math::
   R_{ij} = \frac{s_i + s_j}{d_{ij}}

Then the Davies-Bouldin index is defined as:

.. math::
   DB = \frac{1}{k} \sum_{i=1}^k \max_{i \neq j} R_{ij}


.. topic:: References

 * Davies, David L.; Bouldin, Donald W. (1979).
   :doi:`"A Cluster Separation Measure" <10.1109/TPAMI.1979.4766909>`
   IEEE Transactions on Pattern Analysis and Machine Intelligence.
   PAMI-1 (2): 224-227.

 * Halkidi, Maria; Batistakis, Yannis; Vazirgiannis, Michalis (2001).
   :doi:`"On Clustering Validation Techniques" <10.1023/A:1012801612483>`
   Journal of Intelligent Information Systems, 17(2-3), 107-145.

 * `Wikipedia entry for Davies-Bouldin index
   <https://en.wikipedia.org/wiki/Davies–Bouldin_index>`_.


.. _contingency_matrix:

Contingency Matrix
------------------

Contingency matrix (:func:`sklearn.metrics.cluster.contingency_matrix`)
reports the intersection cardinality for every true/predicted cluster pair.
The contingency matrix provides sufficient statistics for all clustering
metrics where the samples are independent and identically distributed and
one doesn't need to account for some instances not being clustered.

Here is an example::

   >>> from sklearn.metrics.cluster import contingency_matrix
   >>> x = ["a", "a", "a", "b", "b", "b"]
   >>> y = [0, 0, 1, 1, 2, 2]
   >>> contingency_matrix(x, y)
   array([[2, 1, 0],
          [0, 1, 2]])

The first row of output array indicates that there are three samples whose
true cluster is "a". Of them, two are in predicted cluster 0, one is in 1,
and none is in 2. And the second row indicates that there are three samples
whose true cluster is "b". Of them, none is in predicted cluster 0, one is in
1 and two are in 2.

A :ref:`confusion matrix <confusion_matrix>` for classification is a square
contingency matrix where the order of rows and columns correspond to a list
of classes.


Advantages
~~~~~~~~~~

- Allows to examine the spread of each true cluster across predicted
  clusters and vice versa.

- The contingency table calculated is typically utilized in the calculation
  of a similarity statistic (like the others listed in this document) between
  the two clusterings.

Drawbacks
~~~~~~~~~

- Contingency matrix is easy to interpret for a small number of clusters, but
  becomes very hard to interpret for a large number of clusters.

- It doesn't give a single metric to use as an objective for clustering
  optimisation.


.. topic:: References

 * `Wikipedia entry for contingency matrix
   <https://en.wikipedia.org/wiki/Contingency_table>`_

.. _pair_confusion_matrix:

Pair Confusion Matrix
---------------------

The pair confusion matrix
(:func:`sklearn.metrics.cluster.pair_confusion_matrix`) is a 2x2
similarity matrix

.. math::
   C = \left[\begin{matrix}
   C_{00} & C_{01} \\
   C_{10} & C_{11}
   \end{matrix}\right]

between two clusterings computed by considering all pairs of samples and
counting pairs that are assigned into the same or into different clusters
under the true and predicted clusterings.

It has the following entries:

  :math:`C_{00}` : number of pairs with both clusterings having the samples
  not clustered together

  :math:`C_{10}` : number of pairs with the true label clustering having the
  samples clustered together but the other clustering not having the samples
  clustered together

  :math:`C_{01}` : number of pairs with the true label clustering not having
  the samples clustered together but the other clustering having the samples
  clustered together

  :math:`C_{11}` : number of pairs with both clusterings having the samples
  clustered together

Considering a pair of samples that is clustered together a positive pair,
then as in binary classification the count of true negatives is
:math:`C_{00}`, false negatives is :math:`C_{10}`, true positives is
:math:`C_{11}` and false positives is :math:`C_{01}`.

Perfectly matching labelings have all non-zero entries on the
diagonal regardless of actual label values::

   >>> from sklearn.metrics.cluster import pair_confusion_matrix
   >>> pair_confusion_matrix([0, 0, 1, 1], [0, 0, 1, 1])
   array([[8, 0],
          [0, 4]])

::

   >>> pair_confusion_matrix([0, 0, 1, 1], [1, 1, 0, 0])
   array([[8, 0],
          [0, 4]])

Labelings that assign all classes members to the same clusters
are complete but may not always be pure, hence penalized, and
have some off-diagonal non-zero entries::

   >>> pair_confusion_matrix([0, 0, 1, 2], [0, 0, 1, 1])
   array([[8, 2],
          [0, 2]])

The matrix is not symmetric::

   >>> pair_confusion_matrix([0, 0, 1, 1], [0, 0, 1, 2])
   array([[8, 0],
          [2, 2]])

If classes members are completely split across different clusters, the
assignment is totally incomplete, hence the matrix has all zero
diagonal entries::

   >>> pair_confusion_matrix([0, 0, 0, 0], [0, 1, 2, 3])
   array([[ 0,  0],
          [12,  0]])

.. topic:: References

 * :doi:`"Comparing Partitions" <10.1007/BF01908075>`
   L. Hubert and P. Arabie, Journal of Classification 1985
