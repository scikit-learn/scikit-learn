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
    methods accept standard data matrices of shape ``[n_samples, n_features]``.
    These can be obtained from the classes in the :mod:`sklearn.feature_extraction`
    module. For :class:`AffinityPropagation`, :class:`SpectralClustering`
    and :class:`DBSCAN` one can also input similarity matrices of shape
    ``[n_samples, n_samples]``. These can be obtained from the functions
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
     - General-purpose, even cluster size, flat geometry, not too many clusters
     - Distances between points

   * - :ref:`Affinity propagation <affinity_propagation>`
     - damping, sample preference
     - Not scalable with n_samples
     - Many clusters, uneven cluster size, non-flat geometry
     - Graph distance (e.g. nearest-neighbor graph)

   * - :ref:`Mean-shift <mean_shift>`
     - bandwidth
     - Not scalable with ``n_samples``
     - Many clusters, uneven cluster size, non-flat geometry
     - Distances between points

   * - :ref:`Spectral clustering <spectral_clustering>`
     - number of clusters
     - Medium ``n_samples``, small ``n_clusters``
     - Few clusters, even cluster size, non-flat geometry
     - Graph distance (e.g. nearest-neighbor graph)

   * - :ref:`Ward hierarchical clustering <hierarchical_clustering>`
     - number of clusters
     - Large ``n_samples`` and ``n_clusters``
     - Many clusters, possibly connectivity constraints
     - Distances between points

   * - :ref:`Agglomerative clustering <hierarchical_clustering>`
     - number of clusters, linkage type, distance
     - Large ``n_samples`` and ``n_clusters``
     - Many clusters, possibly connectivity constraints, non Euclidean
       distances
     - Any pairwise distance

   * - :ref:`DBSCAN <dbscan>`
     - neighborhood size
     - Very large ``n_samples``, medium ``n_clusters``
     - Non-flat geometry, uneven cluster sizes
     - Distances between nearest points

   * - :ref:`Gaussian mixtures <mixture>`
     - many
     - Not scalable
     - Flat geometry, good for density estimation
     - Mahalanobis distances to  centers

   * - :ref:`Birch`
     - branching factor, threshold, optional global clusterer.
     - Large ``n_clusters`` and ``n_samples``
     - Large dataset, outlier removal, data reduction.
     - Euclidean distance between points

Non-flat geometry clustering is useful when the clusters have a specific
shape, i.e. a non-flat manifold, and the standard euclidean distance is
not the right metric. This case arises in the two top rows of the figure
above.

Gaussian mixture models, useful for clustering, are described in
:ref:`another chapter of the documentation <mixture>` dedicated to
mixture models. KMeans can be seen as a special case of Gaussian mixture
model with equal covariance per component.

.. _k_means:

K-means
=======

The :class:`KMeans` algorithm clusters data by trying to separate samples
in n groups of equal variance, minimizing a criterion known as the
`inertia <inertia>`_ or within-cluster sum-of-squares.
This algorithm requires the number of clusters to be specified.
It scales well to large number of samples and has been used
across a large range of application areas in many different fields.

The k-means algorithm divides a set of :math:`N` samples :math:`X`
into :math:`K` disjoint clusters :math:`C`,
each described by the mean :math:`\mu_j` of the samples in the cluster.
The means are commonly called the cluster "centroids";
note that they are not, in general, points from :math:`X`,
although they live in the same space.
The K-means algorithm aims to choose centroids
that minimise the *inertia*, or within-cluster sum of squared criterion:

.. math:: \sum_{i=0}^{n}\min_{\mu_j \in C}(||x_j - \mu_i||^2)

Inertia, or the within-cluster sum of squares criterion,
can be recognized as a measure of how internally coherent clusters are.
It suffers from various drawbacks:

- Inertia makes the assumption that clusters are convex and isotropic,
  which is not always the case. It responds poorly to elongated clusters,
  or manifolds with irregular shapes.

- Inertia is not a normalized metric: we just know that lower values are
  better and zero is optimal. But in very high-dimensional spaces, Euclidean
  distances tend to become inflated
  (this is an instance of the so-called "curse of dimensionality").
  Running a dimensionality reduction algorithm such as `PCA <PCA>`_
  prior to k-means clustering can alleviate this problem
  and speed up the computations.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_kmeans_assumptions_001.png
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
(use the ``init='kmeans++'`` parameter). This initializes the centroids to be
(generally) distant from each other, leading to provably better results than
random initialization, as shown in the reference.

A parameter can be given to allow K-means to be run in parallel, called
``n_jobs``. Giving this parameter a positive value uses that many processors
(default: 1). A value of -1 uses all available processors, with -2 using one
less, and so on. Parallelization generally speeds up computation at the cost of
memory (in this case, multiple copies of centroids need to be stored, one for
each job).

.. warning::

    The parallel version of K-Means is broken on OS X when `numpy` uses the
    `Accelerate` Framework. This is expected behavior: `Accelerate` can be called
    after a fork but you need to execv the subprocess with the Python binary
    (which multiprocessing does not do under posix).

K-means can be used for vector quantization. This is achieved using the
transform method of a trained model of :class:`KMeans`.

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

 * :ref:`sphx_glr_auto_examples_text_document_clustering.py`: Document clustering using sparse
   MiniBatchKMeans

 * :ref:`sphx_glr_auto_examples_cluster_plot_dict_face_patches.py`


.. topic:: References:

 * `"Web Scale K-Means clustering"
   <http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf>`_
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
used, and the *damping factor*.

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

    r(i, k) \leftarrow s(i, k) - max [ a(i, \acute{k}) + s(i, \acute{k}) \forall \acute{k} \neq k ]

Where :math:`s(i, k)` is the similarity between samples :math:`i` and :math:`k`.
The availability of sample :math:`k`
to be the exemplar of sample :math:`i` is given by:

.. math::

    a(i, k) \leftarrow min [0, r(k, k) + \sum_{\acute{i}~s.t.~\acute{i} \notin \{i, k\}}{r(\acute{i}, k)}]

To begin with, all values for :math:`r` and :math:`a` are set to zero,
and the calculation of each iterates until convergence.

.. _mean_shift:

Mean Shift
==========
:class:`MeanShift` clustering aims to discover *blobs* in a smooth density of
samples. It is a centroid based algorithm, which works by updating candidates
for centroids to be the mean of the points within a given region. These
candidates are then filtered in a post-processing stage to eliminate
near-duplicates to form the final set of centroids.

Given a candidate centroid :math:`x_i` for iteration :math:`t`, the candidate
is updated according to the following equation:

.. math::

    x_i^{t+1} = x_i^t + m(x_i^t)

Where :math:`N(x_i)` is the neighborhood of samples within a given distance
around :math:`x_i` and :math:`m` is the *mean shift* vector that is computed for each
centroid that points towards a region of the maximum increase in the density of points.
This is computed using the following equation, effectively updating a centroid
to be the mean of the samples within its neighborhood:

.. math::

    m(x_i) = \frac{\sum_{x_j \in N(x_i)}K(x_j - x_i)x_j}{\sum_{x_j \in N(x_i)}K(x_j - x_i)}

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

 * `"Mean shift: A robust approach toward feature space analysis."
   <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.76.8968&rep=rep1&type=pdf>`_
   D. Comaniciu and P. Meer, *IEEE Transactions on Pattern Analysis and Machine Intelligence* (2002)


.. _spectral_clustering:

Spectral clustering
===================

:class:`SpectralClustering` does a low-dimension embedding of the
affinity matrix between samples, followed by a KMeans in the low
dimensional space. It is especially efficient if the affinity matrix is
sparse and the `pyamg <http://pyamg.org/>`_ module is installed.
SpectralClustering requires the number of clusters to be specified. It
works well for a small number of clusters but is not advised when using
many clusters.

For two clusters, it solves a convex relaxation of the `normalised
cuts <http://people.eecs.berkeley.edu/~malik/papers/SM-ncut.pdf>`_ problem on
the similarity graph: cutting the graph in two so that the weight of the
edges cut is small compared to the weights of the edges inside each
cluster. This criteria is especially interesting when working on images:
graph vertices are pixels, and edges of the similarity graph are a
function of the gradient of the image.


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

 * :ref:`sphx_glr_auto_examples_cluster_plot_face_segmentation.py`: Spectral clustering
   to split the image of the raccoon face in regions.

.. |face_kmeans| image:: ../auto_examples/cluster/images/sphx_glr_plot_face_segmentation_001.png
    :target: ../auto_examples/cluster/plot_face_segmentation.html
    :scale: 65

.. |face_discretize| image:: ../auto_examples/cluster/images/sphx_glr_plot_face_segmentation_002.png
    :target: ../auto_examples/cluster/plot_face_segmentation.html
    :scale: 65

Different label assignment strategies
---------------------------------------

Different label assignment strategies can be used, corresponding to the
``assign_labels`` parameter of :class:`SpectralClustering`.
The ``"kmeans"`` strategy can match finer details of the data, but it can be
more unstable. In particular, unless you control the ``random_state``, it
may not be reproducible from run-to-run, as it depends on a random
initialization. On the other hand, the ``"discretize"`` strategy is 100%
reproducible, but it tends to create parcels of fairly even and
geometrical shape.

=====================================  =====================================
 ``assign_labels="kmeans"``              ``assign_labels="discretize"``
=====================================  =====================================
|face_kmeans|                          |face_discretize|
=====================================  =====================================


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

:class:`AgglomerativeClustering` can also scale to large number of samples
when it is used jointly with a connectivity matrix, but is computationally
expensive when no connectivity constraints are added between samples: it
considers at each step all the possible merges.

.. topic:: :class:`FeatureAgglomeration`

   The :class:`FeatureAgglomeration` uses agglomerative clustering to
   group together features that look very similar, thus decreasing the
   number of features. It is a dimensionality reduction tool, see
   :ref:`data_reduction`.

Different linkage type: Ward, complete and average linkage
-----------------------------------------------------------

:class:`AgglomerativeClustering` supports Ward, average, and complete
linkage strategies.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_digits_linkage_001.png
    :target: ../auto_examples/cluster/plot_digits_linkage.html
    :scale: 43

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_digits_linkage_002.png
    :target: ../auto_examples/cluster/plot_digits_linkage.html
    :scale: 43

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_digits_linkage_003.png
    :target: ../auto_examples/cluster/plot_digits_linkage.html
    :scale: 43


Agglomerative cluster has a "rich get richer" behavior that leads to
uneven cluster sizes. In this regard, complete linkage is the worst
strategy, and Ward gives the most regular sizes. However, the affinity
(or distance used in clustering) cannot be varied with Ward, thus for non
Euclidean metrics, average linkage is a good alternative.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_digits_linkage.py`: exploration of the
   different linkage strategies in a real dataset.


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
:ref:`raccoon face <sphx_glr_auto_examples_cluster_plot_face_ward_segmentation.py>` example.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_face_ward_segmentation.py`: Ward clustering
   to split the image of a raccoon face in regions.

 * :ref:`sphx_glr_auto_examples_cluster_plot_ward_structured_vs_unstructured.py`: Example of
   Ward algorithm on a swiss-roll, comparison of structured approaches
   versus unstructured approaches.

 * :ref:`sphx_glr_auto_examples_cluster_plot_feature_agglomeration_vs_univariate_selection.py`:
   Example of dimensionality reduction with feature agglomeration based on
   Ward hierarchical clustering.

 * :ref:`sphx_glr_auto_examples_cluster_plot_agglomerative_clustering.py`

.. warning:: **Connectivity constraints with average and complete linkage**

    Connectivity constraints and complete or average linkage can enhance
    the 'rich getting richer' aspect of agglomerative clustering,
    particularly so if they are built with
    :func:`sklearn.neighbors.kneighbors_graph`. In the limit of a small
    number of clusters, they tend to give a few macroscopically occupied
    clusters and almost empty ones. (see the discussion in
    :ref:`sphx_glr_auto_examples_cluster_plot_agglomerative_clustering.py`).

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

Average and complete linkage can be used with a variety of distances (or
affinities), in particular Euclidean distance (*l2*), Manhattan distance
(or Cityblock, or *l1*), cosine distance, or any precomputed affinity
matrix.

* *l1* distance is often good for sparse features, or sparse noise: ie
  many of the features are zero, as in text mining using occurences of
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

In the figure below, the color indicates cluster membership, with large circles
indicating core samples found by the algorithm. Smaller circles are non-core
samples that are still part of a cluster. Moreover, the outliers are indicated
by black points below.

.. |dbscan_results| image:: ../auto_examples/cluster/images/sphx_glr_plot_dbscan_001.png
        :target: ../auto_examples/cluster/plot_dbscan.html
        :scale: 50

.. centered:: |dbscan_results|

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_cluster_plot_dbscan.py`

.. topic:: Implementation

    The algorithm is non-deterministic, but the core samples will
    always belong to the same clusters (although the labels may be
    different). The non-determinism comes from deciding to which cluster a
    non-core sample belongs. A non-core sample can have a distance lower
    than ``eps`` to two core samples in different clusters. By the
    triangular inequality, those two core samples must be more distant than
    ``eps`` from each other, or they would be in the same cluster. The non-core
    sample is assigned to whichever cluster is generated first, where
    the order is determined randomly. Other than the ordering of
    the dataset, the algorithm is deterministic, making the results relatively
    stable between runs on the same data.

    The current implementation uses ball trees and kd-trees
    to determine the neighborhood of points,
    which avoids calculating the full distance matrix
    (as was done in scikit-learn versions before 0.14).
    The possibility to use custom metrics is retained;
    for details, see :class:`NearestNeighbors`.

.. topic:: Memory consumption for large sample sizes

    This implementation is by default not memory efficient because it constructs
    a full pairwise similarity matrix in the case where kd-trees or ball-trees cannot
    be used (e.g. with sparse matrices). This matrix will consume n^2 floats.
    A couple of mechanisms for getting around this are:

    - A sparse radius neighborhood graph (where missing
      entries are presumed to be out of eps) can be precomputed in a memory-efficient
      way and dbscan can be run over this with ``metric='precomputed'``.

    - The dataset can be compressed, either by removing exact duplicates if
      these occur in your data, or by using BIRCH. Then you only have a
      relatively small number of representatives for a large number of points.
      You can then provide a ``sample_weight`` when fitting DBSCAN.

.. topic:: References:

 * "A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases
   with Noise"
   Ester, M., H. P. Kriegel, J. Sander, and X. Xu,
   In Proceedings of the 2nd International Conference on Knowledge Discovery
   and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996

.. _birch:

Birch
=====

The :class:`Birch` builds a tree called the Characteristic Feature Tree (CFT)
for the given data. The data is essentially lossy compressed to a set of
Characteristic Feature nodes (CF Nodes). The CF Nodes have a number of
subclusters called Characteristic Feature subclusters (CF Subclusters)
and these CF Subclusters located in the non-terminal CF Nodes
can have CF Nodes as children.

The CF Subclusters hold the necessary information for clustering which prevents
the need to hold the entire input data in memory. This information includes:

- Number of samples in a subcluster.
- Linear Sum - A n-dimensional vector holding the sum of all samples
- Squared Sum - Sum of the squared L2 norm of all samples.
- Centroids - To avoid recalculation linear sum / n_samples.
- Squared norm of the centroids.

The Birch algorithm has two parameters, the threshold and the branching factor.
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

**Birch or MiniBatchKMeans?**

 - Birch does not scale very well to high dimensional data. As a rule of thumb if
   ``n_features`` is greater than twenty, it is generally better to use MiniBatchKMeans.
 - If the number of instances of data needs to be reduced, or if one wants a
   large number of subclusters either as a preprocessing step or otherwise,
   Birch is more useful than MiniBatchKMeans.


**How to use partial_fit?**

To avoid the computation of global clustering, for every call of ``partial_fit``
the user is advised

 1. To set ``n_clusters=None`` initially
 2. Train all data by multiple calls to partial_fit.
 3. Set ``n_clusters`` to a required value using
    ``brc.set_params(n_clusters=n_clusters)``.
 4. Call ``partial_fit`` finally with no arguments, i.e ``brc.partial_fit()``
    which performs the global clustering.

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_birch_vs_minibatchkmeans_001.png
    :target: ../auto_examples/cluster/plot_birch_vs_minibatchkmeans.html

.. topic:: References:

 * Tian Zhang, Raghu Ramakrishnan, Maron Livny
   BIRCH: An efficient data clustering method for large databases.
   http://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf

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
belong to the same class are more similar that members of different
classes according to some similarity metric.

.. currentmodule:: sklearn.metrics

.. _adjusted_rand_score:

Adjusted Rand index
-------------------

Given the knowledge of the ground truth class assignments ``labels_true``
and our clustering algorithm assignments of the same samples
``labels_pred``, the **adjusted Rand index** is a function that measures
the **similarity** of the two assignments, ignoring permutations and **with
chance normalization**::

  >>> from sklearn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]

  >>> metrics.adjusted_rand_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.24...

One can permute 0 and 1 in the predicted labels, rename 2 to 3, and get
the same score::

  >>> labels_pred = [1, 1, 0, 0, 3, 3]
  >>> metrics.adjusted_rand_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.24...

Furthermore, :func:`adjusted_rand_score` is **symmetric**: swapping the argument
does not change the score. It can thus be used as a **consensus
measure**::

  >>> metrics.adjusted_rand_score(labels_pred, labels_true)  # doctest: +ELLIPSIS
  0.24...

Perfect labeling is scored 1.0::

  >>> labels_pred = labels_true[:]
  >>> metrics.adjusted_rand_score(labels_true, labels_pred)
  1.0

Bad (e.g. independent labelings) have negative or close to 0.0 scores::

  >>> labels_true = [0, 1, 2, 0, 3, 4, 5, 1]
  >>> labels_pred = [1, 1, 0, 0, 2, 2, 2, 2]
  >>> metrics.adjusted_rand_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  -0.12...


Advantages
~~~~~~~~~~

- **Random (uniform) label assignments have a ARI score close to 0.0**
  for any value of ``n_clusters`` and ``n_samples`` (which is not the
  case for raw Rand index or the V-measure for instance).

- **Bounded range [-1, 1]**: negative values are bad (independent
  labelings), similar clusterings have a positive ARI, 1.0 is the perfect
  match score.

- **No assumption is made on the cluster structure**: can be used
  to compare clustering algorithms such as k-means which assumes isotropic
  blob shapes with results of spectral clustering algorithms which can
  find cluster with "folded" shapes.


Drawbacks
~~~~~~~~~

- Contrary to inertia, **ARI requires knowledge of the ground truth
  classes** while is almost never available in practice or requires manual
  assignment by human annotators (as in the supervised learning setting).

  However ARI can also be useful in a purely unsupervised setting as a
  building block for a Consensus Index that can be used for clustering
  model selection (TODO).


.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_cluster_plot_adjusted_for_chance_measures.py`: Analysis of
   the impact of the dataset size on the value of clustering measures
   for random assignments.


Mathematical formulation
~~~~~~~~~~~~~~~~~~~~~~~~

If C is a ground truth class assignment and K the clustering, let us
define :math:`a` and :math:`b` as:

- :math:`a`, the number of pairs of elements that are in the same set
  in C and in the same set in K

- :math:`b`, the number of pairs of elements that are in different sets
  in C and in different sets in K

The raw (unadjusted) Rand index is then given by:

.. math:: \text{RI} = \frac{a + b}{C_2^{n_{samples}}}

Where :math:`C_2^{n_{samples}}` is the total number of possible pairs
in the dataset (without ordering).

However the RI score does not guarantee that random label assignments
will get a value close to zero (esp. if the number of clusters is in
the same order of magnitude as the number of samples).

To counter this effect we can discount the expected RI :math:`E[\text{RI}]` of
random labelings by defining the adjusted Rand index as follows:

.. math:: \text{ARI} = \frac{\text{RI} - E[\text{RI}]}{\max(\text{RI}) - E[\text{RI}]}

.. topic:: References

 * `Comparing Partitions
   <http://link.springer.com/article/10.1007%2FBF01908075>`_
   L. Hubert and P. Arabie, Journal of Classification 1985

 * `Wikipedia entry for the adjusted Rand index
   <https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index>`_

.. _mutual_info_score:

Mutual Information based scores
-------------------------------

Given the knowledge of the ground truth class assignments ``labels_true`` and
our clustering algorithm assignments of the same samples ``labels_pred``, the
**Mutual Information** is a function that measures the **agreement** of the two
assignments, ignoring permutations.  Two different normalized versions of this
measure are available, **Normalized Mutual Information(NMI)** and **Adjusted
Mutual Information(AMI)**. NMI is often used in the literature while AMI was
proposed more recently and is **normalized against chance**::

  >>> from sklearn import metrics
  >>> labels_true = [0, 0, 0, 1, 1, 1]
  >>> labels_pred = [0, 0, 1, 1, 2, 2]

  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.22504...

One can permute 0 and 1 in the predicted labels, rename 2 to 3 and get
the same score::

  >>> labels_pred = [1, 1, 0, 0, 3, 3]
  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.22504...

All, :func:`mutual_info_score`, :func:`adjusted_mutual_info_score` and
:func:`normalized_mutual_info_score` are symmetric: swapping the argument does
not change the score. Thus they can be used as a **consensus measure**::

  >>> metrics.adjusted_mutual_info_score(labels_pred, labels_true)  # doctest: +ELLIPSIS
  0.22504...

Perfect labeling is scored 1.0::

  >>> labels_pred = labels_true[:]
  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)
  1.0

  >>> metrics.normalized_mutual_info_score(labels_true, labels_pred)
  1.0

This is not true for ``mutual_info_score``, which is therefore harder to judge::

  >>> metrics.mutual_info_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.69...

Bad (e.g. independent labelings) have non-positive scores::

  >>> labels_true = [0, 1, 2, 0, 3, 4, 5, 1]
  >>> labels_pred = [1, 1, 0, 0, 2, 2, 2, 2]
  >>> metrics.adjusted_mutual_info_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  -0.10526...


Advantages
~~~~~~~~~~

- **Random (uniform) label assignments have a AMI score close to 0.0**
  for any value of ``n_clusters`` and ``n_samples`` (which is not the
  case for raw Mutual Information or the V-measure for instance).

- **Bounded range [0, 1]**:  Values close to zero indicate two label
  assignments that are largely independent, while values close to one
  indicate significant agreement. Further, values of exactly 0 indicate
  **purely** independent label assignments and a AMI of exactly 1 indicates
  that the two label assignments are equal (with or without permutation).

- **No assumption is made on the cluster structure**: can be used
  to compare clustering algorithms such as k-means which assumes isotropic
  blob shapes with results of spectral clustering algorithms which can
  find cluster with "folded" shapes.


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

.. math:: H(U) = \sum_{i=1}^{|U|}P(i)\log(P(i))

where :math:`P(i) = |U_i| / N` is the probability that an object picked at
random from :math:`U` falls into class :math:`U_i`. Likewise for :math:`V`:

.. math:: H(V) = \sum_{j=1}^{|V|}P'(j)\log(P'(j))

With :math:`P'(j) = |V_j| / N`. The mutual information (MI) between :math:`U`
and :math:`V` is calculated by:

.. math:: \text{MI}(U, V) = \sum_{i=1}^{|U|}\sum_{j=1}^{|V|}P(i, j)\log\left(\frac{P(i,j)}{P(i)P'(j)}\right)

where :math:`P(i, j) = |U_i \cap V_j| / N` is the probability that an object
picked at random falls into both classes :math:`U_i` and :math:`V_j`.

The normalized mutual information is defined as

.. math:: \text{NMI}(U, V) = \frac{\text{MI}(U, V)}{\sqrt{H(U)H(V)}}

This value of the mutual information and also the normalized variant is not
adjusted for chance and will tend to increase as the number of different labels
(clusters) increases, regardless of the actual amount of "mutual information"
between the label assignments.

The expected value for the mutual information can be calculated using the
following equation, from Vinh, Epps, and Bailey, (2009). In this equation,
:math:`a_i = |U_i|` (the number of elements in :math:`U_i`) and
:math:`b_j = |V_j|` (the number of elements in :math:`V_j`).


.. math:: E[\text{MI}(U,V)]=\sum_{i=1}^|U| \sum_{j=1}^|V| \sum_{n_{ij}=(a_i+b_j-N)^+
   }^{\min(a_i, b_j)} \frac{n_{ij}}{N}\log \left( \frac{ N.n_{ij}}{a_i b_j}\right)
   \frac{a_i!b_j!(N-a_i)!(N-b_j)!}{N!n_{ij}!(a_i-n_{ij})!(b_j-n_{ij})!
   (N-a_i-b_j+n_{ij})!}

Using the expected value, the adjusted mutual information can then be
calculated using a similar form to that of the adjusted Rand index:

.. math:: \text{AMI} = \frac{\text{MI} - E[\text{MI}]}{\max(H(U), H(V)) - E[\text{MI}]}

.. topic:: References

 * Strehl, Alexander, and Joydeep Ghosh (2002). "Cluster ensembles – a
   knowledge reuse framework for combining multiple partitions". Journal of
   Machine Learning Research 3: 583–617.
   `doi:10.1162/153244303321897735 <http://strehl.com/download/strehl-jmlr02.pdf>`_.

 * Vinh, Epps, and Bailey, (2009). "Information theoretic measures
   for clusterings comparison". Proceedings of the 26th Annual International
   Conference on Machine Learning - ICML '09.
   `doi:10.1145/1553374.1553511 <https://dl.acm.org/citation.cfm?doid=1553374.1553511>`_.
   ISBN 9781605585161.

 * Vinh, Epps, and Bailey, (2010). Information Theoretic Measures for
   Clusterings Comparison: Variants, Properties, Normalization and
   Correction for Chance, JMLR
   http://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf

 * `Wikipedia entry for the (normalized) Mutual Information
   <https://en.wikipedia.org/wiki/Mutual_Information>`_

 * `Wikipedia entry for the Adjusted Mutual Information
   <https://en.wikipedia.org/wiki/Adjusted_Mutual_Information>`_

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

  >>> metrics.homogeneity_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.66...

  >>> metrics.completeness_score(labels_true, labels_pred) # doctest: +ELLIPSIS
  0.42...

Their harmonic mean called **V-measure** is computed by
:func:`v_measure_score`::

  >>> metrics.v_measure_score(labels_true, labels_pred)    # doctest: +ELLIPSIS
  0.51...

The V-measure is actually equivalent to the mutual information (NMI)
discussed above normalized by the sum of the label entropies [B2011]_.

Homogeneity, completeness and V-measure can be computed at once using
:func:`homogeneity_completeness_v_measure` as follows::

  >>> metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
  ...                                                      # doctest: +ELLIPSIS
  (0.66..., 0.42..., 0.51...)

The following clustering assignment is slightly better, since it is
homogeneous but not complete::

  >>> labels_pred = [0, 0, 0, 1, 2, 2]
  >>> metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
  ...                                                      # doctest: +ELLIPSIS
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

 .. [RH2007] `V-Measure: A conditional entropy-based external cluster evaluation
   measure <http://aclweb.org/anthology/D/D07/D07-1043.pdf>`_
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

  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.47140...

One can permute 0 and 1 in the predicted labels, rename 2 to 3 and get
the same score::

  >>> labels_pred = [1, 1, 0, 0, 3, 3]

  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.47140...

Perfect labeling is scored 1.0::

  >>> labels_pred = labels_true[:]
  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  1.0

Bad (e.g. independent labelings) have zero scores::

  >>> labels_true = [0, 1, 2, 0, 3, 4, 5, 1]
  >>> labels_pred = [1, 1, 0, 0, 2, 2, 2, 2]
  >>> metrics.fowlkes_mallows_score(labels_true, labels_pred)  # doctest: +ELLIPSIS
  0.0

Advantages
~~~~~~~~~~

- **Random (uniform) label assignments have a FMI score close to 0.0**
  for any value of ``n_clusters`` and ``n_samples`` (which is not the
  case for raw Mutual Information or the V-measure for instance).

- **Bounded range [0, 1]**:  Values close to zero indicate two label
  assignments that are largely independent, while values close to one
  indicate significant agreement. Further, values of exactly 0 indicate
  **purely** independent label assignments and a AMI of exactly 1 indicates
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
    http://wildfire.stat.ucla.edu/pdflibrary/fowlkes.pdf

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
  >>> dataset = datasets.load_iris()
  >>> X = dataset.data
  >>> y = dataset.target

In normal usage, the Silhouette Coefficient is applied to the results of a
cluster analysis.

  >>> import numpy as np
  >>> from sklearn.cluster import KMeans
  >>> kmeans_model = KMeans(n_clusters=3, random_state=1).fit(X)
  >>> labels = kmeans_model.labels_
  >>> metrics.silhouette_score(X, labels, metric='euclidean')
  ...                                                      # doctest: +ELLIPSIS
  0.55...

.. topic:: References

 * Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
   Interpretation and Validation of Cluster Analysis". Computational
   and Applied Mathematics 20: 53–65.
   `doi:10.1016/0377-0427(87)90125-7 <http://dx.doi.org/10.1016/0377-0427(87)90125-7>`_.


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

.. _calinski_harabaz_index:

Calinski-Harabaz Index
----------------------

If the ground truth labels are not known, the Calinski-Harabaz index
(:func:`sklearn.metrics.calinski_harabaz_score`) can be used to evaluate the
model, where a higher Calinski-Harabaz score relates to a model with better
defined clusters.

For :math:`k` clusters, the Calinski-Harabaz score :math:`s` is given as the
ratio of the between-clusters dispersion mean and the within-cluster
dispersion:

.. math::
  s(k) = \frac{\mathrm{Tr}(B_k)}{\mathrm{Tr}(W_k)} \times \frac{N - k}{k - 1}

where :math:`B_K` is the between group dispersion matrix and :math:`W_K`
is the within-cluster dispersion matrix defined by:

.. math:: W_k = \sum_{q=1}^k \sum_{x \in C_q} (x - c_q) (x - c_q)^T

.. math:: B_k = \sum_q n_q (c_q - c) (c_q - c)^T

with :math:`N` be the number of points in our data, :math:`C_q` be the set of
points in cluster :math:`q`, :math:`c_q` be the center of cluster
:math:`q`, :math:`c` be the center of :math:`E`, :math:`n_q` be the number of
points in cluster :math:`q`.


  >>> from sklearn import metrics
  >>> from sklearn.metrics import pairwise_distances
  >>> from sklearn import datasets
  >>> dataset = datasets.load_iris()
  >>> X = dataset.data
  >>> y = dataset.target

In normal usage, the Calinski-Harabaz index is applied to the results of a
cluster analysis.

  >>> import numpy as np
  >>> from sklearn.cluster import KMeans
  >>> kmeans_model = KMeans(n_clusters=3, random_state=1).fit(X)
  >>> labels = kmeans_model.labels_
  >>> metrics.calinski_harabaz_score(X, labels)  # doctest: +ELLIPSIS
  560.39...


Advantages
~~~~~~~~~~

- The score is higher when clusters are dense and well separated, which relates
  to a standard concept of a cluster.

- The score is fast to compute


Drawbacks
~~~~~~~~~

- The Calinski-Harabaz index is generally higher for convex clusters than other
  concepts of clusters, such as density based clusters like those obtained
  through DBSCAN.

.. topic:: References

 *  Caliński, T., & Harabasz, J. (1974). "A dendrite method for cluster
    analysis". Communications in Statistics-theory and Methods 3: 1-27.
    `doi:10.1080/03610926.2011.560741 <http://dx.doi.org/10.1080/03610926.2011.560741>`_.
