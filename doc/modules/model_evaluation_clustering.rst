.. currentmodule:: sklearn.metrics

.. _clustering_metrics:

==================
Clustering metrics
==================

Evaluating the performance of a clustering algorithm is not as trivial as
counting the number of errors or the precision and recall of a supervised
classification algorithm. In particular any evaluation metric should not
take the absolute values of the cluster labels into account but rather
if this clustering define separations of the data similar to some ground
truth set of classes or satisfying some assumption such that members
belong to the same class are more similar that members of different
classes according to some similarity metric.

.. _clustering_evaluation:

Clustering with known ground truth
==================================

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

 * `V-Measure: A conditional entropy-based external cluster evaluation
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


Clustering with unknown ground truth
====================================

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

.. _biclustering_evaluation:

Biclustering evaluation
=======================

There are two ways of evaluating a biclustering result: internal and
external. Internal measures, such as cluster stability, rely only on
the data and the result themselves. Currently there are no internal
bicluster measures in scikit-learn. External measures refer to an
external source of information, such as the true solution. When
working with real data the true solution is usually unknown, but
biclustering artificial data may be useful for evaluating algorithms
precisely because the true solution is known.

To compare a set of found biclusters to the set of true biclusters,
two similarity measures are needed: a similarity measure for
individual biclusters, and a way to combine these individual
similarities into an overall score.

To compare individual biclusters, several measures have been used. For
now, only the Jaccard index is implemented:

.. math::
    J(A, B) = \frac{|A \cap B|}{|A| + |B| - |A \cap B|}

where :math:`A` and :math:`B` are biclusters, :math:`|A \cap B|` is
the number of elements in their intersection. The Jaccard index
achieves its minimum of 0 when the biclusters to not overlap at all
and its maximum of 1 when they are identical.

Several methods have been developed to compare two sets of biclusters.
For now, only :func:`consensus_score` (Hochreiter et. al., 2010) is
available:

1. Compute bicluster similarities for pairs of biclusters, one in each
   set, using the Jaccard index or a similar measure.

2. Assign biclusters from one set to another in a one-to-one fashion
   to maximize the sum of their similarities. This step is performed
   using the Hungarian algorithm.

3. The final sum of similarities is divided by the size of the larger
   set.

The minimum consensus score, 0, occurs when all pairs of biclusters
are totally dissimilar. The maximum score, 1, occurs when both sets
are identical.


.. topic:: References:

 * Hochreiter, Bodenhofer, et. al., 2010. `FABIA: factor analysis
   for bicluster acquisition
   <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881408/>`__.
