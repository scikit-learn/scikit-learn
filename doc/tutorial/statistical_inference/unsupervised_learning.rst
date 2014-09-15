============================================================
Unsupervised learning: seeking representations of the data
============================================================

Clustering: grouping observations together
============================================

.. topic:: The problem solved in clustering

    Given the iris dataset, if we knew that there were 3 types of iris, but
    did not have access to a taxonomist to label them: we could try a
    **clustering task**: split the observations into well-separated group
    called *clusters*.

..   
   >>> # Set the PRNG   
   >>> import numpy as np
   >>> np.random.seed(1)

K-means clustering
-------------------

Note that there exist a lot of different clustering criteria and associated
algorithms. The simplest clustering algorithm is 
:ref:`k_means`.

.. image:: ../../auto_examples/cluster/images/plot_cluster_iris_002.png
    :target: ../../auto_examples/cluster/plot_cluster_iris.html
    :scale: 70
    :align: right


:: 

    >>> from sklearn import cluster, datasets
    >>> iris = datasets.load_iris()
    >>> X_iris = iris.data
    >>> y_iris = iris.target

    >>> k_means = cluster.KMeans(n_clusters=3)
    >>> k_means.fit(X_iris) # doctest: +ELLIPSIS
    KMeans(copy_x=True, init='k-means++', ...
    >>> print(k_means.labels_[::10])
    [1 1 1 1 1 0 0 0 0 0 2 2 2 2 2]
    >>> print(y_iris[::10])
    [0 0 0 0 0 1 1 1 1 1 2 2 2 2 2]

.. |k_means_iris_bad_init| image:: ../../auto_examples/cluster/images/plot_cluster_iris_003.png
   :target: ../../auto_examples/cluster/plot_cluster_iris.html
   :scale: 63

.. |k_means_iris_8| image:: ../../auto_examples/cluster/images/plot_cluster_iris_001.png
   :target: ../../auto_examples/cluster/plot_cluster_iris.html
   :scale: 63

.. |cluster_iris_truth| image:: ../../auto_examples/cluster/images/plot_cluster_iris_004.png
   :target: ../../auto_examples/cluster/plot_cluster_iris.html
   :scale: 63

.. warning:: 
   
    There is absolutely no guarantee of recovering a ground truth. First,
    choosing the right number of clusters is hard. Second, the algorithm
    is sensitive to initialization, and can fall into local minima,
    although scikit-learn employs several tricks to mitigate this issue.

    .. list-table::
        :class: centered
        
        * 
        
            - |k_means_iris_bad_init|

            - |k_means_iris_8|

            - |cluster_iris_truth|

        * 
        
            - **Bad initialization**
            
            - **8 clusters**
            
            - **Ground truth**

    **Don't over-interpret clustering results**

.. |lena| image:: ../../auto_examples/cluster/images/plot_lena_compress_001.png
   :target: ../../auto_examples/cluster/plot_lena_compress.html
   :scale: 60

.. |lena_regular| image:: ../../auto_examples/cluster/images/plot_lena_compress_002.png
   :target: ../../auto_examples/cluster/plot_lena_compress.html
   :scale: 60

.. |lena_compressed| image:: ../../auto_examples/cluster/images/plot_lena_compress_003.png
   :target: ../../auto_examples/cluster/plot_lena_compress.html
   :scale: 60

.. |lena_histogram| image:: ../../auto_examples/cluster/images/plot_lena_compress_004.png
   :target: ../../auto_examples/cluster/plot_lena_compress.html
   :scale: 60

.. topic:: **Application example: vector quantization**

    Clustering in general and KMeans, in particular, can be seen as a way
    of choosing a small number of exemplars to compress the information.
    The problem is sometimes known as 
    `vector quantization <http://en.wikipedia.org/wiki/Vector_quantization>`_. 
    For instance, this can be used to posterize an image::

        >>> import scipy as sp
        >>> try:
        ...    lena = sp.lena()
        ... except AttributeError:
        ...    from scipy import misc
        ...    lena = misc.lena()
    	>>> X = lena.reshape((-1, 1)) # We need an (n_sample, n_feature) array
    	>>> k_means = cluster.KMeans(n_clusters=5, n_init=1)
    	>>> k_means.fit(X) # doctest: +ELLIPSIS
    	KMeans(copy_x=True, init='k-means++', ...
    	>>> values = k_means.cluster_centers_.squeeze()
    	>>> labels = k_means.labels_
    	>>> lena_compressed = np.choose(labels, values)
    	>>> lena_compressed.shape = lena.shape

    .. list-table::
      :class: centered 

      *
        - |lena|

        - |lena_compressed|

        - |lena_regular|

        - |lena_histogram|

      *

        - Raw image

        - K-means quantization

        - Equal bins

        - Image histogram


Hierarchical agglomerative clustering: Ward
---------------------------------------------

A :ref:`hierarchical_clustering` method is a type of cluster analysis
that aims to build a hierarchy of clusters. In general, the various approaches
of this technique are either:

  * **Agglomerative** - bottom-up approaches: each observation starts in its
    own cluster, and clusters are iterativelly merged in such a way to
    minimize a *linkage* criterion. This approach is particularly interesting
    when the clusters of interest are made of only a few observations. When
    the number of clusters is large, it is much more computationally efficient
    than k-means.

  * **Divisive** - top-down approaches: all observations start in one
    cluster, which is iteratively split as one moves down the hierarchy.
    For estimating large numbers of clusters, this approach is both slow (due
    to all observations starting as one cluster, which it splits recursively)
    and statistically ill-posed.

Connectivity-constrained clustering
.....................................

With agglomerative clustering, it is possible to specify which samples can be
clustered together by giving a connectivity graph. Graphs in the scikit
are represented by their adjacency matrix. Often, a sparse matrix is used.
This can be useful, for instance, to retrieve connected regions (sometimes
also referred to as connected components) when
clustering an image:

.. image:: ../../auto_examples/cluster/images/plot_lena_ward_segmentation_001.png
    :target: ../../auto_examples/cluster/plot_lena_ward_segmentation.html
    :scale: 40
    :align: right

.. literalinclude:: ../../auto_examples/cluster/plot_lena_ward_segmentation.py
    :lines: 21-45

..
    >>> from sklearn.feature_extraction.image import grid_to_graph
    >>> connectivity = grid_to_graph(*lena.shape)


Feature agglomeration
......................

We have seen that sparsity could be used to mitigate the curse of
dimensionality, *i.e* an insufficient amount of observations compared to the
number of features. Another approach is to merge together similar
features: **feature agglomeration**. This approach can be implemented by
clustering in the feature direction, in other words clustering the
transposed data.

.. image:: ../../auto_examples/cluster/images/plot_digits_agglomeration_001.png
    :target: ../../auto_examples/cluster/plot_digits_agglomeration.html
    :align: right
    :scale: 57

::

   >>> digits = datasets.load_digits()
   >>> images = digits.images
   >>> X = np.reshape(images, (len(images), -1))
   >>> connectivity = grid_to_graph(*images[0].shape)

   >>> agglo = cluster.FeatureAgglomeration(connectivity=connectivity,
   ...                                      n_clusters=32)
   >>> agglo.fit(X) # doctest: +ELLIPSIS
   FeatureAgglomeration(affinity='euclidean', compute_full_tree='auto',...
   >>> X_reduced = agglo.transform(X)

   >>> X_approx = agglo.inverse_transform(X_reduced)
   >>> images_approx = np.reshape(X_approx, images.shape)

.. topic:: ``transform`` and ``inverse_transform`` methods

   Some estimators expose a ``transform`` method, for instance to reduce
   the dimensionality of the dataset.

Decompositions: from a signal to components and loadings
===========================================================

.. topic:: **Components and loadings**

   If X is our multivariate data, then the problem that we are trying to solve
   is to rewrite it on a different observational basis: we want to learn
   loadings L and a set of components C such that *X = L C*.
   Different criteria exist to choose the components

Principal component analysis: PCA
-----------------------------------

:ref:`PCA` selects the successive components that
explain the maximum variance in the signal.

.. |pca_3d_axis| image:: ../../auto_examples/decomposition/images/plot_pca_3d_001.png
   :target: ../../auto_examples/decomposition/plot_pca_3d.html
   :scale: 70

.. |pca_3d_aligned| image:: ../../auto_examples/decomposition/images/plot_pca_3d_002.png
   :target: ../../auto_examples/decomposition/plot_pca_3d.html
   :scale: 70

.. rst-class:: centered

   |pca_3d_axis| |pca_3d_aligned|

The point cloud spanned by the observations above is very flat in one
direction: one of the three univariate features can almost be exactly
computed using the other two. PCA finds the directions in which the data is
not *flat*

When used to *transform* data, PCA can reduce the dimensionality of the
data by projecting on a principal subspace.

.. np.random.seed(0)

::

    >>> # Create a signal with only 2 useful dimensions
    >>> x1 = np.random.normal(size=100)
    >>> x2 = np.random.normal(size=100)
    >>> x3 = x1 + x2
    >>> X = np.c_[x1, x2, x3]

    >>> from sklearn import decomposition
    >>> pca = decomposition.PCA()
    >>> pca.fit(X)
    PCA(copy=True, n_components=None, whiten=False)
    >>> print(pca.explained_variance_)  # doctest: +SKIP
    [  2.18565811e+00   1.19346747e+00   8.43026679e-32]

    >>> # As we can see, only the 2 first components are useful
    >>> pca.n_components = 2
    >>> X_reduced = pca.fit_transform(X)
    >>> X_reduced.shape
    (100, 2)

.. Eigenfaces here?

Independent Component Analysis: ICA
-------------------------------------

:ref:`ICA` selects components so that the distribution of their loadings carries
a maximum amount of independent information. It is able to recover
**non-Gaussian** independent signals:

.. image:: ../../auto_examples/decomposition/images/plot_ica_blind_source_separation_001.png
   :target: ../../auto_examples/decomposition/plot_ica_blind_source_separation.html
   :scale: 70
   :align: center

.. np.random.seed(0)

::

    >>> # Generate sample data
    >>> time = np.linspace(0, 10, 2000)
    >>> s1 = np.sin(2 * time)  # Signal 1 : sinusoidal signal
    >>> s2 = np.sign(np.sin(3 * time))  # Signal 2 : square signal
    >>> S = np.c_[s1, s2]
    >>> S += 0.2 * np.random.normal(size=S.shape)  # Add noise
    >>> S /= S.std(axis=0)  # Standardize data
    >>> # Mix data
    >>> A = np.array([[1, 1], [0.5, 2]])  # Mixing matrix
    >>> X = np.dot(S, A.T)  # Generate observations

    >>> # Compute ICA
    >>> ica = decomposition.FastICA()
    >>> S_ = ica.fit_transform(X)  # Get the estimated sources
    >>> A_ = ica.mixing_.T
    >>> np.allclose(X,  np.dot(S_, A_) + ica.mean_)
    True

