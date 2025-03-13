.. currentmodule:: sklearn

.. _manifold:

=================
Manifold Learning
=================

.. rst-class:: quote

                 | Look for the bare necessities
                 | The simple bare necessities
                 | Forget about your worries and your strife
                 | I mean the bare necessities
                 | Old Mother Nature's recipes
                 | That bring the bare necessities of life
                 |    -- Baloo's song [The Jungle Book]


.. topic:: Note

   For the time being, scikit-learn only provides implementations of a small
   subset of the existing manifold learning algorithms.

.. topic:: Overview

    This module implements data embedding techniques that are often useful for
    data visualization or as preprocessing steps. The user guide includes
    information regarding the optimal selection of these methods for various
    problems, as well as the limitations inherent to each individual algorithm.

    The following embedding techniques are provided:

    Combinatorial algorithms:

    * Locally Linear Embedding performs a non-linear low-dimensional
      embedding, preserving distances within neighborhoods. It can be thought
      of as a series of local Principal Component Analyses which are globally
      compared to find the best non-linear embedding. :class:`locally_linear.LocallyLinearEmbedding`

    * Isomap seeks a lower-dimensional embedding which maintains geodesic
      distances between all points. :class:`manifold.Isomap`

    * t-distributed Stochastic Neighbor Embedding (t-SNE) converts
      affinities of data points to probabilities. The affinities in the original
      space are represented by Gaussian joint probabilities and the affinities
      in the embedded space are represented by Student's t-distributions. This
      allows t-SNE to be particularly sensitive to local structure and has a few
      other advantages over existing techniques:

         * Revealing the structure at many scales on a single map
         * Revealing data that lie in multiple, different, manifolds or clusters
         * Reducing the tendency to crowd points together at the center

      t-SNE has two parameters that have pronounced effects on the resulting
      maps: perplexity and early exaggeration factor.
      :class:`manifold.TSNE`

    * t-distributed Stochastic Neighbor Embedding with Particle Swarm Optimization (t-SNE-PSO)
      uses the same probabilistic formulation as t-SNE but replaces gradient descent with
      Particle Swarm Optimization for the embedding step. Introduced by Allaoui et al. (2025),
      this approach helps avoid local minima and often produces better cluster separation in the
      resulting embedding by using a dynamic update of cognitive and social coefficients to balance
      exploration and exploitation in the optimization landscape.
      :class:`manifold.TSNEPSO`

Spectral techniques:

    * Spectral Embedding, including Laplacian Eigenmaps, a simple approach
      to low-dimensional representation. It computes a low-dimensional
      representation of the data using a spectral decomposition of the graph
      Laplacian. The graph generated can be considered as a discrete
      approximation of the low dimensional manifold in the high dimensional
      space. :class:`manifold.SpectralEmbedding`

    * Multi-dimensional Scaling projects the data in a low-dimensional
      space where the distances respect well the distances in the original
      high-dimensional space. :class:`manifold.MDS`

Non-linear dimensionality reduction through Kernel PCA:

    * Kernel PCA is an extension of PCA achieved by the use of kernel methods.
      :class:`decomposition.KernelPCA`

Dimensionality reduction for neighborhood graphs:

    * Neighborhood Components Analysis is a probabilistic approach to
      dimensionality reduction. The algorithm finds a projection of the
      input space to a lower-dimensional space that preserves distances
      between neighbors with high probability. The method is non-parametric,
      which means that it does not learn an explicit mapping function from
      the training data. :class:`neighbors.NeighborhoodComponentsAnalysis`

Visualizing the stock market structure:

    For a more detailed description of this example, see the
    :ref:`stock market structure example <stock_market>`.

    For an example of dimensionality reduction on text data see
    :ref:`sphx_glr_auto_examples_text_plot_document_clustering.py`.

Manifold learning on handwritten digits: Locally Linear Embedding, Isomap...:

    An illustration of various embeddings on the digits dataset.

    Given a 8x8 image of a digit, we want to recover a "semantic" 2D
    representation of this digit. Each digit (0-9) is a point in the 2D plane,
    and the embedding should separate the different classes.

    See :ref:`sphx_glr_auto_examples_manifold_plot_lle_digits.py` for the
    example.

.. topic:: References

    .. [1] Roweis, S. & Saul, L. Nonlinear dimensionality reduction by
        locally linear embedding.  Science 290:2323 (2000).

    .. [2] Tenenbaum, J.B.; de Silva, V. & Langford, J.C. A global geometric
        framework for nonlinear dimensionality reduction.  Science 290:2319 (2000).

    .. [3] Belkin, M. & Niyogi, P. Laplacian Eigenmaps for Dimensionality
        Reduction and Data Representation.
        Neural Computation, 15 (2003).

    .. [4] Zhang, Z. & Wang, J. MLLE: Modified Locally Linear Embedding
        Using Multiple Weights.
        `<https://citeseerx.ist.psu.edu/doc_view/pid/b0453cba43b1b4823fda04e56673127b9d6265af>`_

    .. [5] Zhang, Z. & Zha, H. Principal manifolds and nonlinear
        dimensionality reduction via tangent space alignment.
        `<https://archive.siam.org/meetings/sdm04/proceedings/sdm04_025.pdf>`_

    .. [6] Van Der Maaten, L.J.P.; Postma, E.O. & van den Herik, H.J. (2009) Dimensionality
        Reduction: A Comparative Review. Journal of Machine Learning Research, 10 (1-41), 66-71.
        `<https://lvdmaaten.github.io/publications/papers/TR_Dimensionality_Reduction_Review_2009.pdf>`_

    .. [7] Kennedy, J. and Eberhart, R., 1995. Particle swarm optimization.
        In Proceedings of ICNN'95 - International Conference on Neural Networks,
        Vol. 4, pp. 1942-1948.

    .. [8] Shi, Y. and Eberhart, R., 1998. A modified particle swarm optimizer.
        In 1998 IEEE International Conference on Evolutionary Computation Proceedings,
        pp. 69-73.
        
    .. [9] Allaoui, M., Belhaouari, S.B., Hedjam, R., Bouanane, K., & Kherfi, M.L. (2025).
        t-SNE-PSO: A PSO-based optimization approach for t-distributed Stochastic Neighbor 
        Embedding. Expert Systems with Applications. 