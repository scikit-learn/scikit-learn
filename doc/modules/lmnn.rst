.. _lmnn:

=============================
Large Margin Nearest Neighbor
=============================

.. sectionauthor:: John Chiotellis <johnyc.code@gmail.com>

.. currentmodule:: sklearn.neighbors

Large Margin Nearest Neighbor (LMNN, :class:`LargeMarginNearestNeighbor`) is a distance metric learning algorithm which aims to improve the accuracy of nearest neighbors classification compared to the standard Euclidean distance.

.. |illustration_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_illustration_001.png
   :target: ../auto_examples/neighbors/plot_lmnn_illustration.html
   :scale: 50

.. |illustration_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_illustration_002.png
   :target: ../auto_examples/neighbors/plot_lmnn_illustration.html
   :scale: 50

.. centered:: |illustration_1| |illustration_2|


For each training sample, the algorithm fixes :math:`k` "target neighbors", namely the :math:`k`-nearest training samples (as measured by the Euclidean distance) that share the same label. Given these target neighbors, LMNN learns a linear transformation of the data by optimizing a trade-off between two goals. The first one is to make each (transformed) point closer to its target neighbors than to any differently-labeled point by a large margin, thereby enclosing the target neighbors in a sphere around the reference sample. Data samples from different classes that violate this margin are called "impostors". The second goal is to minimize the distances of each sample to its target neighbors, which can be seen as some sort of regularization.

Combined with a nearest neighbors classifier (:class:`KNeighborsClassifier`), this method is attractive for classification because it can naturally handle multi-class problems without any increase in the model size, and only a single parameter (``n_neighbors``) has to be selected by the user before training.

Large Margin Nearest Neighbor classification has been shown to work well in practice for data sets of varying size and difficulty. In contrast to related methods such as LDA, LMNN does not make any assumptions about the class distributions. The nearest neighbor classification can naturally produce highly irregular decision boundaries.

.. |classification_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_classification_001.png
   :target: ../auto_examples/neighbors/plot_lmnn_classification.html
   :scale: 50

.. |classification_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_classification_002.png
   :target: ../auto_examples/neighbors/plot_lmnn_classification.html
   :scale: 50

.. centered:: |classification_1| |classification_2|


The plot shows decision boundaries for nearest neighbor classification and large margin nearest neighbor classification.


Dimensionality reduction
========================

:class:`LargeMarginNearestNeighbor` can be used to perform supervised dimensionality reduction. The input data are projected onto a linear subspace consisting of the directions which minimize the LMNN objective. The desired dimensionality can be set using the parameter ``n_features_out``.
For instance, the following shows a comparison of dimensionality reduction with :class:`PCA`, :class:`LinearDiscriminantAnalysis` and :class:`LargeMarginNearestNeighbor` on the Olivetti dataset, a dataset with size :math:`n_{samples} = 400` and :math:`n_{features} = 64 \times 64 = 4096`. The data set is splitted in a training and test set of equal size. For evaluation the 3-nearest neighbor classification accuracy is computed on the 2-dimensional embedding found by each method. Each data sample belongs to one of 40 classes.

.. |dim_reduction_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_dim_reduction_001.png
   :target: ../auto_examples/neighbors/plot_lmnn_dim_reduction.html
   :scale: 33%

.. |dim_reduction_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_dim_reduction_002.png
   :target: ../auto_examples/neighbors/plot_lmnn_dim_reduction.html
   :scale: 33%

.. |dim_reduction_3| image:: ../auto_examples/neighbors/images/sphx_glr_plot_lmnn_dim_reduction_003.png
   :target: ../auto_examples/neighbors/plot_lmnn_dim_reduction.html
   :scale: 33%

.. centered:: |dim_reduction_1| |dim_reduction_2| |dim_reduction_3|


Mathematical formulation
========================

LMNN learns a linear transformation matrix :math:`L` of size ``(n_features_out, n_features)``.
The objective function consists of two competing terms, the pull loss that pulls target neighbors closer to their reference sample and the push loss that pushes impostors away:

.. math::

    \varepsilon_{\text{pull}} (L) = \sum_{i, j \rightsquigarrow i} ||L(x_i - x_j)||^2,

.. math::

    \varepsilon_{\text{push}} (L) = \sum_{i, j \rightsquigarrow i}
    \sum_{l} (1 - y_{il}) [1 + || L(x_i - x_j)||^2 - || L
    (x_i - x_l)||^2]_+,

where :math:`y_{il} = 1` if :math:`y_i = y_l` and :math:`0` otherwise, :math:`[x]_+ = \max(0, x)` is the hinge loss, and :math:`j \rightsquigarrow i` means that the :math:`j^{th}` sample is a target neighbor of the :math:`i^{th}` sample.

LMNN solves the following (nonconvex) minimization problem:

.. math::

    \min_L \varepsilon(L) = (1 - \mu) \varepsilon_{\text{pull}} (L) +
    \mu \varepsilon_{\text{push}} (L) \text{, } \quad \mu \in [0,1].

The parameter :math:`\mu` is a trade-off between penalizing large distances to target neighbors and penalizing margin violations by impostors. In practice, the two terms are weighted equally.

To use this model for classification, one can simply fit a nearest neighbors classifier (:class:`KNeighborsClassifier`) on the data transformed by the linear transformation :math:`L` found by LMNN (this is implemented by :func:`LargeMarginNearestNeighbor.transform`).


Mahalanobis distance
--------------------

LMNN can be seen as learning a (squared) Mahalanobis distance metric:

.. math::

    || L(x_i - x_j)||^2 = (x_i - x_j)^TM(x_i - x_j),

where :math:`M = L^T L` is a symmetric positive semi-definite matrix of size ``(n_features, n_features)``. The objective function of LMNN can be rewritten and solved with respect to :math:`M` directly. This results in a convex but constrained problem (since :math:`M` must be symmetric positive semi-definite). Please, see references for more details.


Implementation
==============

This implementation follows closely the MATLAB implementation found at https://bitbucket.org/mlcircus/lmnn which solves the unconstrained problem. It finds a linear transformation :math:`L` by optimization with L-BFGS instead of solving the constrained problem that finds the globally optimal distance metric. Different from the paper, the problem solved by this implementation is with the *squared* hinge loss (to make the problem differentiable). The parameter :math:`\mu` is fixed to :math:`0.5`.

See the examples below and the doc string of :meth:`LargeMarginNearestNeighbor.fit` for further information.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_neighbors_plot_lmnn_classification.py` 
 * :ref:`sphx_glr_auto_examples_neighbors_plot_lmnn_dim_reduction.py`


.. topic:: References:

   * | `"Distance Metric Learning for Large Margin Nearest Neighbor Classification"
       <http://jmlr.csail.mit.edu/papers/volume10/weinberger09a/weinberger09a.pdf>`_,
     | Weinberger, Kilian Q., and Lawrence K. Saul, Journal of Machine Learning Research, Vol. 10, Feb. 2009, pp. 207-244.

   * `Wikipedia entry on Large Margin Nearest Neighbor
     <https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor>`_

