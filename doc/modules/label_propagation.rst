.. _semi_supervised:

===================================================
Semi-Supervised
===================================================

Label Propagation
=================

`sklearn.semi_supervised.label_propagation` contains a few variations of semi-supervised
graph inference algorithms. In the semi-supervised classification setting, the
learning algorithm is fed both labeled and unlabeled data. The algorithm can better
learn the total structure of the data by knowing how the unlabeled data is distributed.
These algorithms can perform well when we have a very small amount of labeled points
and a large amount of unlabeled points.

A few features available in this model:
  * Can be used for classification and regression tasks
  * Kernel methods to project data into alternate dimensional spaces

.. currentmodule:: sklearn.semi_supervised
.. topic:: Input labels

    It is important to assign an identifier to unlabeled points along with the
    labeled data when training the model with the `fit` method. The identifier
    that this implementation uses the integer value :math:`-1`.

This module provides two label propagation models: :class:`LabelPropagation` and
:class:`LabelSpreading`. Both work by constructing a similarity graph over all
items in the input dataset. They differ in modifications to the similarity
matrix that graph and the clamping effect on the label distributions.

Clamping allows the algorithm to change the weight of the true ground labeled
data to some degree. The :class:`LabelPropagation` algorithm performs hard
clamping of input labels, which means :math:`\alpha=1`. This clamping factor
can be relaxed, to say :math:`\alpha=0.8`, which means that we will always
retain 80 percent of our original label distribution, but the algorithm gets to
change it's confidence of the distribution within 20 percent.

:class:`LabelPropagation` uses the raw similarity matrix constructed from the
data with no modifications. In contrast, :class:`LabelSpreading` minimizes a
loss function that has regularization properties. The algorithm iterates on
a modified version of the original graph and normalizes the edge weights by
computing the normalized graph Laplacian matrix. This procedure is also used in
:class:`sklearn.cluster.SpectralClustering`.

Label propagation models have two built-in kernel methods. Choice of kernel
effects both scalability and performance of the algorithms. The following are
available:

  * rbf (:math:`\exp(-\gamma |x-y|^2), \gamma > 0`). :math:`\gamma` is
    specified by keyword gamma.

  * knn (:math:`1[x' \in kNN(x)]`). :math:`k` is specified by keyword
    n_neighbors.

RBF kernel will produce a fully connected graph which is represented in memory
by a dense matrix. This matrix may be very large and combined with the cost of
performing a full matrix multiplication calculation for each iteration of the
algorithm can lead to prohibitively long running times. On the other hand,
the KNN kernel will produce a much more memory friendly sparse matrix
which can drastically reduce running times.

.. topic:: Examples

  * :ref:`example_semi_supervised_plot_label_propagation_versus_svm_iris.py`
  * :ref:`example_semi_supervised_plot_label_propagation_structure.py`


.. topic:: References

    [1] Yoshua Bengio, Olivier Delalleau, Nicolas Le Roux. In Semi-Supervised
    Learning (2006), pp. 193-216

    [2] Olivier Delalleau, Yoshua Bengio, Nicolas Le Roux. Efficient
    Non-Parametric Function Induction in Semi-Supervised Learning. AISTAT 2005
