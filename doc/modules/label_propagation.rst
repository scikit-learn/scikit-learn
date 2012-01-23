.. _label_propagation:

===================================================
Label Propagation
===================================================

`sklearn.semisupervised.label_propagation` contains a few variations of semi-supervised
graph inference algorithms. In the semi-supervised classification setting, the
learning algorithm is fed both labeled and unlabeled data. With the addition of
unlabeled data in the training model, the algorithm can better learn the total
structure of the data. These algorithms generally do very well in practice even
when faced with far fewer labeled points than ordinary classification models.

A few features available in this model:
  * Can be used for classification and regression tasks
  * Kernel methods to project data into alternate dimensional spaces

.. currentmodule:: sklearn.label_propagation
.. topic:: Input labels for semi-supervised learning
    It is important to assign an identifier to unlabeled points along with the
    labeled data when training the model with the `fit` method.

This module provides two label propagation models: :class:`LabelPropagation` and
:class:`LabelSpreading`. Both work by forming a fully connected graph for each
item in the input dataset. They differ only in the definition of the matrix
that represents the graph and the clamp effect on the label distributions.
:class:`LabelPropagation` is far more intuitive than :class:`LabelSpreading`
which is motivatived by deeper mathematics.

Clamping
========

Clamping allows the algorithm to change the weight of the true ground labeled
data to some degree. The :class:`LabelPropagation` algorithm performs hard
clamping of input labels, which means :math:`\alpha=1`. This clamping factor
can be relaxed, to say :math:`\alpha=0.8`, which means that we will always
retain 80 percent of our original label distribution, but the algorithm gets to
change it's confidence of the distribution within 20 percent.

Kernel Functions
================

Label propagation models have two built-in kernel methods. Choice of kernel
effects both scalability and performance of the algorithms. The following are
available:

  * rbf (:math:`exp(-\gamma |x-x'|^2), \gamma > 0`). :math:`\gamma` is
    specified by keyword gamma.

  * knn (:math:`1[x' \in kNN(x)]`). :math:`k` is specified by keyword
    n_neighbors.

RBF kernel will produce a fully connected graph which is represented in memory
by a dense matrix. This matrix may be very large and combined with the cost of
performing a full matrix multiplication calculation for each iteration of the
algorithm can lead to prohibitively long running times. On the other hand,
the KNN kernel will produce a much more memory friendly sparse matrix
which can drastically reduce running times.

Examples
========
  * :ref:`example_label_propagation_plot_label_propagation_versus_svm_iris.py`
  * :ref:`example_label_propagation_structure.py`


References
==========
[1] Yoshua Bengio, Olivier Delalleau, Nicolas Le Roux. In Semi-Supervised
Learning (2006), pp. 193-216
