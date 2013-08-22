
.. _multiclass:

====================================
Multiclass and multilabel algorithms
====================================

.. currentmodule:: sklearn.multiclass

The :mod:`sklearn.multiclass` module implements *meta-estimators* to perform
``multiclass`` and ``multilabel`` classification. Those meta-estimators are
meant to turn a binary classifier or a regressor into a multi-class/label
classifier.

  - **Multiclass classification** means a classification task with more than
    two classes; e.g., classify a set of images of fruits which may be oranges,
    apples, or pears. Multiclass classification makes the assumption that each
    sample is assigned to one and only one label: a fruit can be either an
    apple or a pear but not both at the same time.

  - **Multilabel classification** assigns to each sample a set of target
    labels. This can be thought as predicting properties of a data-point
    that are not mutually exclusive, such as topics that are relevant for a
    document. A text might be about any of religion, politics, finance or
    education at the same time or none of these.

  - **Multioutput-multiclass classification** and **multi-task classification**
    means that an estimators have to handle
    jointly several classification tasks. This is a generalization
    of the multi-label classification task, where the set of classification
    problem is restricted to binary classification, and of the multi-class
    classification task. *The output format is a 2d numpy array.*

    The set of labels can be different for each output variable.
    For instance a sample could be assigned "pear" for an output variable that
    takes possible values in a finite set of species such as "pear", "apple",
    "orange" and "green" for a second output variable that takes possible values
    in a finite set of colors such as "green", "red", "orange", "yellow"...

    This means that any classifiers handling multi-output
    multiclass or multi-task classification task
    supports the multi-label classification task as a special case.
    Multi-task classification is similar to the multi-output
    classification task with different model formulations. For
    more information, see the relevant estimator documentation.

Estimators in this module are meta-estimators. For example, it is possible to
use these estimators to turn a binary classifier or a regressor into a
multiclass classifier. It is also possible to use these estimators with
multiclass estimators in the hope that their generalization error or runtime
performance improves.

You don't need to use these estimators unless you want to experiment with
different multiclass strategies: all classifiers in scikit-learn support
multiclass classification out-of-the-box. Below is a summary of the
classifiers supported by scikit-learn grouped by strategy:

  - Inherently multiclass: :ref:`Naive Bayes <naive_bayes>`,
    :class:`sklearn.lda.LDA`,
    :ref:`Decision Trees <tree>`, :ref:`Random Forests <forest>`,
    :ref:`Nearest Neighbors <neighbors>`.
  - One-Vs-One: :class:`sklearn.svm.SVC`.
  - One-Vs-All: all linear models except :class:`sklearn.svm.SVC`.

Some estimators also support multioutput-multiclass classification
tasks :ref:`Decision Trees <tree>`, :ref:`Random Forests <forest>`,
:ref:`Nearest Neighbors <neighbors>`.


.. warning::

    For the moment, no metric supports the multioutput-multiclass
    classification task.

Multilabel classification format
================================

In multilabel learning, the joint set of binary classification task
is expressed with either a sequence of sequences or a label binary indicator
array.

In the sequence of sequences format, each set of labels is represented as
a sequence of integer, e.g. ``[0]``, ``[1, 2]``. An empty set of labels is
then expressed as ``[]``, and a set of samples as ``[[0], [1, 2], []]``.
In the label indicator format, each sample is one row of a 2d array of
shape (n_samples, n_classes) with binary values: the one, i.e. the non zero
elements, corresponds to the subset of labels. Our previous example is
therefore expressed as ``np.array([[1, 0, 0], [0, 1, 1], [0, 0, 0])``
and an empty set of labels would be represented by a row of zero elements.


In the preprocessing module, the transformer
:class:`sklearn.preprocessing.label_binarize` and the function
:func:`sklearn.preprocessing.LabelBinarizer`
can help you to convert the sequence of sequences format to the label
indicator format.

  >>> from sklearn.datasets import make_multilabel_classification
  >>> from sklearn.preprocessing import LabelBinarizer
  >>> X, Y = make_multilabel_classification(n_samples=5, random_state=0)
  >>> Y
  ([0, 1, 2], [4, 1, 0, 2], [4, 0, 1], [1, 0], [3, 2])
  >>> LabelBinarizer().fit_transform(Y)
  array([[1, 1, 1, 0, 0],
         [1, 1, 1, 0, 1],
         [1, 1, 0, 0, 1],
         [1, 1, 0, 0, 0],
         [0, 0, 1, 1, 0]])

.. warning::

    - The sequence of sequences format will disappear in a near future.
    - All estimators or functions support both multilabel format.


One-Vs-The-Rest
===============

This strategy, also known as **one-vs-all**, is implemented in
:class:`OneVsRestClassifier`.  The strategy consists in fitting one classifier
per class. For each classifier, the class is fitted against all the other
classes. In addition to its computational efficiency (only `n_classes`
classifiers are needed), one advantage of this approach is its
interpretability. Since each class is represented by one and one classifier
only, it is possible to gain knowledge about the class by inspecting its
corresponding classifier. This is the most commonly used strategy and is a fair
default choice.

Multiclass learning
-------------------

Below is an example of multiclass learning using OvR::

  >>> from sklearn import datasets
  >>> from sklearn.multiclass import OneVsRestClassifier
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> OneVsRestClassifier(LinearSVC(random_state=0)).fit(X, y).predict(X)
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

Multilabel learning
-------------------

:class:`OneVsRestClassifier` also supports multilabel classification.
To use this feature, feed the classifier a list of tuples containing
target labels, like in the example below.


.. figure:: ../auto_examples/images/plot_multilabel_1.png
    :target: ../auto_examples/plot_multilabel.html
    :align: center
    :scale: 75%


.. topic:: Examples:

    * :ref:`example_plot_multilabel.py`


One-Vs-One
==========

:class:`OneVsOneClassifier` constructs one classifier per pair of classes.
At prediction time, the class which received the most votes is selected.
Since it requires to fit `n_classes * (n_classes - 1) / 2` classifiers,
this method is usually slower than one-vs-the-rest, due to its
O(n_classes^2) complexity. However, this method may be advantageous for
algorithms such as kernel algorithms which don't scale well with
`n_samples`. This is because each individual learning problem only involves
a small subset of the data whereas, with one-vs-the-rest, the complete
dataset is used `n_classes` times.

Multiclass learning
-------------------

Below is an example of multiclass learning using OvO::

  >>> from sklearn import datasets
  >>> from sklearn.multiclass import OneVsOneClassifier
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> OneVsOneClassifier(LinearSVC(random_state=0)).fit(X, y).predict(X)
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])


Error-Correcting Output-Codes
=============================

Output-code based strategies are fairly different from one-vs-the-rest and
one-vs-one. With these strategies, each class is represented in a euclidean
space, where each dimension can only be 0 or 1. Another way to put it is
that each class is represented by a binary code (an array of 0 and 1). The
matrix which keeps track of the location/code of each class is called the
code book. The code size is the dimensionality of the aforementioned space.
Intuitively, each class should be represented by a code as unique as
possible and a good code book should be designed to optimize classification
accuracy. In this implementation, we simply use a randomly-generated code
book as advocated in [2]_ although more elaborate methods may be added in the
future.

At fitting time, one binary classifier per bit in the code book is fitted.
At prediction time, the classifiers are used to project new points in the
class space and the class closest to the points is chosen.

In :class:`OutputCodeClassifier`, the `code_size` attribute allows the user to
control the number of classifiers which will be used. It is a percentage of the
total number of classes.

A number between 0 and 1 will require fewer classifiers than
one-vs-the-rest. In theory, ``log2(n_classes) / n_classes`` is sufficient to
represent each class unambiguously. However, in practice, it may not lead to
good accuracy since ``log2(n_classes)`` is much smaller than n_classes.

A number greater than than 1 will require more classifiers than
one-vs-the-rest. In this case, some classifiers will in theory correct for
the mistakes made by other classifiers, hence the name "error-correcting".
In practice, however, this may not happen as classifier mistakes will
typically be correlated. The error-correcting output codes have a similar
effect to bagging.


Multiclass learning
-------------------

Below is an example of multiclass learning using Output-Codes::

  >>> from sklearn import datasets
  >>> from sklearn.multiclass import OutputCodeClassifier
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> clf = OutputCodeClassifier(LinearSVC(random_state=0),
  ...                            code_size=2, random_state=0)
  >>> clf.fit(X, y).predict(X)
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1,
         1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

.. topic:: References:

    .. [1] "Solving multiclass learning problems via error-correcting output codes",
        Dietterich T., Bakiri G.,
        Journal of Artificial Intelligence Research 2,
        1995.

    .. [2] "The error coding method and PICTs",
        James G., Hastie T.,
        Journal of Computational and Graphical statistics 7,
        1998.

    .. [3] "The Elements of Statistical Learning",
        Hastie T., Tibshirani R., Friedman J., page 606 (second-edition)
        2008.
