.. _model_evaluation:

===================
Model evaluation
===================

.. TODO

  Metrics
  =======


Random classification accuracy
=================================

.. currentmodule:: sklearn.random_classifier

When doing classification, a simple sanity check consists in comparing
one's classifier against the accuracy of a purely random classifier.
:class:`RandomClassifier` implements three different strategies for
random classification. `stratified` generates predictions by respecting
the training set's class distribution. `most_frequent` always predicts the
most frequent label in the training set (useful for binary classification).
`uniform` generates predictions uniformly at random. For example, let's compare
the accuracy of `SVC` and `most_frequent`::

  >>> from sklearn.datasets import load_iris
  >>> from sklearn.random_classifier import RandomClassifier
  >>> from sklearn.svm import SVC
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> y[y != 1] = -1
  >>> clf = SVC(kernel='linear', C=1).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.73...
  >>> clf = RandomClassifier(sampling='most_frequent', random_state=0).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.66...

When the accuracy of a classifier is too close to random classification, it
probably means that something went wrong: features are not helpful, a
hyparameter is not correctly tuned, the classifier is suffering from class
imbalance, etc... In the above example, changing the kernel to `rbf`,
the accuracy goes from 0.73 to almost 1.0.
