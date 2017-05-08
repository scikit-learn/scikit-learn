.. currentmodule:: sklearn.dummy

.. _dummy_estimators:

Dummy estimators
================


When doing supervised learning, a simple sanity check consists of comparing
one's estimator against simple rules of thumb. :class:`DummyClassifier`
implements several such simple strategies for classification:

- ``stratified`` generates random predictions by respecting the training
  set class distribution.
- ``most_frequent`` always predicts the most frequent label in the training set.
- ``prior`` always predicts the class that maximizes the class prior
  (like ``most_frequent`) and ``predict_proba`` returns the class prior.
- ``uniform`` generates predictions uniformly at random.
- ``constant`` always predicts a constant label that is provided by the user.
   A major motivation of this method is F1-scoring, when the positive class
   is in the minority.

Note that with all these strategies, the ``predict`` method completely ignores
the input data!

To illustrate :class:`DummyClassifier`, first let's create an imbalanced
dataset::

  >>> from sklearn.datasets import load_iris
  >>> from sklearn.model_selection import train_test_split
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> y[y != 1] = -1
  >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

Next, let's compare the accuracy of ``SVC`` and ``most_frequent``::

  >>> from sklearn.dummy import DummyClassifier
  >>> from sklearn.svm import SVC
  >>> clf = SVC(kernel='linear', C=1).fit(X_train, y_train)
  >>> clf.score(X_test, y_test) # doctest: +ELLIPSIS
  0.63...
  >>> clf = DummyClassifier(strategy='most_frequent',random_state=0)
  >>> clf.fit(X_train, y_train)
  DummyClassifier(constant=None, random_state=0, strategy='most_frequent')
  >>> clf.score(X_test, y_test)  # doctest: +ELLIPSIS
  0.57...

We see that ``SVC`` doesn't do much better than a dummy classifier. Now, let's
change the kernel::

  >>> clf = SVC(kernel='rbf', C=1).fit(X_train, y_train)
  >>> clf.score(X_test, y_test)  # doctest: +ELLIPSIS
  0.97...

We see that the accuracy was boosted to almost 100%.  A cross validation
strategy is recommended for a better estimate of the accuracy, if it
is not too CPU costly. For more information see the :ref:`cross_validation`
section. Moreover if you want to optimize over the parameter space, it is highly
recommended to use an appropriate methodology; see the :ref:`grid_search`
section for details.

More generally, when the accuracy of a classifier is too close to random, it
probably means that something went wrong: features are not helpful, a
hyperparameter is not correctly tuned, the classifier is suffering from class
imbalance, etc...

:class:`DummyRegressor` also implements four simple rules of thumb for regression:

- ``mean`` always predicts the mean of the training targets.
- ``median`` always predicts the median of the training targets.
- ``quantile`` always predicts a user provided quantile of the training targets.
- ``constant`` always predicts a constant value that is provided by the user.

In all these strategies, the ``predict`` method completely ignores
the input data.
