Tutorial
========

Load a dataset
--------------

>>> from scikits.learn import datasets
>>> iris = datasets.load('iris')

A dataset is a dictionary-like object that holds all the samples and
some metadata about the samples. You can access the underlying data
with members .data and .target.


Prediction
----------
Suppose some given data points each belong to one of two classes, and
the goal is to decide which class a new data point will be in. In
``scikits.learn`` this is done with an *estimator*. An *estimator* is
just a plain Python class that implements the methods fit(X, Y) and
predict(T).

An example of predictor is the class
``scikits.learn.neighbors.Neighbors``(XXX ref). The constructor of a predictor
takes as arguments the parameters of the model. In this case, our only
parameter is k, the number of neighbors to consider.

>>> from scikits.learn import neighbors
>>> clf = neighbors.Neighbors(k=3)

The predictor now must be fitted to the model, that is, it must
`learn` from the model. This is done by passing our training set to
the ``fit`` method.

>>> clf.fit(iris.data, iris.target) #doctest: +ELLIPSIS
<scikits.learn.neighbors.Neighbors instance at 0x...>

Now you can predict new values

>>> print clf.predict([[0, 0, 0, 0]])
[[ 0.]]


Regression
----------
In the regression problem, classes take continous values.
Linear Regression. TODO
