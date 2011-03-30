===========
Naive Bayes
===========

.. currentmodule:: scikits.learn.naive_bayes


**Naive Bayes** algorithms are a set of supervised learning methods
based on applying Baye's theorem with strong (naive) independence
assumptions.

The advantage of Naive Bayes approaches are:

   - It requires a small amount of training data to estimate the
     parameters necessary for classification.

   - In spite of their naive design and apparently over-simplified
     assumptions, naive Bayes classifiers have worked quite well in
     many complex real-world situations.

   - The decoupling of the class conditional feature distributions
     means that each distribution can be independently estimated as a
     one dimensional distribution. This in turn helps to alleviate
     problems stemming from the curse of dimensionality.


Gaussian Naive Bayes
--------------------

:class:`GNB` implements the Gaussian Naive Bayes algorithm for classification.



.. topic:: Examples:

 * :ref:`example_naive_bayes.py`,
