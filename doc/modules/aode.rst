.. _aode:

===========
Averaged One-Dependance Estimator
===========

.. currentmodule:: sklearn.aode


AODE achieves highly accurate classification by averaging over 
all of a small space of alternative naive-Bayes-like models 
that have weaker (and hence less detrimental) independence assumptions 
than naive Bayes. AODE only works on categorical data. Given a 
class variable :math:`y` and a dependent feature vector :math:`x_1` 
through :math:`x_n`, Bayes' theorem states the following relationship:

.. math::

   P(y \mid x_1, \dots, x_n) = \frac{P(x_1, \dots x_n, y)}
                                    {P(x_1, \dots, x_n)}

Since :math:`P(x_1, \dots, x_n)` is constant given the input,
we can use the following classification rule:

.. math::

   P(y | x_1, \dots, x_n) \propto P(y, x_1, \dots, x_n)

For any :math:`1 \leq i \leq n`

.. math::

   P(y, x_1, \dots, x_n) = P(y,x_i)\,P(x_1, \dots, x_n|y,x_i),

Under the assumption that :math:`x_1, x_2, \dots, x_n` are independent 
given :math:`y` and :math:`x_i`, it follows that:

.. math::

   P(y, x_1, \dots, x_n) = P(y,x_i) \prod_{j=1}^{n} P(x_j|y,x_i)

AODE estimates :math:`P(y | x_1, \dots, x_n)` by averaging this
for all values of :math:`i`

.. math::

    P(y|x_1, \dots, x_n)=\frac{\sum_{i=1}^{n}P(y,x_i) \prod_{j=1}^{n} P(x_j|y,x_i)}{\sum_{y' \in Y}\sum_{i=1}^{n}P(y',x_i) \prod_{j=1}^{n} P(x_j|y',x_i)} 

In order to avoid inaccurate estimates, we exclude variables where there 
are fewer than :math:`m` samples of with the same value of :math:`x_i` in the training data. 
This :math:`m` is set to 1 by default.

.. math::

    P(y|x_1, \dots, x_n)=\frac{\sum_{i:1<i<n∧F(x_i)≥m}P(y,x_i) \prod_{j=1}^{n} P(x_j|y,x_i)}{\sum_{y' \in Y}\sum_{i:1<i<n∧F(x_i)≥m}P(y',x_i) \prod_{j=1}^{n} P(x_j|y',x_i)} 

Where :math:`F(x_i)` is a count of the number of training examples having attribute-value :math:`x_i`
If every variable is excluded, AODE defaults to naive Bayes.

The smoothing parameter accounts for features not present in the learning 
samples and prevents zero probabilities in further computations. Setting 
is called Laplace smoothing, while is called Lidstone smoothing.

AODE has a time complexity of :math:`O(tn^2)` for training, and :math:`O(kn^2)`, per example.
:math:`t` is the number of training examples, :math:`n` is the number of attributes, and :math:`k` in the number 
of classes. The space complexity of AODE is :math:`O(k(nv)^2)`, where :math:`v` is 
the average number of distinct values for an attribute.

.. topic:: References:

 * G.I. Webb (2005). `Not so naive bayes: Aggregating one-dependence estimators.
   <https://link.springer.com/article/10.1007/s10994-005-4258-6>`_
   Proc. FLAIRS.

