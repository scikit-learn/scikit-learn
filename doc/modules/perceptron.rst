.. Copyright (c) 2011 University of Amsterdam
   Copyright (c) 2009 Leif Johnson <leif@leifjohnson.net>
   
   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   
   Adapted for scikit-learn by Lars Buitinck, ILPS, University of Amsterdam.

===========
Perceptrons
===========

.. currentmodule:: scikits.learn.perceptron

Perceptrons are supervised discriminative linear models for classification.
The basic perceptron is simply a vector w in feature vector space that is
normal to a linear separating surface going through the origin. Data points
that fall on one side of the surface are in class (+), and data points on the
other side are in class (-). To determine the side of the surface on which some
feature vector x falls, we compute the dot product of x and w and return the
sign of the resulting scalar. Learning takes place in the perceptron whenever
it makes an incorrect classification from labeled data (error-driven learning):
x is simply added to w to minimally correct the output of the perceptron for
this data point.

It can be shown that the perceptron eventually converges when run on linearly
separable data. However, for the sake of efficiency and because it never
converges on datasets that are not linearly separable, the perceptrons in this
module only run for an `n_iter` number of iterations, to be set by the user.

This module generalizes the perceptron to K > 2 labels by maintaining a
decision surface for each of the K possible labels, a strategy known as "one
versus all." Each surface can be thought of as separating feature vectors in
one label from feature vectors in all other labels. The label with the greatest
resulting scalar value is returned as the prediction for that feature vector.


Online learning
===============

In addition to the usual `fit`-style batch learning, the `perceptron` module
implements an online learning interface through the methods `partial_setup`
and `partial_fit`. The former must be called to initialize online learning;
the second can then be called as often as required to learn from partial
datasets. This is useful if the dataset is too large to fit in memory at once.

The online learning API can in principle *not* be mixed with the regular batch
learning API. Also note that the averaged perceptron might give slightly
different results when using the online and batch APIs, because the batch API
triggers some optimizations.


Averaged perceptrons
====================

If the `averaged` parameter is true, this classifier maintains a running sum
of past weight matrices, weighted by the number of iterations that each weight
matrix survived before making an error. At classification time, the averaged
perceptron uses both the current weight matrix and the weighted sum of past
matrices to make its decision.

This is an approximation to the "voted perceptron" algorithm described by
Freund and Schapire (1999). The averaging approach improves on the basic
perceptron algorithm by providing a "large margin" approach to handling
datasets that are not linearly separable.

By averaging instead of voting, we gain the benefit that previous weight
vectors are effectively summarized in two `n_labels` * `n_features` matrices,
saving memory and computation time compared to the voted perceptron.


.. topic:: References:

   * Y. Freund and R.E. Schapire (1999),
     `"Large margin classification using the perceptron algorithm"
     <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.48.8200&rep=rep1&type=pdf>`_.
     Machine Learning 37(3):277-296

