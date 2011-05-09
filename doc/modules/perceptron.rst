===========
Perceptrons
===========

.. currentmodule:: scikits.learn.perceptron

Perceptrons are supervised discriminative linear models for classification.
The basic perceptron is simply a vector W in event space that is normal to a
linear separating surface going through the origin. Data points that fall on
one side of the surface are in class (+), and data points on the other side
are in class (-). To determine the side of the surface on which some event X
falls, we compute the dot product of X and W and return the sign of the
resulting scalar. Learning takes place in the perceptron whenever it makes an
incorrect classification from labeled data : X is simply added to W to
minimally correct the output of the perceptron for this data point. This might
invalidate classifications for other events in the data set, but it is simple.

This module generalizes the perceptron to K > 2 labels by maintaining a
decision surface for each of the K possible labels. Each surface can be
thought of as separating events in one label from events in all other labels.
At classification time, we compute the kernel value of X and each of the Ws.
The label with the greatest resulting scalar value is returned as the
prediction for that event.

This module also incorporates the "kernel trick" to generalize the notion of
the dot product to spaces of different dimensionality. Instead of computing the
dot product between X and W, we use a kernel function that accepts X and W as
inputs and returns a scalar value indicating their similarity in some mapped
comparison space.


Online learning
===============
In addition to the usual `fit`-style batch learning, the `perceptron` module
implements an online learning interface through the methods `partial_setup`
and `partial_fit`. The former must be called to initialize online learning;
the second can then be called as often as required to learn from partial
datasets. This is useful if the dataset is too large to fit in memory at once.


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


Sparse perceptrons
==================
The :class:`SparsePerceptron` stores its weight vectors in sparse arrays.
It assumes that all events are binary-valued. To use this class with a
non-binary event space, you will need to chop up the event space into binary
events.

In addition, this class may be configured to retain only the top-weighted
features for each label, making this a sort of beam search version of the
general Perceptron. This reduces the space requirements for the classifier,
at the cost of lower accuracy.


.. topic:: References:

   * Y. Freund and R.E. Schapire (1999),
     `"Large margin classification using the perceptron algorithm"
     <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.48.8200&rep=rep1&type=pdf>`_,
     Machine Learning 37(3):277-296

