# Copyright (c) 2011 University of Amsterdam
# Copyright (c) 2009 Leif Johnson <leif@leifjohnson.net>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Adapted for scikit-learn by Lars Buitinck, ILPS, University of Amsterdam.

'''Basic Perceptron library.

Perceptrons are supervised discriminative linear models for classification.
The basic perceptron is simply a vector W in event space that is normal to a
linear separating surface going through the origin. Data points that fall on
one side of the surface are in outcome (+), and data points on the other side
are in outcome (-). To determine the side of the surface on which some event X
falls, we compute the dot product of X and W and return the sign of the
resulting scalar. Learning takes place in the perceptron whenever it makes an
incorrect classification from labeled data : X is simply added to W to
minimally correct the output of the perceptron for this data point. This might
invalidate classifications for other events in the data set, but it is simple.

This particular library generalizes the perceptron to K > 2 outcomes by
maintaining a decision surface for each of the K possible outcomes. Each
surface can be thought of as separating events in one outcome from events in
all other outcomes. At classification time, we compute the kernel value of X
and each of the Ws. The outcome with the greatest resulting scalar value is
returned as the prediction for that event.

This library also incorporates the "kernel trick" to generalize the notion of
the dot product to spaces of different dimensionality. Instead of computing the
dot product between X and W, we use a kernel function that accepts X and W as
inputs and returns a scalar value indicating their similarity in some mapped
comparison space.
'''

from .base import BaseEstimator
from collections import defaultdict
import numpy as np


def RadialBasisKernel(gamma):
    '''This kernel returns the rbf between the input vectors.'''

    def __init__(self, gamma):
        self._minus_gamma = gamma

    def __call__(self, weights, event):
        delta = np.linalg.norm(event - weights)
        return np.exp(self._minus_gamma * delta * delta)


class PolynomialKernel(object):
    '''This kernel returns the dot product raised to some power.'''
    def __init__(self, degree, alpha=1.0):
        self._degree = degree
        self._alpha = alpha

    def __call__(self, weights, event):
        return (np.dot(weights, event) + self._alpha) ** self._degree


class Perceptron(BaseEstimator):
    '''Classifier implementing single-layer perceptron with error-driven
    learning.

    Parameters
    ----------
    kernel : callable, optional
        Kernel function for mapping non-linear data.
    '''

    def __init__(self, kernel=None):
        self._kernel = kernel or np.dot

    def predict(self, X):
        """Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """

        X = np.atleast_2d(X)
        y = np.empty(X.shape[0])

        for i, x in enumerate(X):
            y[i] = self._classify(x)[0]
        return y

    def _classify(self, event):
        '''Classify an event into one of our possible outcomes.

        event: A np vector of the dimensionality of our input space.

        Returns a pair of (outcome, score) for the most likely outcome.
        '''
        return self._max_outcome(self._weights, event)

    def _fit(self, X, y, n_iter, averaged):
        X = np.asanyarray(X)
        y = np.asanyarray(y)

        n_samples, n_features = X.shape
        n_labels = len(np.unique(y))

        self._weights = np.zeros((n_labels, n_features), 'd')
        if averaged:
            self._iterations = np.zeros((n_labels, ), 'i')
            self._survived = np.zeros((n_labels, ), 'i')
            self._history = np.zeros(self._weights.shape, 'd')
            self._acc = np.zeros(self._weights.shape, 'd')

        for i in xrange(n_iter):
            for j in xrange(n_samples):
                self._learn(X[j], y[j])

        return self

    def fit(self, X, y, n_iter=1):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.
        n_iter : int, optional
            Number of iterations to perform. Default 1.

        Returns
        -------
        self
        """

        return self._fit(X, y, n_iter=n_iter, averaged=False)

    def _learn(self, event, outcome):
        '''Adjust the hyperplane based on a classification attempt.

        event: A np vector of the dimensionality of our input space.
        outcome: An integer in [0, outcome_size) indicating the correct outcome
          for this event.
        '''
        pred, _ = self._classify(event)
        if pred != outcome:
            self._update_weights(pred, event, -1)
            self._update_weights(outcome, event, 1)

    def _update_weights(self, class_index, event, delta):
        '''Update the weights for an index based on an event.

        class_index: The index of a weight vector to update.
        event: An event vector to use for the update.
        delta: Weight the update by this value.
        '''
        self._weights[class_index] += delta * event

    def _max_outcome(self, weights, event):
        '''Return the maximum scoring outcome for a set of weights.

        weights: An array of weight vectors, one per outcome.
        event: An event vector.
        '''
        scores = np.array([self._kernel(w, event) for w in weights])
        index = scores.argmax()
        # Return the score as a probability computed locally among our possible
        # outcomes. Subtracting the maximum raw score before exponentiating
        # stabilizes the result numerically, while also encouraging low-scoring
        # values to err by squashing towards 0.
        return index, 1.0 / np.exp(scores - scores[index]).sum()


class AveragedPerceptron(Perceptron):
    '''Classifier based on a weighted sum of perceptrons.

    This perceptron algorithm performs similarly to the basic perceptron when
    learning from labeled data: whenever the predicted outcome for an event
    differs from the true outcome, the weights of the perceptron are updated to
    classify this new event correctly.

    However, in addition to updating the weights of the perceptron in response
    to errors during learning, the AveragedPerceptron also makes a copy of the
    old weight matrix and adds it to a running sum of all past weight matrices.
    The sums are weighted by the number of iterations that each constituent
    weight matrix survived before making an error.

    At classification time, the AveragedPerceptron uses both the current weight
    matrix and the weighted sum of past matrices to make its decision.

    This is equivalent to the "voted perceptron" algorithm described by Freund
    and Schapire (1999). The averaging approach improves on the basic perceptron
    algorithm by providing a "large margin" approach to handling datasets that
    are not linearly separable.
    '''

    def __init__(self, kernel=None):
        super(AveragedPerceptron, self).__init__(kernel)

    def _classify(self, event):
        return self._max_outcome(self._history, event)

    def fit(self, X, y, n_iter=1):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.
        n_iter : int, optional
            Number of iterations to perform. Default 1.

        Returns
        -------
        self
        """

        return self._fit(X, y, n_iter=n_iter, averaged=True)

    def _learn(self, event, outcome):
        self._iterations[outcome] += 1
        pred, score = self._max_outcome(self._weights, event)
        if pred == outcome:
            self._survived[outcome] += 1
        else:
            self._update_history(pred)
            self._update_history(outcome)
            self._update_weights(pred, event, -1)
            self._update_weights(outcome, event, 1)

    def _update_history(self, class_index):
        '''Update the history for a particular class.'''
        s = self._survived[class_index]
        if s > 0:
            acc = self._acc[class_index]
            acc += s * self._weights[class_index]
            self._history[class_index] = acc / self._iterations[class_index]
            self._survived[class_index] = 0


class SparseDot(object):
    '''Kernel for SparseAveragedPerceptron'''
    @staticmethod
    def __call__(weights, event):
        '''Return the dot product of a sparse vector and an event.'''
        return sum(weights.get(f, 0) for f in event)


class SparseAveragedPerceptron(AveragedPerceptron):
    '''A voted perceptron using sparse vectors for storing weights.

    This class achieves sparseness by assuming that all events are
    binary-valued. To use this class with a non-binary event space, you will
    need to manually chop up the event space into binary events. For example,
    suppose event 'foo' can take on N possible discrete values ; you might
    convert 'foo' into N separate events, 'foo=1' .. 'foo=N'.

    In addition, the classifier only retains the top-weighted features for
    each class, making this a sort of beam search version of the general
    Perceptron. This reduces the space requirements for the classifier, at the
    cost of lower accuracy.
    '''

    def __init__(self, beam_width=1000):
        '''Use sparse vectors to store weights and events.'''
        super(SparseAveragedPerceptron, self).__init__()
        self._kernel = SparseDot()
        self._beam_width = beam_width

    def fit(self, X, y, n_iter=1):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.
        n_iter : int, optional
            Number of iterations to perform. Default 1.

        Returns
        -------
        self
        """

        X = np.asanyarray(X)
        y = np.asanyarray(y)

        n_samples, n_features = X.shape
        n_labels = len(np.unique(y))

        self._weights = [defaultdict(float) for i in xrange(n_labels)]
        self._history = [defaultdict(float) for i in xrange(n_labels)]
        self._acc = [defaultdict(float) for i in xrange(n_labels)]
        self._iterations = np.zeros((n_samples, ), 'i')
        self._survived = np.zeros((n_samples, ), 'i')

        for i in xrange(n_iter):
            for j in xrange(n_samples):
                self._learn(X[j], y[j])

        return self

    def _update_history(self, class_index):
        s = self._survived[class_index]
        if s > 0:
            i = self._iterations[class_index]
            a = self._acc[class_index]
            h = self._history[class_index]
            for f, w in self._weights[class_index].iteritems():
                a[f] += s * w
                h[f] = a[f] / i
            self._prune(h)
            self._survived[class_index] = 0

    def _update_weights(self, class_index, event, delta):
        w = self._weights[class_index]
        for f in event:
            w[f] += delta
        self._prune(w)

    def _prune(self, weights):
        '''Prune the weights in a sparse vector to our beam width.'''
        if len(weights) < 1.3 * self._beam_width:
            return
        fws = sorted(weights.iteritems(), key=lambda x: -abs(x[1]))
        for f, _ in fws[self._beam_width:]:
            del weights[f]
