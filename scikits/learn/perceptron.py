# Copyright (c) 2009 Leif Johnson <leif@leifjohnson.net>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''Basic Perceptron library.

Perceptrons are discriminative models for classifying data points (also called
"events") into one of a number of discrete categories (also called "outcomes").
The basic Perceptron is simply a vector W in event space that is normal to a
linear separating surface going through the origin. Data points that fall on one
side of the surface are in outcome (+), and data points on the other side are in
outcome (-). To determine the side of the surface on which some event X falls,
we compute the dot product of X and W and return the sign of the resulting
scalar. Learning takes place in the Perceptron whenever it makes an incorrect
classification from labeled data : X is simply added to W to minimally correct
the output of the Perceptron for this data point. This might invalidate
classifications for other events in the data set, but it is simple.

This particular library generalizes the Perceptron to K > 2 outcomes by
maintaining a decision surface for each of the K possible outcomes. Each surface
can be thought of as separating events in one outcome from events in all other
outcomes. At classification time, we compute the kernel value of X and each of
the Ws. The outcome with the greatest resulting scalar value is returned as the
prediction for that event.

This library also incorporates the "kernel trick" to generalize the notion of
the dot product to spaces of different dimensionality. Instead of computing the
dot product between X and W, we use a Kernel function that accepts X and W as
inputs and returns a scalar value indicating their similarity in some mapped
comparison space.
'''

import numpy
from collections import defaultdict


def radial_basis_kernel(gamma):
    '''This kernel returns the rbf between the input vectors.'''
    def run(weights, event):
        delta = numpy.linalg.norm(event - weights)
        return numpy.exp(-gamma * delta * delta)
    return run


def polynomial_kernel(degree, alpha=1.0):
    '''This kernel returns the dot product raised to some power.'''
    def run(weights, event):
        return (numpy.dot(weights, event) + alpha) ** degree
    return run


class Perceptron(object):
    '''A Perceptron is a discriminative machine learning algorithm.'''

    def __init__(self, event_size, outcome_size, kernel=None):
        '''Initialize this Perceptron.

        event_size: An integer indicating the dimensionality of the event space.
        outcome_size: An integer indicating the number of classes in the data.
        kernel: A Kernel object that we can use to compute the distance between
          a data point and our weight vector.
        '''
        assert outcome_size >= 2
        assert event_size >= 1
        self._weights = numpy.zeros((outcome_size, event_size), 'd')
        self._kernel = kernel or numpy.dot

    def classify(self, event):
        '''Classify an event into one of our possible outcomes.

        event: A numpy vector of the dimensionality of our input space.

        Returns a pair of (outcome, score) for the most likely outcome.
        '''
        return Perceptron._max_outcome(self._weights, event, self._kernel)

    def learn(self, event, outcome):
        '''Adjust the hyperplane based on a classification attempt.

        event: A numpy vector of the dimensionality of our input space.
        outcome: An integer in [0, outcome_size) indicating the correct outcome
          for this event.
        '''
        pred, _ = self.classify(event)
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

    @staticmethod
    def _max_outcome(weights, event, kernel):
        '''Return the maximum scoring outcome for a set of weights.

        weights: An array of weight vectors, one per outcome.
        event: An event vector.
        kernel: A kernel function to combine weights with the event.
        '''
        scores = numpy.array([kernel(w, event) for w in weights])
        index = scores.argmax()
        # Return the score as a probability computed locally among our possible
        # outcomes. Subtracting the maximum raw score before exponentiating
        # stabilizes the result numerically, while also encouraging low-scoring
        # values to err by squashing towards 0.
        return index, 1.0 / numpy.exp(scores - scores[index]).sum()


class AveragedPerceptron(Perceptron):
    '''A weighted sum of individual Perceptrons.

    This Perceptron algorithm performs similarly to the basic Perceptron when
    learning from labeled data : Whenever the predicted outcome for an event
    differs from the true outcome, the weights of the Perceptron are updated to
    classify this new event correctly.

    However, in addition to updating the weights of the Perceptron in response
    to errors during learning, the AveragedPerceptron also makes a copy of the
    old weight matrix and adds it to a running sum of all past weight matrices.
    The sums are weighted by the number of iterations that each constituent
    weight matrix survived before making an error.

    At classification time, the AveragedPerceptron uses both the current weight
    matrix and the weighted sum of past matrices to make its decision.

    This is equivalent to the "voted perceptron" algorithm described by Freund
    and Schapire (1999). The averaging approach improves on the basic Perceptron
    algorithm by providing a "large margin" approach to handling datasets that
    are not linearly separable.
    '''

    def __init__(self, event_size, outcome_size, kernel=None):
        super(AveragedPerceptron, self).__init__(
            event_size, outcome_size, kernel)
        self._iterations = numpy.zeros((outcome_size, ), 'i')
        self._survived = numpy.zeros((outcome_size, ), 'i')
        self._history = numpy.zeros(self._weights.shape, 'd')
        self._acc = numpy.zeros(self._weights.shape, 'd')

    def classify(self, event):
        return Perceptron._max_outcome(self._history, event, self._kernel)

    def learn(self, event, outcome):
        self._iterations[outcome] += 1
        pred, score = Perceptron._max_outcome(self._weights, event, self._kernel)
        if pred == outcome:
            self._survived[outcome] += 1
            return

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

    def __init__(self, event_size, outcome_size, beam_width=1000):
        '''Use sparse vectors to store weights and events.'''
        super(SparseAveragedPerceptron, self).__init__(event_size, outcome_size)
        self._weights = [defaultdict(float) for _ in xrange(outcome_size)]
        self._history = [defaultdict(float) for _ in xrange(outcome_size)]
        self._acc = [defaultdict(float) for _ in xrange(outcome_size)]
        self._kernel = SparseAveragedPerceptron._dot
        self._beam_width = beam_width

    @staticmethod
    def _dot(weights, event):
        '''Return the dot product of a sparse vector and an event.'''
        return sum(weights.get(f, 0) for f in event)

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
