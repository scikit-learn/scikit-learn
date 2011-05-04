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

'''Single-layer (averaged) perceptrons.

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
'''

from .base import BaseEstimator
from collections import defaultdict
import numpy as np


def RadialBasisKernel(gamma):
    '''This kernel returns the rbf between the input vectors.

    Defined as exp(-gamma * delta**2)'''

    def __init__(self, gamma):
        self._minus_gamma = -gamma

    def __call__(self, weights, x):
        delta = np.linalg.norm(x - weights)
        return np.exp(self._minus_gamma * delta * delta)


class PolynomialKernel(object):
    '''This kernel returns the dot product raised to some power.

    Also known as homogeneous kernel.

    Parameters
    ----------
    degree : float
        Degree of polynomial. Should be >1.
    alpha : float, optional
        Added to dot product before exponentiation, default 1.
    '''
    def __init__(self, degree, alpha=1.0):
        self._degree = degree
        self._alpha = alpha

    def __call__(self, weights, x):
        return (np.dot(weights, x) + self._alpha) ** self._degree


class Perceptron(BaseEstimator):
    '''Classifier implementing single-layer perceptron with error-driven
    learning.

    Parameters
    ----------
    kernel : callable, optional
        Kernel function for mapping non-linear data.
    '''

    def __init__(self, kernel=None, learning_rate=1.):
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
            y[i] = self._classify(x)
        return y

    def _classify(self, x):
        '''Classify an event into one of our possible labels.

        x: A np vector of the dimensionality of our input space.

        Returns the label for the most likely outcome.
        '''
        return self._max_outcome(self._weights, x)[0]

    def _fit(self, X, y, learning_rate, n_iter, shuffle, averaged):
        if shuffle:                 # copy X and y
            X = np.array(X)
            y = np.array(y)
        else:
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

        self._run_iters(n_iter, X, y, learning_rate, shuffle)

        return self

    def _run_iters(self, n, X, y, learning_rate, shuffle):
        '''Run n iterations of training'''
        n_samples = len(y)
        order = range(n_samples) if shuffle else xrange(n_samples)

        for i in xrange(n):
            if shuffle:
                np.random.shuffle(order)
            for j in order:
                self._learn(X[j], y[j], learning_rate)

    def fit(self, X, y, learning_rate=1., n_iter=1, shuffle=False):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.
        learning_rate : float, optional
            Learning rate: multiplied into weight changes, default 1.
        n_iter : int, optional, default 1
            Number of iterations to perform.
        shuffle : bool, optional, default False
            Randomize input sequence between iterations.

        Returns
        -------
        self
        """

        return self._fit(X, y, learning_rate=learning_rate, n_iter=n_iter,
                         shuffle=shuffle, averaged=False)

    def _learn(self, x, label, rate):
        '''Adjust the hyperplane based on a classification attempt.

        Parameters
        ----------
        x : array
        label : int
        rate : float

        Returns the prediction for x and whether an update has occurred (bool)
        '''

        # always predict as ordinary perceptron
        pred = Perceptron._classify(self, x)
        rate = .5 * rate    # we're going to update twice
        must_update = (pred != label)
        if must_update:
            self._update_weights(pred, x, -rate)
            self._update_weights(label, x, rate)
        return (pred, must_update)

    def _update_weights(self, class_index, x, delta):
        '''Update the weights for an index based on an event.

        class_index: The index of a weight vector to update.
        x: An event vector to use for the update.
        delta: Weight the update by this value.
        '''
        self._weights[class_index] += delta * x

    def _max_outcome(self, weights, x):
        '''Return the maximum scoring label for a set of weights and its score

        weights: An array of weight vectors, one per label.
        x: An event vector.
        '''
        scores = np.array([self._kernel(w, x) for w in weights])
        index = scores.argmax()
        # Return the score as a probability computed locally among our possible
        # outcomes. Subtracting the maximum raw score before exponentiating
        # stabilizes the result numerically, while also encouraging low-scoring
        # values to err by squashing towards 0.
        return index, 1.0 / np.exp(scores - scores[index]).sum()


class AveragedPerceptron(Perceptron):
    '''Classifier based on a weighted sum of perceptrons.

    This perceptron algorithm performs similarly to the basic perceptron when
    learning from labeled data: whenever the predicted label for an event
    differs from the true label, the weights of the perceptron are updated to
    classify this new event correctly.

    However, in addition to updating the weights of the perceptron in response
    to errors during learning, the AveragedPerceptron also makes a copy of the
    old weight matrix and adds it to a running sum of all past weight matrices.
    The sums are weighted by the number of iterations that each constituent
    weight matrix survived before making an error.

    At classification time, the AveragedPerceptron uses both the current weight
    matrix and the weighted sum of past matrices to make its decision.

    This is an approximation to the "voted perceptron" algorithm described by
    Freund and Schapire (1999). The averaging approach improves on the basic
    perceptron algorithm by providing a "large margin" approach to handling
    datasets that are not linearly separable.
    '''

    def __init__(self, kernel=None):
        super(AveragedPerceptron, self).__init__(kernel)

    def _classify(self, x):
        return self._max_outcome(self._history, x)[0]

    def fit(self, X, y, learning_rate=1., n_iter=1, shuffle=False):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.
        learning_rate : float, optional
            Learning rate: multiplied into weight changes, default 1.
        n_iter : int, optional
            Number of iterations to perform. Default 1.
        shuffle : bool, optional, default False
            Randomize input sequence between iterations.

        Returns
        -------
        self
        """

        return self._fit(X, y, learning_rate=learning_rate, n_iter=n_iter,
                         shuffle=shuffle, averaged=True)

    def _learn(self, x, label, rate):
        self._iterations[label] += 1
        (pred, update) = super(AveragedPerceptron, self)._learn(x, label, rate)
        if update:
            self._update_history(pred)
            self._update_history(label)
        else:
            self._survived[label] += 1
        return (pred, update)

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
    def __call__(weights, x):
        '''Return the dot product of a sparse vector and an event.'''
        return sum(weights.get(f, 0) for f in x)


class SparseAveragedPerceptron(AveragedPerceptron):
    '''A voted perceptron using sparse vectors for storing weights.

    This class achieves sparseness by assuming that all events are
    binary-valued. To use this class with a non-binary event space, you will
    need to manually chop up the event space into binary events. For example,
    suppose event 'foo' can take on N possible discrete values ; you might
    convert 'foo' into N separate events, 'foo=1' .. 'foo=N'.

    In addition, this class may be configured to retain only the top-weighted
    features for each label, making this a sort of beam search version of the
    general Perceptron. This reduces the space requirements for the classifier,
    at the cost of lower accuracy.
    '''

    def __init__(self):
        '''Use sparse vectors to store weights and events.'''
        super(SparseAveragedPerceptron, self).__init__(kernel=SparseDot())

    def fit(self, X, y, beam_width=None, learning_rate=1., n_iter=1,
            shuffle=False):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values.

        beam_width : int, optional
            Width of beam for beam search. This parameter is given to fit
            because it implies pruning of the feature set.
        learning_rate : float, optional
            Learning rate: multiplied into weight changes, default 1.
        n_iter : int, optional, default 1
            Number of iterations to perform.
        shuffle : bool, optional, default False
            Randomize input sequence between iterations.

        Returns
        -------
        self
        """

        if shuffle:                 # copy X and y
            X = np.array(X)
            y = np.array(y)
        else:
            X = np.asanyarray(X)
            y = np.asanyarray(y)

        n_samples, n_features = X.shape
        n_labels = len(np.unique(y))

        self._beam_width = beam_width

        self._weights = [defaultdict(float) for i in xrange(n_labels)]
        self._history = [defaultdict(float) for i in xrange(n_labels)]
        self._acc = [defaultdict(float) for i in xrange(n_labels)]
        self._iterations = np.zeros((n_samples, ), 'i')
        self._survived = np.zeros((n_samples, ), 'i')

        self._run_iters(n_iter, X, y, learning_rate, shuffle)

        for i in xrange(n_iter):
            for j in xrange(n_samples):
                self._learn(X[j], y[j], learning_rate)

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

    def _update_weights(self, class_index, x, delta):
        w = self._weights[class_index]
        for f in x:
            w[f] += delta
        self._prune(w)

    def _prune(self, weights):
        '''Prune the weights in a sparse vector to our beam width.'''
        beamwidth = self._beam_width
        if beamwidth is not None and len(weights) < 1.3 * beamwidth:
            fws = sorted(weights.iteritems(), key=lambda x: -abs(x[1]))
            for f, _ in fws[beamwidth:]:
                del weights[f]
