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
from .utils import safe_asanyarray
from .utils.extmath import safe_sparse_dot
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
        return (safe_sparse_dot(weights, x) + self._alpha) ** self._degree


class Perceptron(BaseEstimator):
    '''Single-layer (averaged) perceptron classifier with error-driven
    learning.

    If the `averaged` parameter is true, this classifier maintains a running
    sum of past weight matrices, weighted by the number of iterations that each
    weight matrix survived before making an error. At classification time, the
    averaged perceptron uses both the current weight matrix and the weighted
    sum of past matrices to make its decision.

    This is an approximation to the "voted perceptron" algorithm described by
    Freund and Schapire (1999). The averaging approach improves on the basic
    perceptron algorithm by providing a "large margin" approach to handling
    datasets that are not linearly separable.

    Parameters
    ----------
    averaged : bool, optional, default False
        Train as averaged perceptron.
    kernel : callable, optional
        Kernel function for mapping non-linear data.
    learning_rate : float, optional
        Learning rate: multiplied into weight changes, default 1.
    n_iter : int, optional, default 1
        Number of iterations to perform per (partial) training set.
    shuffle : bool, optional, default False
        Randomize input sequence between iterations.
    '''

    def __init__(self, averaged=False, kernel=None, learning_rate=1.,
                 n_iter=1, shuffle=False):
        self.averaged = averaged
        self.learning_rate = learning_rate
        self.kernel = kernel or safe_sparse_dot
        self.n_iter = n_iter
        self.shuffle = shuffle

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

        w = self._history if self.averaged else self._weights
        for i, x in enumerate(X):
            y[i] = self._classify(x, w)
        return y

    def _classify(self, x, weights):
        '''Classify x into one of our possible labels.
        Returns the label for the most likely outcome.
        '''
        return self._max_outcome(weights, x)[0]

    def fit(self, X, y):
        """Fit classifier according to inputs X with labels y (batch learning)

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self
        """
        X = safe_asanyarray(X)
        y = safe_asanyarray(y)
        n_samples, n_features = X.shape
        n_labels = len(np.unique(y))

        return self.partial_setup(n_features, n_labels).partial_fit(X, y)

    def partial_setup(self, n_features, n_labels):
        """Setup classifier for online learning.

        Must be run before partial_fit.

        Parameters
        ----------
        n_features : int
            Number of features in (length of) training vectors.
        n_labels : int
            Number of output labels.

        Returns
        -------
        self
        """
        self._n_features = n_features
        self._n_labels = n_labels

        self._weights = np.zeros((n_labels, n_features), 'd')
        if self.averaged:
            self._iterations = np.zeros((n_labels, ), 'i')
            self._survived = np.zeros((n_labels, ), 'i')
            self._history = np.zeros(self._weights.shape, 'd')
            self._acc = np.zeros(self._weights.shape, 'd')

        return self

    def partial_fit(self, X, y):
        """Partially fit classifier according to inputs X with labels y.

        Can be run multiple times to perform online learning. Must be run
        after partial_setup.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self
        """
        if self.shuffle:            # copy X and y
            X = np.array(X)
            y = np.array(y)
        else:
            X = safe_asanyarray(X)
            y = safe_asanyarray(y)

        assert X.shape[1] == self._n_features
        assert len(np.unique(y)) <= self._n_labels

        n_samples = len(y)
        assert X.shape[0] == n_samples
        order = range(n_samples) if self.shuffle else xrange(n_samples)

        for i in xrange(self.n_iter):
            if self.shuffle:
                np.random.shuffle(order)
            for j in order:
                self._learn(X[j], y[j])

        return self

    def _learn(self, x, label):
        lrn = self._learn_averaged if self.averaged else self._learn_ordinary
        lrn(x, label)

    def _learn_averaged(self, x, label):
        '''Learn as averaged perceptron.'''

        self._iterations[label] += 1
        (pred, update) = self._learn_ordinary(x, label)
        if update:
            self._update_history(pred)
            self._update_history(label)
        else:
            self._survived[label] += 1
        #return (pred, update)

    def _learn_ordinary(self, x, label):
        '''Learn as ordinary perceptron.

        Parameters
        ----------
        x : array
        label : int

        Returns the prediction for x and whether an update has occurred (bool)
        '''

        # always predict as ordinary perceptron
        pred = self._classify(x, self._weights)
        rate = .5 * self.learning_rate      # we're going to update twice
        must_update = (pred != label)
        if must_update:
            self._update_weights(pred, x, -rate)
            self._update_weights(label, x, rate)
        return (pred, must_update)

    def _update_history(self, label):
        '''Update the history for a particular class (averaged perceptron).'''
        s = self._survived[label]
        if s > 0:
            acc = self._acc[label]
            acc += s * self._weights[label]
            self._history[label] = acc / self._iterations[label]
            self._survived[label] = 0

    def _update_weights(self, label, x, delta):
        '''Update the weights for an index based on an event.

        label: The index of a weight vector to update.
        x: An event vector to use for the update.
        delta: Weight the update by this value.
        '''
        self._weights[label] += delta * x

    def _max_outcome(self, weights, x):
        '''Return the maximum scoring label for a set of weights and its score

        weights: An array of weight vectors, one per label.
        x: An event vector.
        '''
        scores = np.array([self.kernel(w, x) for w in weights])
        index = scores.argmax()
        # Return the score as a probability computed locally among our possible
        # outcomes. Subtracting the maximum raw score before exponentiating
        # stabilizes the result numerically, while also encouraging low-scoring
        # values to err by squashing towards 0.
        return index, 1.0 / np.exp(scores - scores[index]).sum()


class SparseDot(object):
    '''Kernel for SparsePerceptron'''
    @staticmethod
    def __call__(weights, x):
        '''Return the dot product of a sparse vector and an event.'''
        return sum(weights.get(f, 0) for f in x)


class SparsePerceptron(Perceptron):
    '''(Averaged) perceptron using sparse weight vectors.

    This class achieves sparseness by assuming that all events are
    binary-valued. To use this class with a non-binary event space, you will
    need to manually chop up the event space into binary events. For example,
    suppose event 'foo' can take on N possible discrete values ; you might
    convert 'foo' into N separate events, 'foo=1' .. 'foo=N'.

    In addition, this class may be configured to retain only the top-weighted
    features for each label, making this a sort of beam search version of the
    general Perceptron. This reduces the space requirements for the classifier,
    at the cost of lower accuracy.

    Parameters
    ----------
    averaged : bool, optional, default False
        Train as averaged perceptron.
    beam_width : int, optional
        Width of beam for beam search. This parameter is given to fit
        because it implies pruning of the feature set.
    learning_rate : float, optional
        Learning rate: multiplied into weight changes, default 1.
    n_iter : int, optional, default 1
        Number of iterations to perform per (partial) training set.
    shuffle : bool, optional, default False
        Randomize input sequence between iterations.
    '''

    def __init__(self, averaged=False, beam_width=None, learning_rate=1.,
                 n_iter=1, shuffle=False):
        '''Use sparse vectors to store weights and events.'''
        super(SparsePerceptron, self).__init__(averaged=averaged,
                                               kernel=SparseDot(),
                                               learning_rate=learning_rate,
                                               n_iter=n_iter,
                                               shuffle=shuffle)
        self.beam_width=beam_width

    def fit(self, X, y):
        """Fit classifier according to inputs X with labels y

        Can be run multiple times for online learning.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self
        """
        X = safe_asanyarray(X)
        y = safe_asanyarray(y)
        n_samples, n_features = X.shape
        n_labels = len(np.unique(y))

        return self.partial_setup(n_features, n_labels).partial_fit(X, y)

    def partial_setup(self, n_features, n_labels):
        """Setup classifier for online learning.

        Must be run before partial_fit.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self
        """
        self._n_features = n_features
        self._n_labels = n_labels

        self._weights = [defaultdict(float) for i in xrange(n_labels)]
        if self.averaged:
            self._history = [defaultdict(float) for i in xrange(n_labels)]
            self._acc = [defaultdict(float) for i in xrange(n_labels)]
            self._iterations = np.zeros((n_labels, ), 'i')
            self._survived = np.zeros((n_labels, ), 'i')

        return self

    def _update_history(self, label):
        s = self._survived[label]
        if s > 0:
            i = self._iterations[label]
            a = self._acc[label]
            h = self._history[label]
            for f, w in self._weights[label].iteritems():
                a[f] += s * w
                h[f] = a[f] / i
            self._prune(h)
            self._survived[label] = 0

    def _update_weights(self, label, x, delta):
        w = self._weights[label]
        for f in x:
            w[f] += delta
        self._prune(w)

    def _prune(self, weights):
        '''Prune the weights in a sparse vector to our beam width.'''
        if self.beam_width is not None \
          and len(weights) < 1.3 * self.beam_width:
            fws = sorted(weights.iteritems(), key=lambda x: -abs(x[1]))
            for f, _ in fws[self.beam_width:]:
                del weights[f]
