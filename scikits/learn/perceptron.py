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

'''Single-layer (averaged) perceptrons.'''

from .base import BaseEstimator
from .utils import safe_asanyarray
import numpy as np


class Perceptron(BaseEstimator):
    '''Single-layer (averaged) perceptron classifier with error-driven
    learning.

    Parameters
    ----------
    averaged : bool, default False
        Train as averaged perceptron.
    learning_rate : float, default 1.
        Learning rate: multiplied into weight changes.
    n_iter : int, default 1
        Number of iterations to perform per (partial) training set.
    shuffle : bool, default False
        Randomize input sequence between iterations.
    '''

    def __init__(self, averaged=False, learning_rate=1., n_iter=1,
                 shuffle=False):
        self.averaged = averaged
        self.learning_rate = learning_rate
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
        w = self._history if self.averaged else self._weights
        bias = self._biashist if self.averaged else self._bias

        return (np.dot(X, w.T) + bias).argmax(axis=1)

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
        if self.shuffle:
            Xy = np.concatenate((X, np.atleast_2d(y).T), axis=1)
            X = Xy[:, :-1]
            y = Xy[:, -1]
        else:
            X = safe_asanyarray(X)
            y = safe_asanyarray(y)
        n_samples, n_features = X.shape
        n_labels = len(np.unique(y))

        self._weights = np.zeros((n_labels, n_features), 'd')
        self._bias = np.zeros(n_labels, 'd')
        if self.averaged:
            self._history = np.zeros(self._weights.shape, 'd')
            self._biashist = np.zeros(n_labels, 'd')

        for i in xrange(self.n_iter):
            if self.shuffle:
                np.random.shuffle(Xy)
            for j in xrange(n_samples):
                pred = (np.dot(X[j], self._weights.T) + self._bias).argmax()
                if pred != y[j]:
                    delta = self._update_weights(X[j], pred, y[j])
                    if self.averaged:
                        delta *= self.n_iter * n_samples - (i * n_samples + j)
                        self._history[pred] -= X[j] * delta
                        self._history[y[j]] += X[j] * delta
                        self._biashist[pred] -= delta
                        self._biashist[y[j]] += delta
        if self.averaged:
            self._history /= self.n_iter * n_samples

        return self

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
        self._bias = np.zeros(n_labels, 'd')
        if self.averaged:
            self._iterations = np.zeros((n_labels, ), 'i')
            self._survived = np.zeros((n_labels, ), 'i')
            self._history = np.zeros(self._weights.shape, 'd')
            self._biashist = np.zeros(n_labels, 'd')
            self._acc = np.zeros(self._weights.shape, 'd')
            self._biasacc = np.zeros(n_labels, 'd')

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
        if self.shuffle:
            Xy = np.concatenate((X, np.atleast_2d(y).T), axis=1)
            X = Xy[:, :-1]
            y = Xy[:, -1]
        else:
            X = safe_asanyarray(X)
            y = safe_asanyarray(y)

        assert X.shape[1] == self._n_features
        assert len(np.unique(y)) <= self._n_labels

        n_samples = len(y)
        assert X.shape[0] == n_samples

        lrn = self._learn_averaged if self.averaged else self._learn_ordinary

        for i in xrange(self.n_iter):
            if self.shuffle:
                np.random.shuffle(Xy)
            for j in xrange(n_samples):
                lrn(X[j], y[j])

        return self

    def _learn_averaged(self, x, label):
        '''Learn as averaged perceptron.'''

        self._iterations[label] += 1
        (pred, update) = self._learn_ordinary(x, label)
        if update:
            self._update_history(pred)
            self._update_history(label)
        else:
            self._survived[label] += 1

    def _learn_ordinary(self, x, label):
        '''Learn as ordinary perceptron.

        Parameters
        ----------
        x : array
        label : int

        Returns the prediction for x and whether an update has occurred (bool)
        '''

        # always predict as ordinary perceptron
        pred = (np.dot(x, self._weights.T) + self._bias).argmax()
        must_update = (pred != label)
        if must_update:
            self._update_weights(x, pred, label)
        return (pred, must_update)

    def _update_history(self, label):
        '''Update running average for averaged perceptron (online only).'''
        s = self._survived[label]
        if s > 0:
            acc = self._acc[label]
            acc += s * self._weights[label]
            biasacc = self._biasacc[label]
            biasacc += s * self._bias[label]
            self._history[label] = acc / self._iterations[label]
            self._biashist[label] = biasacc / self._iterations[label]
            self._survived[label] = 0

    def _update_weights(self, x, pred, label):
        # .5 since we're going to update twice
        delta = .5 * self.learning_rate
        xdelta = x * delta
        self._weights[pred] -= xdelta
        self._weights[label] += xdelta
        self._bias[pred] -= delta
        self._bias[label] += delta
        return delta
