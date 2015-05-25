"""
===========================================================================
Impact of increasing the number of hidden neurons in random neural networks
===========================================================================

This illustrates how the random neural network behaves when increasing
the number of hidden neurons. Larger number of hidden neurons increases
training score, but might reduce the testing score as a result of overfitting.

The example generates a plot showing the how training and testing scores change
with the number of hidden neurons on a small dataset.

"""
print(__doc__)


# Author: Issam H. Laradji
# License: BSD 3 clause

import numpy as np

from sklearn.neural_network import RandomBasisFunction
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline

###############################################################################
# Generate sample data
n_samples_train, n_samples_test = 100, 50
n_features = 50

np.random.seed(0)

coef = np.random.randn(n_features)
X = np.random.randn(n_samples_train + n_samples_test, n_features)
y = np.dot(X, coef)

# Split train and test data
X_train, X_test = X[:n_samples_train], X[n_samples_train:]
y_train, y_test = y[:n_samples_train], y[n_samples_train:]

###############################################################################
# Compute train and test errors
n_hidden_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

rnn = make_pipeline(RandomBasisFunction(), Ridge(alpha=0))

train_errors = list()
test_errors = list()
for n_hidden in n_hidden_list:
    rnn.set_params(randombasisfunction__n_activated_features=n_hidden)
    rnn.fit(X_train, y_train)

    train_errors.append(rnn.score(X_train, y_train))
    test_errors.append(rnn.score(X_test, y_test))


###############################################################################
# Plot results functions

import pylab as pl

pl.plot(n_hidden_list, train_errors, label='Train')
pl.plot(n_hidden_list, test_errors, label='Test')

pl.legend(loc='lower left')
pl.title("Random neural network on training vs. testing scores")
pl.xlabel('number of neurons in the hidden layer ')
pl.ylabel('The $R^2$ score')

pl.ylim([0.1, 1.01])

pl.show()
