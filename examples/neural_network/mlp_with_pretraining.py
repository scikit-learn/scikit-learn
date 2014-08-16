"""
=======================================================================
Pre-training Multi-layer Perceptron using Restricted Boltzmann Machines
=======================================================================

This compares the performance of multi-layer perceptron (MLP) with and without
pre-training using Restricted Boltzmann Machines. MLP without pre-training
initializes the coefficient and intercept components using scaled, random
distribution. For pre-training, an RBM trains on the dataset and the resultant
parameters are given to MLP as initial coefficient and intercept parameters.
This example justifies the hypothesis that pre-training sometimes allow MLP
to converge in a better local minima. However, it is not always the case that
pre-training is beneficial. In fact, as the labeled training set grows large
the less beneficial pre-training is likely to be. Moreover, it can even result
in a performance less than if the weights were randomly initialized.
"""

from __future__ import print_function

print(__doc__)

import numpy as np

from sklearn.cross_validation import train_test_split
from sklearn.datasets import load_digits
from sklearn.neural_network import BernoulliRBM, MultilayerPerceptronClassifier


random_state = 0
n_hidden = 50

# Load Data
digits = load_digits()
X, y = np.asarray(digits.data, 'float32'), digits.target
X = (X - np.min(X, 0)) / (np.max(X, 0) + 0.0001)  # 0-1 scaling
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size=0.2,
                                                    random_state=random_state)
# Cross-validate multi-layer perceptron without pre-training
mlp = MultilayerPerceptronClassifier(n_hidden=n_hidden,
                                     random_state=random_state)

mlp.fit(X_train, y_train)
score_without_pretraining = mlp.score(X_test, y_test)

# Cross-validate multi-layer perceptron with rbm pre-training
rbms = [BernoulliRBM(n_components=n_hidden, random_state=random_state,
                     learning_rate=0.001, n_iter=50),
        BernoulliRBM(n_components=10, random_state=random_state,
                     learning_rate=0.001, n_iter=50)]

mlp = MultilayerPerceptronClassifier(n_hidden=n_hidden,
                                     random_state=random_state,
                                     warm_start=rbms)

# Train multi-layer perceptron and get the score
mlp.fit(X_train, y_train)
score_with_pretraining = mlp.score(X_test, y_test)

print("Testing accuracy of mlp without pretraining: %.3f" %
      (score_without_pretraining))
print("Testing accuracy of mlp with pretraining: %.3f" %
      (score_with_pretraining))
