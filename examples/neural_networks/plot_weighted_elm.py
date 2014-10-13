"""
==================================
Weighted Extreme Learning Machines
==================================

Plot decision functions of extreme learning machines with different class
weights. Assigning larger weight to a class will push the decision function
away from that class to have more of its samples correctly classified.
Such scheme is useful for imbalanced data so that underrepresented classes
are emphasized and therefore not ignored by the classifier.

"""
print(__doc__)

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neural_network import ELMClassifier


def plot_decision_function(clf, axis, title):
    xx, yy = np.meshgrid(np.linspace(-5, 5, 500),
                         np.linspace(-5, 5, 500))

    # plot the decision function for each datapoint on the grid
    Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    axis.imshow(Z, interpolation='nearest',
                extent=(xx.min(), xx.max(), yy.min(), yy.max()),
                aspect='auto', origin='lower', cmap=plt.cm.RdBu_r, vmin=-3,
                vmax=3)
    axis.contour(xx, yy, Z, levels=[0], linewidths=2, linetypes='--',
                 color='b')
    axis.scatter(X[:, 0], X[:, 1], s=20, c=y, cmap=plt.cm.Paired)
    axis.axis('off')
    axis.set_title(title)


rng = np.random.RandomState(0)
n_samples_1 = 1000
n_samples_2 = 100
X = np.r_[1.5 * rng.randn(n_samples_1, 2),
          0.5 * rng.randn(n_samples_2, 2) + [2, 2]]
y = [0] * (n_samples_1) + [1] * (n_samples_2)

n_hidden = 100

activation = 'tanh'
clf_weightless = ELMClassifier(n_hidden=n_hidden, activation=activation,
                               C=10e5, class_weight=None)
clf_weight_auto = ELMClassifier(n_hidden=n_hidden, C=10e5,
                                activation=activation, class_weight='auto')
clf_weight_1000 = ELMClassifier(n_hidden=n_hidden, C=10e5,
                                activation=activation, class_weight={1: 1000})

clf_weightless.fit(X, y)
clf_weight_auto.fit(X, y)
clf_weight_1000.fit(X, y)

_, axes = plt.subplots(1, 3, figsize=(10, 4))

plot_decision_function(clf_weightless, axes[0], 'class_weight=None')
plot_decision_function(clf_weight_auto, axes[1], 'class_weight="auto"')
plot_decision_function(clf_weight_1000, axes[2], 'class_weight={1:1000}')

plt.subplots_adjust(left=0, bottom=0, right=1, top=0.9, wspace=0, hspace=0)

plt.show()
