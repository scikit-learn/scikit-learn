"""
=========================================================================
Comparing randomized search and grid search for hyperparameter estimation
=========================================================================

This example shows how the :class:`sklearn.model_selection.UnsupervisedRandomizedSearch`
object can be used to optimize the hyperparameters of an unsupervised learning 
algorithm.

The set of hyperparameters should be evaluated on a scoring function 
appropriate to the type of algorithm. In this example t-distributed Stochastic 
Neighbor Embedding (t-SNE) is used to visualise the structure of 3-dimensional 
data in two dimensions. The scoring function is "trustworthiness" which 
expresses to what extent the local structure is retained in lower dimensions.

"""
print(__doc__)

import time

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from scipy.stats import randint as sp_randint

from sklearn import datasets
from sklearn.manifold.t_sne import TSNE, trustworthiness
from sklearn.model_selection import UnsupervisedRandomizedSearch


def trustworthiness_scorer(estimator, X, y=None):
    return trustworthiness(X, estimator.embedding_)


n_points = 250
X, color = datasets.samples_generator.make_s_curve(n_points, random_state=999)
n_neighbors = 10
n_components = 2

fig = plt.figure(figsize=(8, 10))
plt.suptitle("Manifold Learning with %i points, %i neighbors"
             % (n_points, n_neighbors), fontsize=14)

ax = fig.add_subplot(221, projection='3d')
plt.scatter(X[:, 0], X[:, 2], c=color, cmap=plt.cm.Spectral)

tsne = TSNE(random_state=0)

np.random_seed = 0
learning_rate = sp_randint(100, 1000)
perplexity = sp_randint(1, 100)

scorer = trustworthiness_scorer
urs = UnsupervisedRandomizedSearch(
    tsne, {'learning_rate': learning_rate,
           'perplexity': perplexity},
    scoring=scorer,
    n_iter=15,
    random_state=2)

t0 = time.time()
urs.fit(X)
t1 = time.time()
print("Randomized search t-SNE: %.2g sec" % (t1 - t0))
print(urs.best_params_)

t0 = time.time()
Y = urs.best_estimator_.fit_transform(X)
t1 = time.time()

ax = fig.add_subplot(222)
plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=plt.cm.Spectral)
plt.title("t-SNE (%.2g sec)" % (t1 - t0))
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

plt.show()