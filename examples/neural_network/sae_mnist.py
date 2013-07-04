from sklearn.datasets import fetch_mldata
import random
import numpy as np
from sklearn.linear_model import SGDClassifier
from sklearn.neural_network import SAE

mnist = fetch_mldata('MNIST original')
X, y = mnist.data, mnist.target
random.seed(100)
indices = np.array(random.sample(range(70000), 1000))
X, y = X[indices].astype('float64'), y[indices]
# For SAE, feature values in the range [0, 1] is necessary
X = X / 255
sae = SAE(
    optimization_method='l-bfgs-b',
    verbose=True,
    max_iter=400,
    n_hidden=191,
    random_state=3)
sae_features = sae.fit_transform(X)
clf = SGDClassifier(random_state=3)
clf.fit(X, y)
#Should get a score of 0.943
print 'SGD on raw pixels score: ', clf.score(X, y)
clf.fit(sae_features, y)
#Should get a score of 0.982
print 'SGD on extracted features score: ', clf.score(sae_features, y)


