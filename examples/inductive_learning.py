"""
==============================================
Inductive Learning with Scikit Learn
==============================================

Clustering is expensive, especially when our dataset contains millions of 
datapoints. Recomputing the clusters everytime we receive some new data 
is thus in many cases, intractable. With more data, there is also the 
possibility of degrading the previous clustering. 

One solution to this problem, is to first infer the target classes using 
some unsupervised learning algorithm and then fit a classifier on the 
inferred targets, treating it as a supervised problem. This is known as 
Transductive learning.


"""
print(__doc__)

from sklearn.base import clone, BaseEstimator
from sklearn.utils.metaestimators import if_delegate_has_method

class InductiveClusterer(BaseEstimator):
    def __init__(self, clusterer, classifier):
        self.clusterer = clusterer
        self.classifier = classifier

    def fit(self, X, y=None):
        self.clusterer_ = clone(self.clusterer)
        self.classifier_ = clone(self.classifier)
        y = self.clusterer_.fit_predict(X)
        self.classifier_.fit(X, y)
        return self

    @if_delegate_has_method(delegate='classifier')
    def predict(self, X):
        return self.classifier_.predict(X)

    @if_delegate_has_method(delegate='classifier')
    def decision_function(self, X):
        return self.classifier_.decision_function(X)


# Generate Synthetic Data
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn import datasets

n_samples = 5000

colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])

blobs = datasets.make_blobs(n_samples=3*n_samples, random_state=8)
noise = np.random.rand(n_samples,2)

noise = StandardScaler().fit_transform(noise)
dataset = blobs

X, y = dataset
# normalize dataset for easier parameter selection
X = StandardScaler().fit_transform(X)
X = np.concatenate((X, noise), axis=0)
plt.scatter(X[:, 0], X[:, 1], color="black", s=2)
plt.show()

from sklearn import svm
from sklearn import cluster, datasets

# Declaring a Clustering and Classification Model
dbscan = cluster.DBSCAN(eps=0.1,min_samples=n_samples/50)
clf = svm.SVC( kernel='rbf', decision_function_shape="ovr", degree=3)

# Declaring the inductive learning model
inductiveLearner = InductiveClusterer(dbscan, clf)

inductiveLearner.fit(X)

# Inferring class on a new random dataset
X_new = StandardScaler().fit_transform(np.random.rand(n_samples*2,2))
y_pred = inductiveLearner.predict(X_new)
plt.scatter(X_new[:, 0], X_new[:, 1], color=colors[y_pred].tolist(), s=5)

plt.plot()
plt.show()
