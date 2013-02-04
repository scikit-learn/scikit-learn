"""
===================================================
Faces recognition benchmark using EigenFace and LaplacianFace
===================================================
"""

import logging

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.datasets import fetch_lfw_people, olivetti_faces
from sklearn.decomposition import PCA, RandomizedPCA
from sklearn.manifold.lpp import LocalityPreservingProjection as LPP
from sklearn.neighbors.classification import KNeighborsClassifier
from sklearn.cross_validation import train_test_split

class _BaseFaceRecognizer(BaseEstimator, ClassifierMixin, TransformerMixin):
    def __init__(self, transformer=None, classifier=None):
        if transformer:
            self.transformer_ = transformer
        else:
            self.transformer_ = PCA(n_components=0.8)
        
        if classifier:
            self.classifier_ = classifier
        else:
            self.classifier_ = KNeighborsClassifier(1)


    def transform(self, faces):
        return self.transformer_.transform(faces)

    def fit(self, faces, labels):
        self._faces = faces.copy()

        try:
            self.transformer_.fit(faces, labels)
        except TypeError:
            self.transformer_.fit(faces)
            
        features = self.transform(faces)
        self._features = features.copy()
        self._labels = list(labels)
        self.classifier_.fit(features, labels)
        return self
           
    def predict(self, faces):
        features = self.transform(faces)
        return self.classifier_.predict(features)
    
class EigenFace(_BaseFaceRecognizer):
    def __init__(self, n_components=0.8, copy=True, whiten=False):
        transformer = RandomizedPCA(n_components, copy, whiten)
        _BaseFaceRecognizer.__init__(self, transformer=transformer)

class LaplacianFace(_BaseFaceRecognizer):
    def __init__(self, n_neighbors=None, n_components=2, kernel_param=10.0):
        transformer = LPP(n_neighbors=n_neighbors,
                n_components=n_components, kernel_param=kernel_param)
        _BaseFaceRecognizer.__init__(self, transformer=transformer)
        

#-----------------main-----------------------

# Display progress logs on stdout
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


###############################################################################
# Download the data, if not already on disk and load it as numpy arrays

faces = fetch_lfw_people(min_faces_per_person=70, resize=0.4)
#faces = olivetti_faces.fetch_olivetti_faces()

# introspect the images arrays to find the shapes (for plotting)
n_samples, h, w = faces.images.shape

# fot machine learning we use the 2 data directly (as relative pixel
# positions info is ignored by this model)
X = faces.data
n_features = X.shape[1]

# the label to predict is the id of the person
y = faces.target
target_names = faces.target

print "Total dataset size:"
print "n_samples: %d" % n_samples
print "n_features: %d" % n_features


scores = [[], []]
n_components = 50
clfs = [EigenFace(n_components),
        LaplacianFace(5, n_components,kernel_param="auto")]


for i in xrange(10):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25)
    for j in xrange(len(clfs)):
        clfs[j].fit(X_train, y_train)
        scores[j].append(clfs[j].score(X_test, y_test))

mean = np.mean(scores, 1) * 100
print "Precision"              
print "EigenFace : %.1f%%, LaplacianFace : %.1f%%" % (mean[0], mean[1])
