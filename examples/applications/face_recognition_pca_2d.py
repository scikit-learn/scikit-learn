"""
===================================================
Faces recognition example using 2DPCA and Random Forests
===================================================

This dataset contains a set of face images taken between April 1992 and
April 1994 at AT&T Laboratories Cambridge. It is the Olivetti faces :

  http://cs.nyu.edu/~roweis/data/olivettifaces.mat (4MB)

.. Olivetti: http://www.cs.nyu.edu/~roweis/

Expected results for the 40 people in the dataset::


precision    recall  f1-score   support

     Face 0       1.00      0.40      0.57         5
     Face 1       1.00      1.00      1.00         5
     Face 2       0.80      0.80      0.80         5
     Face 3       1.00      1.00      1.00         5
     Face 4       1.00      0.80      0.89         5
     Face 5       0.83      1.00      0.91         5
     Face 6       0.71      1.00      0.83         5
     Face 7       0.80      0.80      0.80         5
     Face 8       1.00      1.00      1.00         5
     Face 9       1.00      0.80      0.89         5
    Face 10       1.00      1.00      1.00         5
    Face 11       0.67      0.80      0.73         5
    Face 12       1.00      0.80      0.89         5
    Face 13       1.00      1.00      1.00         5
    Face 14       0.83      1.00      0.91         5
    Face 15       1.00      1.00      1.00         5
    Face 16       1.00      1.00      1.00         5
    Face 17       0.83      1.00      0.91         5
    Face 18       1.00      1.00      1.00         5
    Face 19       1.00      0.80      0.89         5
    Face 20       0.71      1.00      0.83         5
    Face 21       1.00      0.80      0.89         5
    Face 22       1.00      1.00      1.00         5
    Face 23       0.83      1.00      0.91         5
    Face 24       0.83      1.00      0.91         5
    Face 25       0.80      0.80      0.80         5
    Face 26       1.00      1.00      1.00         5
    Face 27       1.00      1.00      1.00         5
    Face 28       1.00      1.00      1.00         5
    Face 29       1.00      1.00      1.00         5
    Face 30       1.00      1.00      1.00         5
    Face 31       1.00      0.60      0.75         5
    Face 32       1.00      1.00      1.00         5
    Face 33       1.00      1.00      1.00         5
    Face 34       1.00      0.80      0.89         5
    Face 35       1.00      0.80      0.89         5
    Face 36       0.83      1.00      0.91         5
    Face 37       1.00      1.00      1.00         5
    Face 38       0.83      1.00      0.91         5
    Face 39       0.80      0.80      0.80         5

avg / total       0.93      0.92      0.91       200



"""
print __doc__

from time import time
import numpy as np


from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import fetch_olivetti_faces
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.decomposition import PCA2D


##############################################################################
# Download the data, if not already on disk and load it as numpy arrays

olivetti_faces = fetch_olivetti_faces()

# introspect the images arrays to find the shapes
n_samples, h, w = olivetti_faces.images.shape

# we use the 2 data directly
X = olivetti_faces.images


# the label to predict is the id of the person
y = olivetti_faces.target

target_names = np.array(["Face %d" % x for x in range(0,
                        np.size(olivetti_faces.target), 1)])
# target_names = np.array(map(str,
#               (np.arange(np.size(olivetti_faces.target)).reshape(
#              olivetti_faces.target.shape))))
n_classes = target_names.shape[0]

print "Total dataset size:"
print "n_samples: %d" % n_samples
print "features: %d X %d" % (h, w)
print "n_classes: %d" % n_classes


##############################################################################
# Split into a training set and a test set using a stratified shuffle split

# split into a training and testing set
sss = StratifiedShuffleSplit(y, 1, test_size=0.5, random_state=0)
for train_index, test_index in sss:
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]


##############################################################################
# Compute a 2DPCA (eigenfaces) on the face dataset (treated as unlabeled
# dataset): unsupervised feature extraction / dimensionality reduction

n_row_components = 0.9
n_column_components = 0.9

print("Reducing images by keeping %d %% of the variance on the row"
      % (100 * n_row_components))
print "and %d %% of the variance on the column from %d faces" % (
      100 * n_column_components, X_train.shape[0])
t0 = time()
pca = PCA2D(n_row_components=n_row_components, n_column_components=
            n_column_components, row_whiten=True,
            column_whiten=True, epsilon=0).fit(X_train)
print "done in %0.3fs" % (time() - t0)


print "Projecting the input data on the eigenfaces orthonormal basis"
t0 = time()
X_train_pca = pca.transform(X_train)
X_test_pca = pca.transform(X_test)
print "done in %0.3fs" % (time() - t0)

# Converting the images to 1D for classification
X_train_pca = X_train_pca.reshape(X_train_pca.shape[0], X_train_pca.shape[1] *
                                  X_train_pca.shape[2])

X_test_pca = X_test_pca.reshape(X_test_pca.shape[0], X_test_pca.shape[1] *
                                X_test_pca.shape[2])

X_test = X_test.reshape(X_test.shape[0], h * w)
X_train = X_train.reshape(X_train.shape[0], h * w)


##############################################################################
# Train a Random Forests classification model

print "Fitting the classifier to the training set"
t0 = time()
param_grid = {
    'n_estimators': [10, 50, 90, 300, 500, 700, 1000, 1500, 2000],
}
clf = GridSearchCV(RandomForestClassifier(), param_grid)
clf = clf.fit(X_train_pca, y_train)
print "done in %0.3fs" % (time() - t0)
print "Best estimator found by grid search:"
print clf.best_estimator_


##############################################################################
# Quantitative evaluation of the model quality on the test set

print "Predicting the people names on the testing set"
t0 = time()
y_pred = clf.predict(X_test_pca)
print "done in %0.3fs" % (time() - t0)

print classification_report(y_test, y_pred, target_names=target_names)
print confusion_matrix(y_test, y_pred, labels=range(n_classes))

print "Size of the components retained %d X %d" % (pca.n_row_components_,
                                                   pca.n_column_components_)
