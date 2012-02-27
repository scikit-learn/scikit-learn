"""Face recognition using PCA (eigenfaces) and SVM"""
# License: Simplified BSD

import os
from gzip import GzipFile

import sys
import numpy as np
import pylab as pl

from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import RandomizedPCA
from sklearn.svm import SVC

# Load dataset in memory
folder_name = sys.argv[1]
faces_filename = os.path.join(folder_name, "faces.npy.gz")
filenames_filename = os.path.join(folder_name, "face_filenames.txt")

faces = np.load(GzipFile(faces_filename))
face_filenames = [l.strip() for l in file(filenames_filename).readlines()]

# normalize each picture by centering brightness
faces -= faces.mean(axis=1)[:, np.newaxis]


# Index category names into integers suitable for scikit-learn
#
# Here we do a little dance to convert file names in integer indices
# (class indices in machine learning talk) that are suitable to be used
# as a target for training a classifier. Note the use of an array with
# unique entries to store the relation between class index and name,
# often called a 'Look Up Table' (LUT).
# Also, note the use of 'searchsorted' to convert an array in a set of
# integers given a second array to use as a LUT.
person_names = np.array([f.rsplit('_', 1)[0] for f in face_filenames])

# A unique integer per category
target_names = np.unique(person_names)

# Turn the person_names in their corresponding integer label
target = np.searchsorted(target_names, person_names)

# Subsample the dataset to restrict to the most frequent person_names
selected_target = np.argsort(np.bincount(target))[-5:]
most_frequent_mask = np.array([item in selected_target for item in target])

X = faces[most_frequent_mask]
y = target[most_frequent_mask]

n_samples, n_features = X.shape

print "Dataset size:"
print "n_samples: %d" % n_samples
print "n_features: %d" % n_features

# TASK: Split the dataset into a training and test set: 75% / 25%
split = n_samples * 3 / 4
X_train, X_test = X[:split], X[split:]
y_train, y_test = y[:split], y[split:]


# Part 1: unsupervised feature extraction

n_components = 150

# TASK: Compute a PCA (eigenfaces) on the face dataset (treated as unlabeled
# dataset)
print "Extracting the top %d eigenfaces" % n_components
pca = RandomizedPCA(n_components=n_components, whiten=True).fit(X_train)
eigenfaces = pca.components_.T.reshape((n_components, 64, 64))
# and project the input data on the eigenfaces orthonormal basis
# to do unsupervised feature extraction / dimensionality reduction
X_train_pca = pca.transform(X_train)
X_test_pca = pca.transform(X_test)


# Part 2: train a SVM classification model

param_grid = {
 'C': [1000, 10000, 100000],
 'gamma': [0.0001, 0.001, 0.01, 0.1],
}

# TASK: build a grid search for a RBF kernel SVM model named 'clf'
print "Fitting the classifier to the training set"
clf = GridSearchCV(SVC(kernel='rbf'), param_grid,
                   fit_params={'class_weight': 'auto'},
                   n_jobs=-1)
clf = clf.fit(X_train_pca, y_train)
print "Best estimator found by grid search:"
print clf.best_estimator


# Quantitative evaluation of the model quality on the test set

y_pred = clf.predict(X_test_pca)
print classification_report(y_test, y_pred, labels=selected_target,
                            target_names=target_names[selected_target])

print confusion_matrix(y_test, y_pred, labels=selected_target)


# Qualitative evaluation of the predictions using matplotlib

n_row = 3
n_col = 4


def title(y_pred, y_test, target_names, i):
    pred_name = target_names[y_pred[i]].rsplit('_', 1)[-1]
    true_name = target_names[y_test[i]].rsplit('_', 1)[-1]
    return 'predicted: %s\ntrue:      %s' % (pred_name, true_name)


pl.figure(figsize=(2 * n_col, 2.3 * n_row))
pl.subplots_adjust(bottom=0, left=.01, right=.99, top=.95, hspace=.15)
for i in range(n_row * n_col):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(X_test[i].reshape((64, 64)), cmap=pl.cm.gray)
    pl.title(title(y_pred, y_test, target_names, i), size=12)
    pl.xticks(())
    pl.yticks(())

pl.show()


