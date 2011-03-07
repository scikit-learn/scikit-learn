"""Face recognition using PCA (eigenfaces) and SVM


This exercises has a lot of boilerplate code to extrat the most represented
faces:

TODO: once the LFW dataset loader is stable enough in scikit-learn 0.8, use it
here instead.
"""
# License: Simplified BSD

import os
from gzip import GzipFile

import sys
import numpy as np
import pylab as pl

from scikits.learn.grid_search import GridSearchCV
from scikits.learn.metrics import classification_report
from scikits.learn.metrics import confusion_matrix
from scikits.learn.pca import RandomizedPCA
from scikits.learn.svm import SVC

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


# Split the dataset into a training and test set

# TODO: define variables X_train, X_test, y_train, y_test by splitting the data
# using a 75% / 25% ratio between training set and test set


# Compute a PCA (eigenfaces) on the face dataset (treated as unlabeled
# dataset): unsupervised feature extraction / dimensionality reduction

n_components = 150

# TODO: project the training and test set on the top n_components extracted
# a truncated PCA performed on the training set


# Train a SVM classification model

# TODO: implement a grid search for the best SVM with gaussian kernel
# by searching through the following parameters

param_grid = {
 'C': [1, 5, 10, 100],
 'gamma': [0.0001, 0.001, 0.01, 0.1],
}

# TODO store the best model in a variable named 'clf'

# Uncomment the following once the above is implemented
#
## Quantitative evaluation of the model quality on the test set
#
#y_pred = clf.predict(X_test_pca)
#print classification_report(y_test, y_pred, labels=selected_target,
#                            class_names=target_names[selected_target])
#
#print confusion_matrix(y_test, y_pred, labels=selected_target)
#
#
## Qualitative evaluation of the predictions using matplotlib
#
#n_row = 3
#n_col = 4
#
#def title(y_pred, y_test, target_names, i):
#    pred_name = target_names[y_pred[i]].rsplit('_', 1)[-1]
#    true_name = target_names[y_test[i]].rsplit('_', 1)[-1]
#    return 'predicted: %s\ntrue:      %s' % (pred_name, true_name)
#
#
#pl.figure(figsize=(2 * n_col, 2.3 * n_row))
#pl.subplots_adjust(bottom=0, left=.01, right=.99, top=.95, hspace=.15)
#for i in range(n_row * n_col):
#    pl.subplot(n_row, n_col, i + 1)
#    pl.imshow(X_test[i].reshape((64, 64)), cmap=pl.cm.gray)
#    pl.title(title(y_pred, y_test, target_names, i), size=12)
#    pl.xticks(())
#    pl.yticks(())
#
#pl.show()


