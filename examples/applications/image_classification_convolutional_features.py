print __doc__

import os, sys
import math
import cPickle
from gzip import GzipFile
from time import time

import numpy as np

from scikits.learn.grid_search import GridSearchCV
from scikits.learn.metrics import classification_report
from scikits.learn.metrics import confusion_matrix
from scikits.learn.feature_extraction.image import ConvolutionalKMeansEncoder
from scikits.learn.svm import SVC
from scikits.learn.svm import LinearSVC
from scikits.learn.preprocessing import Scaler
from scikits.learn.cluster import k_init
from scikits.learn.cluster.k_means_ import all_pairs_l2_distance_squared

X_train = []
y_train = []
n_examples_to_use = int(sys.argv[1])

folder_name = 'cifar-10-batches-py'
for filename in sorted(os.listdir(folder_name)):
    filepath = os.path.join(folder_name, filename)
    if filename.startswith('data_batch_'):
        dataset = cPickle.load(file(filepath, 'rb'))
        X_train.append(dataset['data'])
        y_train.append(dataset['labels'])
    elif filename == 'test_batch':
        dataset = cPickle.load(file(filepath, 'rb'))
        X_test = np.asarray(dataset['data'], dtype=np.float32)
        y_test = dataset['labels']
    elif filename == 'batch.meta':
        dataset = cPickle.load(file(filepath, 'rb'))
        label_neams = dataset['label_names']

X_train = np.asarray(np.concatenate(X_train), dtype=np.float32)
y_train = np.concatenate(y_train)

# reshape pictures to there natural dimension
X_train = X_train.reshape((X_train.shape[0], 3, 32, 32)).transpose(0, 2, 3, 1)
X_test = X_test.reshape((X_test.shape[0], 3, 32, 32)).transpose(0, 2, 3, 1)

print "scaling images to centered, unit variance vectors"
scaler = Scaler().fit(X_train)
X_train = scaler.transform(X_train, copy=False)
X_test = scaler.transform(X_test,copy=False)

extractor = cPickle.load(open('extractor.pkl'))

# for each image position, extract features from the entire dataset

X_train_features = extractor.transform(X_train[:n_examples_to_use])
X_test_features = extractor.transform(X_test[:n_examples_to_use])

np.save('X_train_features.npy', X_train_features)
np.save('X_test_features.npy', X_test_features)
np.save('y_train_labels.npy', y_train[:n_examples_to_use])
np.save('y_test_labels.npy', y_test[:n_examples_to_use])
