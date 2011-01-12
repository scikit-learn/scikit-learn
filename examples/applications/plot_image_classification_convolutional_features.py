"""
==================================================================
Image classification example using convolutional features and SVMs
==================================================================

The dataset used in this example is the CIFAR-10 dataset:

  http://www.cs.toronto.edu/~kriz/cifar.html

This implementation uses an unsupervised feature extraction scheme
to extract a dictionnary of 400 small (6, 6)-shaped filters to be
convolationally applied to the input images as described in:

  An Analysis of Single-Layer Networks in Unsupervised Feature Learning
  Adam Coates, Honglak Lee and Andrew Ng. In NIPS*2010 Workshop on
  Deep Learning and Unsupervised Feature Learning.

Expected results:

  TODO

"""
print __doc__

import os
import sys
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
from scikits.learn.linear_model import SGDClassifier


def load_cifar10():
    # Download the data, if not already on disk
    url = "http://www.cs.toronto.edu/~kriz/cifar-10-python.tar.gz"
    archive_name = url.rsplit('/', 1)[1]
    folder_name = 'cifar-10-batches-py'

    if not os.path.exists(folder_name):
        if not os.path.exists(archive_name):
            import urllib
            print "Downloading data, please Wait (163MB)..."
            print url
            opener = urllib.urlopen(url)
            open(archive_name, 'wb').write(opener.read())
            print

        import tarfile
        print "Decompressiong the archive: " + archive_name
        tarfile.open(archive_name, "r:gz").extractall()
        print

    ################################################################################
    # Load dataset in memory

    X_train = []
    y_train = []

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

    #n_samples = X_train.shape[0]

    # reshape pictures to there natural dimension
    X_train = X_train.reshape((X_train.shape[0], 3, 32, 32)).transpose(0, 2, 3, 1)
    X_test = X_test.reshape((X_test.shape[0], 3, 32, 32)).transpose(0, 2, 3, 1)

    ## convert to graylevel images for now
    #X_train = X_train.mean(axis=-1)
    #X_test = X_test.mean(axis=-1)
    #pl.imshow(X_train[0], interpolation='nearest'); pl.show()

    # scale dataset
    print "scaling images to centered, unit variance vectors"
    scaler = Scaler().fit(X_train)
    X_train = scaler.transform(X_train,copy=False)
    X_test = scaler.transform(X_test,copy=False)

    return X_train, y_train, X_test, y_test

def train_convolutional_kmeans(cifar10, n_centers=400, n_components=30, n_drop_components=4,
        save_images=True, n_images_to_train_from=10000, n_EM_steps=30):
    X_train, y_train, X_test, y_test = cifar10

    extractor = ConvolutionalKMeansEncoder(
        n_centers=n_centers, # kmeans centers: convolutional filters
        patch_size=8,# size of the side of one filter
        whiten=True, # perform whitening or not before kmeans
        n_drop_components=n_drop_components,  #n  of leading eigenvecs to ignore
        n_components=n_components, # number of singular vectors to keep when whitening
        max_iter=n_EM_steps, 
        n_init=1,   # take best fit of this many trials
        #kmeans_init_algo=lambda X,k,rng:k_init(X,k,rng=rng, n_samples_max=2000),
        kmeans_init_algo='random', #I'm guessing smart init in high dimensions irrelevant
        verbose=1)


    print "training convolutional whitened kmeans feature extractor..."
    t0 = time()
    # restrict training size for faster runtime as a demo
    extractor.fit(X_train[:n_images_to_train_from])
    print "done in %0.3fs" % (time() - t0)

    if extractor.whiten:
        vr = extractor.pca.explained_variance_ratio_
        print "explained variance ratios for %d kept PCA components:" % vr.shape[0]
        print vr
        if extractor.n_drop_components:
            print ".. but DROPPING variance in the leading %i components" % (
                    extractor.n_drop_components,)
    print "kmeans remaining inertia: %0.3fe6" % (extractor.inertia_ / 1e6)


    if save_images:
        extractor.tile_kernels(scale_each=True).save('kernels.png')
        extractor.tile_patches(scale_each=True).save('patches.png')
        extractor.tile_patches_unpca(scale_each=True).save('patches_unpca.png')

    return extractor


################################################################################
# DRIVERS MEANT TO BE CALLED FROM COMMAND LINE (SEE __name__=='__main__' BELOW)
#
def train_kmeans(save_extractor='extractor.pkl', n_centers=400, n_components=30,
        n_drop_components=4, n_EM_steps=30):
    print 'Training convolutional k-means'
    # Qualitative evaluation of the extracted filters
    cifar10 = load_cifar10()
    extractor = train_convolutional_kmeans(cifar10, n_centers=n_centers, 
            n_components=n_components, n_drop_components=n_drop_components, 
            save_images=True, n_EM_steps=n_EM_steps)

    # delete some big useless objects
    del extractor.patches_
    del extractor.patches_unpca_
    if save_extractor:
        cPickle.dump(extractor, open(save_extractor, 'wb'))

def features_from_saved_extractor(save_extractor='extractor.pkl', n_examples_to_use=50000,
        save_prefix="kmeans"):
    print 'Extracting convolutional k-means features from CIFAR-10'
    # This is in a function on its own because it can take a while (20 minutes).
    extractor = cPickle.load(open(save_extractor))

    X_train, y_train, X_test, y_test = load_cifar10()

    # for each image position, extract features from the entire dataset

    X_train_features = extractor.transform(X_train[:n_examples_to_use])
    X_test_features = extractor.transform(X_test[:n_examples_to_use])

    np.save('%s_X_train_features.npy'%save_prefix, X_train_features)
    np.save('%s_X_test_features.npy'%save_prefix, X_test_features)
    np.save('%s_y_train_labels.npy'%save_prefix, y_train[:n_examples_to_use])
    np.save('%s_y_test_labels.npy'%save_prefix, y_test[:n_examples_to_use])

def classify_features(n_examples_to_use=50000, alpha=.0001,n_iter=20,
        save_prefix="kmeans"):
    classif = SGDClassifier(
            loss='hinge',
            penalty='l2',
            alpha=alpha,
            shuffle=True,
            n_iter=n_iter,
            n_jobs=1)

    print classif

    X_train = np.load('%s_X_train_features.npy'%save_prefix)[:n_examples_to_use]
    X_test = np.load('%s_X_test_features.npy'%save_prefix)[:n_examples_to_use]
    y_train = np.load('%s_y_train_labels.npy'%save_prefix)[:n_examples_to_use]
    y_test = np.load('%s_y_test_labels.npy'%save_prefix)[:n_examples_to_use]

    print 'loaded data of shape', X_train.shape, y_train.shape
    print 'loaded data of shape', X_test.shape, y_test.shape
    print "scaling features to centered, unit variance vectors"

    scaler = Scaler().fit(X_train)
    X_train = scaler.transform(X_train, copy=False)
    X_test = scaler.transform(X_test,copy=False)

    print 'training svm'
    classif.fit(X_train.reshape((X_train.shape[0], -1)),y_train)


    pred_train = classif.predict(X_train.reshape((X_train.shape[0],-1)))
    pred_test = classif.predict(X_test.reshape((X_test.shape[0],-1)))

    print 'train accuracy', (pred_train == y_train).mean()
    print 'test accuracy', (pred_test == y_test).mean()


if __name__ == '__main__':
    # simle command-line base calling syntax:
    # example test_features using little data: 
    #   $ python <this file> classify_features  5000 .01 10
    cmd = sys.argv[1]
    args = [eval(arg) for arg in sys.argv[2:]]
    sys.exit(globals()[cmd](*args))

