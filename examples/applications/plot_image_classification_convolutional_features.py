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

import logging
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)

import numpy as np

from scikits.learn.grid_search import GridSearchCV
from scikits.learn.metrics import classification_report
from scikits.learn.metrics import confusion_matrix
from scikits.learn.feature_extraction.image import ConvolutionalKMeansEncoder
from scikits.learn.svm import SVC
from scikits.learn.svm import LinearSVC
from scikits.learn.preprocessing import Scaler
from scikits.learn.cluster.k_means_ import k_init
from scikits.learn.linear_model import SGDClassifier

def load_cifar10(keep_color_shift=True):
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

    if not keep_color_shift:
        # remove average rgb value from each image
        X_train -= X_train.mean(axis=2).mean(axis=1).reshape((X_train.shape[0], 1, 1, 3))
        X_test -= X_test.mean(axis=2).mean(axis=1).reshape((X_test.shape[0], 1, 1, 3))
    print 'loaded images as type', X_train.dtype

    if 0:
        # not scaling dataset because I'm afraid it might mess up
        # the local contrast normalization
        print "scaling images to centered, unit variance vectors"
        scaler = Scaler().fit(X_train)
        X_train = scaler.transform(X_train,copy=False)
        X_test = scaler.transform(X_test,copy=False)

    return X_train, y_train, X_test, y_test

def train_convolutional_kmeans(cifar10, 
        n_centers, 
        n_components, 
        n_drop_components=0,
        save_images=True,
        n_images_to_train_from=10000,
        n_EM_steps=30,
        center_mode='all'):
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
        verbose=1,
        center_mode=center_mode,
        )


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
        extractor.tile_pca_components(scale_each=True).save('pca.png')

    return extractor


################################################################################
# DRIVERS MEANT TO BE CALLED FROM COMMAND LINE (SEE __name__=='__main__' BELOW)
def debug_end_to_end():
    cifar10 = load_cifar10()
    extractor = train_convolutional_kmeans(cifar10, n_centers=20, 
            n_components=10,
            n_drop_components=2, 
            save_images=False,
            n_EM_steps=14,
            n_images_to_train_from=100)

    n_examples_to_use=10
    X_train, y_train, X_test, y_test = cifar10
    X_train_features = extractor.transform(X_train[:n_examples_to_use])
    X_test_features = extractor.transform(X_test[:n_examples_to_use])

#
def train_kmeans(save_extractor='extractor.pkl', 
        n_centers=400,
        n_components=80,
        n_drop_components=2,
        n_EM_steps=80,
        center_mode='all',
        keep_color_shift=False):
    print 'Training convolutional k-means'
    # Qualitative evaluation of the extracted filters
    cifar10 = load_cifar10(keep_color_shift)
    extractor = train_convolutional_kmeans(cifar10, n_centers=n_centers, 
            n_components=n_components, n_drop_components=n_drop_components, 
            save_images=True, n_EM_steps=n_EM_steps, center_mode=center_mode)

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

def classify_features(n_examples_to_use=50000, alpha=.0001, n_iter=20,
        save_prefix="kmeans", scale_input=True):

    #
    # I don't know what to say --
    # I can't figure out how to get good classification results
    # using scikits linear classifiers(!?)
    #

    shuffle = False
    normalize = False
    
    if 0:
        classif = SGDClassifier(
                loss='hinge',
                penalty='l2',
                alpha=1e-5,
                n_iter=50,
                n_jobs=1)
        shuffle = True
        normalize = True
    elif 1:
        classif = LinearSVC(C=0.01)
        
    else:
        classif = TheanoSGDClassifier(n_classes=10,
                learnrate=.005,
                l2_regularization=.1,
                center_and_normalize=True,
                anneal_epoch=50,
                n_epochs=500,
                validset_max_examples=0, #400 filters sees no overfitting
                min_feature_std=0.001)
        scale_input=True # the classifier does it

    print classif

    X_train = np.load('%s_X_train_features.npy'%save_prefix)[:n_examples_to_use]
    X_test = np.load('%s_X_test_features.npy'%save_prefix)[:n_examples_to_use]
    y_train = np.load('%s_y_train_labels.npy'%save_prefix)[:n_examples_to_use]
    y_test = np.load('%s_y_test_labels.npy'%save_prefix)[:n_examples_to_use]

    print 'loaded data of shape', X_train.shape, y_train.shape, X_train.dtype, y_train.dtype
    print 'loaded data of shape', X_test.shape, y_test.shape, X_test.dtype, y_test.dtype
    
    X_train = X_train.reshape((X_train.shape[0],-1))
    X_test = X_test.reshape((X_test.shape[0],-1))

    if scale_input:
        print "scaling features to centered, unit variance vectors"
        scaler = Scaler().fit(X_train)
        X_train = scaler.transform(X_train, copy=False)
        X_test = scaler.transform(X_test,copy=False)

    if shuffle:
        print "shuffling training data"
        idx = np.arange(len(y_train))
        np.random.seed(13)
        np.random.shuffle(idx)
        X_train = X_train[idx]
        y_train = y_train[idx]

    if normalize:
        print "normalize to avg norm 1"
        avg_norm = np.mean(np.sqrt(np.sum(X_train ** 2, axis=1))[:,np.newaxis])
        X_train /= avg_norm
        X_test /= avg_norm
        print "new avg norm:", np.mean(np.sqrt(np.sum(X_train ** 2, axis=1))[:,np.newaxis])

    print 'training classifier'
    classif.fit(X_train,y_train)


    pred_train = classif.predict(X_train)
    pred_test = classif.predict(X_test)

    print 'train accuracy', (pred_train == y_train).mean()
    print 'test accuracy', (pred_test == y_test).mean()


if 0:
    # Example of a BaseEstimator implemented with Theano.
    # This could use the GPU, except that 
    # a) linear regression isn't really worth it, and
    # b) the multi_hinge_margin Op is only implemented for the CPU.
    #
    import numpy as np
    try:
        import theano
        from theano import tensor
    except ImportError:
        print('Failed to import Theano - see installation instructions at '
                'http://www.deeplearning.net/software/theano/') 
        raise
    try:
        from pylearn.shared.layers.logreg import LogisticRegression
        import pylearn.gd.sgd
        from pylearn.formulas.costs import multi_hinge_margin
    except ImportError:
        print('Failed to import pylearn - clone it from https://hg.assembla.com/pylearn')
        raise

    class TheanoSGDClassifier(object):
        def __init__(self,
                n_classes,
                batchsize=100, 
                learnrate = 0.005,
                l1_regularization = 0.0,
                l2_regularization = 0.0,
                min_feature_std =0.3,
                n_epochs = 100,
                anneal_epoch=20,
                center_and_normalize=False,
                validset_fraction=.2,
                validset_max_examples=5000,
                copy_X=True,
                loss_fn='hinge',
                ):
            # add arguments to class
            self.__dict__.update(locals()); del self.self

        def fit(self, X, y):
            batchsize = self.batchsize

            n_valid = int(min(self.validset_max_examples, self.validset_fraction * X.shape[0]))
            # increase to a multiple of batchsize
            while n_valid % batchsize:
                n_valid += 1

            n_train = X.shape[0] - n_valid

            # decrease to a multiple of batchsize
            while n_train % batchsize:
                n_train -= 1

            if self.center_and_normalize and self.copy_X:
                X = X.copy()

            train_features = X[:n_train]
            valid_features = X[n_train:]
            train_labels = y[:n_train]
            valid_labels = y[n_train:]

            if self.center_and_normalize:
                print("Computing mean and std.dev")

                #this loop seems more memory efficient than numpy
                m= np.zeros(train_features.shape[1]) 
                msq= np.zeros(train_features.shape[1])
                for i in xrange(train_features.shape[0]):
                    alpha = 1.0 / (i+1)
                    v = train_features[i]
                    m = alpha * v + (1-alpha)*m
                    msq = alpha * v*v + (1-alpha)*msq

                self.X_mean_ = theano.shared(m.astype(X.dtype))
                self.X_std_ = theano.shared(
                        np.maximum(
                            self.min_feature_std, 
                            np.sqrt(msq - m*m)).astype(X.dtype))

                X -= self.X_mean_.get_value()
                X /= self.X_std_.get_value()

            x_i = tensor.matrix(dtype=X.dtype)
            y_i = tensor.vector(dtype=y.dtype)
            lr = tensor.scalar(dtype=X.dtype)

            feature_logreg = LogisticRegression.new(x_i,
                    n_in = train_features.shape[1], n_out=self.n_classes,
                    dtype=x_i.dtype)

            if self.loss_fn=='log':
                traincost = feature_logreg.nll(y_i).sum()
            elif self.loss_fn=='hinge':
                raw_output = tensor.dot(feature_logreg.input, feature_logreg.w)+feature_logreg.b
                traincost = multi_hinge_margin(raw_output, y_i).sum()
            else:
                raise NotImplementedError(self.loss_fn)
            traincost = traincost + abs(feature_logreg.w).sum() * self.l1_regularization
            traincost = traincost + (feature_logreg.w**2).sum() * self.l2_regularization
            train_logreg_fn = theano.function([x_i, y_i, lr], 
                    [feature_logreg.nll(y_i).mean(),
                        feature_logreg.errors(y_i).mean()],
                    updates=pylearn.gd.sgd.sgd_updates(
                        params=feature_logreg.params,
                        grads=tensor.grad(traincost, feature_logreg.params),
                        stepsizes=[lr/batchsize,lr/(10*batchsize)]))

            test_logreg_fn = theano.function([x_i, y_i],
                    feature_logreg.errors(y_i))

            if self.center_and_normalize:
                feature_logreg_test = LogisticRegression(
                        (x_i - self.X_mean_)/self.X_std_,
                        feature_logreg.w,
                        feature_logreg.b)
                self.predict_fn_ = theano.function([x_i], feature_logreg_test.argmax)
            else:
                self.predict_fn_ = theano.function([x_i], feature_logreg.argmax)

            best_epoch = -1
            best_epoch_valid = -1
            best_epoch_train = -1
            best_epoch_test = -1
            valid_rate=-1
            test_rate=-1
            train_rate=-1

            for epoch in xrange(self.n_epochs):
                # validate
                # Marc'Aurelio, you crazy!!
                # the division by batchsize is done in the cost function
                e_lr = np.float32(self.learnrate / max(1.0, np.floor(max(1.,
                    (epoch+1)/float(self.anneal_epoch))-2)))

                if n_valid:
                    l01s = []
                    for i in xrange(n_valid/batchsize):
                        x_i = valid_features[i*batchsize:(i+1)*batchsize]
                        y_i = valid_labels[i*batchsize:(i+1)*batchsize]

                        #lr=0.0 -> no learning, safe for validation set
                        l01 = test_logreg_fn((x_i), y_i)
                        l01s.append(l01)
                    valid_rate = 1-np.mean(l01s)
                    #print('Epoch %i validation accuracy: %f'%(epoch, valid_rate))

                    if valid_rate > best_epoch_valid:
                        best_epoch = epoch
                        best_epoch_test = test_rate
                        best_epoch_valid = valid_rate
                        best_epoch_train = train_rate

                    print('Epoch=%i best epoch %i valid %f test %f best train %f current train %f'%(
                        epoch, best_epoch, best_epoch_valid, best_epoch_test, best_epoch_train, train_rate))
                    if epoch > self.anneal_epoch and epoch > 2*best_epoch:
                        break
                else:
                    print('Epoch=%i current train %f'%( epoch, train_rate))

                #train
                l01s = []
                nlls = []
                for i in xrange(n_train/batchsize):
                    x_i = train_features[i*batchsize:(i+1)*batchsize]
                    y_i = train_labels[i*batchsize:(i+1)*batchsize]
                    nll, l01 = train_logreg_fn((x_i), y_i, e_lr)
                    nlls.append(nll)
                    l01s.append(l01)
                train_rate = 1-np.mean(l01s)
                #print('Epoch %i train accuracy: %f'%(epoch, train_rate))

        def predict(self, X):
            return self.predict_fn_(X)


if __name__ == '__main__':
    # simle command-line base calling syntax:
    # example test_features using little data: 
    #   $ python <this file> classify_features  5000 .01 10
    cmd = sys.argv[1]
    args = [eval(arg) for arg in sys.argv[2:]]
    sys.exit(globals()[cmd](*args))

