import sys
import numpy

from scikits.learn.preprocessing import Scaler
from scikits.learn import svm
from scikits.learn.linear_model import SGDClassifier

if 0:
    n_examples_to_use = 10000
    classif = svm.SVC(
            kernel='linear',
            C=float(sys.argv[1]), #.03 was good
            cache_size=2000,
            )
else:
    n_examples_to_use = 50000
    classif = SGDClassifier(
            loss='hinge',
            penalty='l2',
            alpha=float(sys.argv[1]),
            shuffle=True,
            n_iter=20,
            n_jobs=1)

print classif

X_train = numpy.load('X_train_features.npy')[:n_examples_to_use]
X_test = numpy.load('X_test_features.npy')[:n_examples_to_use]
y_train = numpy.load('y_train_labels.npy')[:n_examples_to_use]
y_test = numpy.load('y_test_labels.npy')[:n_examples_to_use]

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

