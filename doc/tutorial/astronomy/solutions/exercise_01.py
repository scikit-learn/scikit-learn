"""
Astronomy Tutorial: exercise 1

Classification of photometric sources

usage: python exercise_01.py datadir

  - datadir is $TUTORIAL_DIR/data/sdss_colors
    This directory should contain the files:
       - sdssdr6_colors_class_train.npy
       - sdssdr6_colors_class.200000.npy

Description:
In the tutorial, we used a Naive Bayes Classifier to separate Quasars
And Stars.  In this exercise, we will extend this classification scheme
using Gaussian Mixture Models.

The Gaussian Naive Bayes method starts by fitting an N-dimensional gaussian
distribution to each class of data.  When a test point is evaluated, the
relative log-likelihood from each distribution is used to predict the most
likely value.  We're going to extend this by fitting a sum of gaussians to
each distribution.

There are several places in this file with code to be filled-in as part of
the exercise.  Each of these is labeled TODO below.
"""
import os, sys
import numpy as np
import pylab as pl
from sklearn.mixture import gmm
from sklearn import metrics

try:
    datadir = sys.argv[1]
except:
    print __doc__
    sys.exit()

#----------------------------------------------------------------------
# Load data files
train_data = np.load(os.path.join(datadir,
                                  'sdssdr6_colors_class_train.npy'))
test_data = np.load(os.path.join(datadir,
                                 'sdssdr6_colors_class.200000.npy'))

# set the number of training points: using all points leads to a very
# long running time.  We'll start with 10000 training points.  This
# can be increased if desired.
Ntrain = 10000
#Ntrain = len(train_data)

np.random.seed(0)
np.random.shuffle(train_data)
train_data = train_data[:Ntrain]

#----------------------------------------------------------------------
# Split training data into training and cross-validation sets
N_crossval = Ntrain / 5
train_data = train_data[:-N_crossval]
crossval_data = train_data[-N_crossval:]

#----------------------------------------------------------------------
# Set up data
#
X_train = np.zeros((train_data.size, 4), dtype=float)
X_train[:, 0] = train_data['u-g']
X_train[:, 1] = train_data['g-r']
X_train[:, 2] = train_data['r-i']
X_train[:, 3] = train_data['i-z']
y_train = (train_data['redshift'] > 0).astype(int)
Ntrain = len(y_train)

X_crossval = np.zeros((crossval_data.size, 4), dtype=float)
X_crossval[:, 0] = crossval_data['u-g']
X_crossval[:, 1] = crossval_data['g-r']
X_crossval[:, 2] = crossval_data['r-i']
X_crossval[:, 3] = crossval_data['i-z']
y_crossval = (crossval_data['redshift'] > 0).astype(int)
Ncrossval = len(y_crossval)

#======================================================================
# Recreating Gaussian Naive Bayes
#
#   Here we will use Gaussian Mixture Models to duplicate our Gaussian
#   Naive Bayes results from earlier.  You'll create two sklearn.gmm.GMM()
#   classifier instances, named `clf_0` and `clf_1`.  Each should be
#   initialized with a single component, and diagonal covariance.
#   (hint: look at the doc string for sklearn.gmm.GMM to see how to set
#   this up).  The results should be compared to Gaussian Naive Bayes
#   to check if they're correct.
#
#   Objects to create:
#    - clf_0 : trained on the portion of the training data with y == 0
#    - clf_1 : trained on the portion of the training data with y == 1

#{{{ compute clf_0, clf_1
clf_0 = gmm.GMM(1, 'diag')
i0 = (y_train == 0)
clf_0.fit(X_train[i0])

clf_1 = gmm.GMM(1, 'diag')
i1 = (y_train == 1)
clf_1.fit(X_train[i1])
#}}}

# next we must construct the prior.  The prior is the fraction of training
# points of each type.
# 
# variables to compute:
#  - prior0 : fraction of training points with y == 0
#  - prior1 : fraction of training points with y == 1

#{{{ compute prior0, prior1
num0 = i0.sum()
num1 = i1.sum()

prior0 = num0 / float(Ntrain)
prior1 = num1 / float(Ntrain)
#}}}

# Now we use the prior and the classifiation to compute the log-likelihoods
#  of the cross-validation points.  The log likelihood is given by
#
#    logL(x) = clf.score(x) + log(prior)
#
#  You can use the function np.log() to compute the logarithm of the prior.
#  variables to compute:
#    logL : array, shape = (2, Ncrossval)
#            logL[0] is the log-likelihood for y == 0
#            logL[1] is the log-likelihood for y == 1
logL = None

#{{{ compute logL
logL = np.zeros((2, Ncrossval))
logL[0] = clf_0.score(X_crossval) + np.log(prior0)
logL[1] = clf_1.score(X_crossval) + np.log(prior1)
#}}}

# the predicted value for each sample is the index with the largest
# log-likelihood.
y_pred = np.argmax(logL, 0)

# now we print the results.  We'll use the built-in classification
# report function in sklearn.metrics.  This computes the precision,
# recall, and f1-score for each class.

print "------------------------------------------------------------"
print "One-component Gaussian Mixture:"
print "  results for cross-validation set:"
print metrics.classification_report(y_crossval, y_pred,
                                    target_names=['stars', 'QSOs'])



#----------------------------------------------------------------------
#  Run Gaussian Naive Bayes to double-check that our results are correct.
#  Because of rounding errors, it will not be exact, but the results should
#  be very close.
from sklearn.naive_bayes import GaussianNB
gnb = GaussianNB()
gnb.fit(X_train, y_train)
y_pred = gnb.predict(X_crossval)

print "------------------------------------------------------------"
print "Gaussian Naive Bayes"
print "  results for cross-validation set:"
print "  (results should be within ~0.01 of above results)"
print metrics.classification_report(y_crossval, y_pred,
                                    target_names=['stars', 'QSOs'])

#======================================================================
#  Parameter optimization:
#
#   Now take some time to experiment with the covariance type and the
#   number of components, to see if you can optimize the F1 score
#
#   Note that for a large number of components, the fit can take a long
#   time, and will be dependent on the starting position.  Use the
#   documentation string of GMM to determine the options for covariance.
#
#   It may be helpful to use only a subset of the training data while
#   experimenting with these parameter values.  This is called
#   "Meta-parameter optimization".  It can be accomplished automatically,
#   but here we are doing it by hand for learning purposes.
y_pred = None

#{{{ compute y_pred for cross-validation data
clf_0 = gmm.GMM(5, 'full', random_state=0)
i0 = (y_train == 0)
clf_0.fit(X_train[i0])

clf_1 = gmm.GMM(5, 'full', random_state=0)
i1 = (y_train == 1)
clf_1.fit(X_train[i1])

logL = np.zeros((2, Ncrossval))
logL[0] = clf_0.score(X_crossval) + np.log(prior0)
logL[1] = clf_1.score(X_crossval) + np.log(prior1)

y_pred = np.argmax(logL, 0)
#}}}

print "------------------------------------------------------------"
print "GMM with tweaked parameters:"
print "  results for cross-validation set"
print metrics.classification_report(y_crossval, y_pred,
                                    target_names=['stars', 'QSOs'])

#----------------------------------------------------------------------
# Test Data
# once you have maximized the cross-validation, you can apply the estimator
# to your test data, and check how it compares to the predicted results
# from the researcher who compiled it.

X_test = np.zeros((test_data.size, 4), dtype=float)
X_test[:, 0] = test_data['u-g']
X_test[:, 1] = test_data['g-r']
X_test[:, 2] = test_data['r-i']
X_test[:, 3] = test_data['i-z']
y_pred_literature = (test_data['label'] == 0).astype(int)
Ntest = len(y_pred_literature)

# here you should compute y_pred for the test data, using the classifiers
# clf_0 and clf_1 which you already trained above.

y_pred = None

#{{{ compute y_pred for test data
logL = np.zeros((2, Ntest))
logL[0] = clf_0.score(X_test) + np.log(prior0)
logL[1] = clf_1.score(X_test) + np.log(prior1)
y_pred = np.argmax(logL, 0)
#}}}

print "------------------------------------------------------------"
print "Comparison of current results with published results"
print "  results for test set"
print "    (treating published results as the 'true' result)"
print metrics.classification_report(y_pred_literature, y_pred,
                                    target_names=['stars', 'QSOs'])
