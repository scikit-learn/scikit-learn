"""
Astronomy Tutorial: exercise 2

Photometric redshift determination

usage: python exercise_02.py datadir

  - datadir is $TUTORIAL_DIR/data/sdss_photoz
    This directory should contain the files:
       - sdss_photoz.npy

Here we will take a closer look at the photometric redshift problem discussed
in section 5 of the tutorial.  Using the decision tree classifier, we'll take
a look at the 4-color observations of just over 400,000 points.

The point of this exercise is to answer the question: how can we get the rms
error down to below 0.1?  Would it be a better use of telescope time to
observe more objects, or to observe additional features of the objects
in the data set?  We'll use the techniques discussed in section 3 of the
tutorial.
"""
import os, sys
import numpy as np
import pylab as pl

from sklearn.tree import DecisionTreeRegressor
from sklearn import metrics

try:
    datadir = sys.argv[1]
except:
    print __doc__
    sys.exit()

def compute_rms_error(y_pred, y_true):
    """Compute the rms error between the arrays y_pred and y_true"""
    return np.sqrt(metrics.mean_squared_error(y_pred, y_true))

def compute_outlier_fraction(y_pred, y_true, cutoff=0.2):
    """Compute the outlier rate between the arrays y_pred and y_true"""
    return np.sum((abs(y_pred - y_true) > cutoff)) * 1. / len(y_pred)

#------------------------------------------------------------
# load data and compute colors
data = np.load(os.path.join(datadir, 'sdss_photoz.npy'))

# here we'll truncate the data to 50,000 points.  This will allow the code
# below to be run quickly while it's being written.  When you're satisfied
# that the code is ready to go, you can comment out this line.
data = data[:50000]

print '%i points' % data.shape[0]
u, g, r, i, z = [data[f] for f in 'ugriz']

X = np.zeros((len(data), 4))
X[:, 0] = u - g
X[:, 1] = g - r
X[:, 2] = r - i
X[:, 3] = i - z

y = data['redshift']

#------------------------------------------------------------
# divide into training, cross-validation, and test samples
Ntot = len(y)

Ncv = Ntot / 5
Ntest = Ntot / 5
Ntrain = Ntot - Ncv - Ntest

X_train = X[:Ntrain]
y_train = y[:Ntrain]

X_cv = X[Ntrain:Ntrain + Ncv]
y_cv = y[Ntrain:Ntrain + Ncv]

X_test = X[Ntrain + Ncv:]
y_test = y[Ntrain + Ncv:]

#------------------------------------------------------------
# plot the Decision Tree error as a function of max_depth parameter
#
#  This is the first main part of the exercise.  This is photometric
#  redshift determination using DecisionTreeRegressor.  Here you'll plot
#  the training error and cross-validation error as a function of the
#  meta-parameter 'max_depth'.
#
#  You will create three arrays: max_depth_array, train_error, and cv_error.
#  Use at least 10 different values of max_depth, and compute the training
#  and cross-validation error associated with each of them.
#
#  note that the error can be computed with the function compute_rms_error()

max_depth_array = []
train_error = []
cv_error = []

# TODO:  compute the arrays max_depth_array, train_error, and cv_error

pl.figure()
pl.plot(max_depth_array, cv_error, label='cross-val error')
pl.plot(max_depth_array, train_error, label='training error')

pl.legend()
pl.xlabel('max depth')
pl.ylabel('error')

# select the value of max_depth which led to the best results
max_depth = max_depth_array[np.argmin(cv_error)]
print "max_depth = %i" % max_depth

#------------------------------------------------------------
# plot the Decision Tree error as a function of number of samples
#
#  This is the second main part of the exercise.  Here you'll plot the
#  training error and cross-validation error as a function of the
#  number of training samples.
#
#  You will create three arrays: n_samples_array, train_error, and cv_error.
#  Use at least 40 different values of n_samples, and compute the training
#  and cross-validation error associated with each of them.
#
#  Make sure that when computing the training error for each number of
#  samples, you use the same samples that the model was trained on.

n_samples_array = []
train_error = []
cv_error = []

# TODO:  compute the arrays n_samples_array, train_error, and cv_error

pl.figure()
pl.plot(n_samples_array, cv_error, label='cross-val error')
pl.plot(n_samples_array, train_error, label='training error')

pl.legend()
pl.xlabel('number of samples')
pl.ylabel('error')

#----------------------------------------------------------------------
# Use the whole dataset:
#  If you have been running your code on only a part of the dataset,
#  now that you have it working, you can run it on the full dataset
#  (note: this will take a long time to execute!)  You can do this by
#  commenting out the line
#     data = data[:50000]
#  above.  How does this change the results?


#------------------------------------------------------------
# Catastrophic Outliers
#  Though the rms error is one useful measure of the performance of an
#  algorithm, astronomers are often more interested in reducing the
#  'catastrophic outlier' rate.  Catastrophic outliers are points which
#  are given redshifts very far from the true value.  For accuracy of
#  cosmological results, this is often more important than the overall
#  rms error.
#
#  Here, you can re-implement te above tasks, plotting the catastrophic
#  outlier rate as a function of the max_depth parameter, and as a function
#  of the number of training points.  This can be accomplished either by
#  copying and pasting the above code here, or by modifying the above code.
#
#  To compute the catastrophic error rate, you can use the function
#  compute_outlier_fraction()

# TODO:  repeat the above two plots using catastrophic error rate

#----------------------------------------------------------------------
# Analyze the results
#
#  Compare your results to the discussion of bias and variance in section
#  3.  How do you think these results could be improved?  Is it better to
#  spend telescope time increasing the size of the training set, or would
#  it be better to measure more features of the objects we already have?
#  Does this recommendation change if the astronomer is interested in
#  minimizing the number of catastrophic outliers rather than the rms error?

pl.show()
