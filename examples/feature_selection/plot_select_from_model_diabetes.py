"""
===================================================
Feature selection using SelectFromModel and LassoCV
===================================================

Use SelectFromModel meta-transformer along with Lasso to select the best
couple of features from the diabetes dataset.

Diabetes dataset consists of 10 variables (features) collected from 442
diabetes patients. This example shows how to use SelectFromModel and LassoCv to
find the best two features predictiong disease progression after one year from
the baseline.

Authors: Manoj Kumar <mks542@nyu.edu>
         Maria Telenczuk
License: BSD 3 clause
"""
print(__doc__)

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import load_diabetes
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LassoCV

##############################################################################
# Load the data
# ---------------------------------------------------------
#
# First, let's load the diabetes dataset which is available from within
# sklearn. Then, we will look the features used of the diabates patients:

diabetes = load_diabetes()
X = diabetes.data
y = diabetes.target

feature_names = data.feature_names
print(features_names)

##############################################################################
# Find importance of the features
# ---------------------------------------------------------
#
# Since the L1 norm promotes sparsity of features only a subset of the provided
# variables should be selected. To decide which of the features are the most
# important we are going to use LassoCV estimator:

clf = LassoCV().fit(X, y)
importance = clf.coef_
print(importance)


idx_important = importance.argsort()[-2:][::-1]


clf = LassoCV()
sfm = SelectFromModel(clf, threshold='mean')
sfm.fit(X, y)
importance = sfm.estimator_.coef_
print(importance)

# The numbers with the highest values show the hightest importance of the
# feature 
# can choose the minimum threshold directly (here 0.25). 
# Features will be considered unimportant and will be removed if their features
# coef_ values are below this threshold.
# the transform() function will remove features to those considered important

clf = LassoCV()
sfm = SelectFromModel(clf, threshold=0.25)
sfm.fit(X, y)
n_features = sfm.transform(X).shape[1]
import pdb; pdb.set_trace()

##############################################################################
# 
# ---------------------------------------------------------
#
# Reset the threshold till the number of features equals two.
# Note that the attribute can be set directly instead of repeatedly
# fitting the metatransformer.

while n_features > 2:
    sfm.threshold += 0.1
    X_transform = sfm.transform(X)
    n_features = X_transform.shape[1]

##############################################################################
# 
# ---------------------------------------------------------
#
# Plot the selected two features from X.

plt.title(
    "Features from diabets using SelectFromModel with "
    "threshold %0.3f." % sfm.threshold)
feature1 = X_transform[:, 0]
feature2 = X_transform[:, 1]
plt.plot(feature1, feature2, 'r.')
plt.xlabel("Feature number 1")
plt.ylabel("Feature number 2")
plt.ylim([np.min(feature2), np.max(feature2)])
plt.show()
