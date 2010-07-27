"""
Recursive feature elimination
=======================================================================

A recursive feature elimination is performed prior to SVM classification.
"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold, GridSearchCV
from scikits.learn import datasets
from scikits.learn.rfe import RFE



################################################################################
# Loading the Digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

################################################################################
# Create the RFE object and compute a cross-validated score

svc = SVC(kernel="linear", C=1)
rfe = RFE(estimator=svc, n_features=10, percentage=0.1)
clf = SVC(kernel="linear",C=1)

y_pred = []
y_true = []
for train, test in StratifiedKFold(y, 2):
    Xtrain, ytrain, Xtest, ytest = X[train], y[train], X[test], y[test]

    ### Fit and predict rfe
    support_ = rfe.fit(X[train], y[train]).support_
    y_pred.append(clf.fit(Xtrain[:,support_],ytrain).predict(
          Xtest[:,support_]))
    ranking_ = rfe.ranking_


y_true = np.concatenate(y_true)
y_pred = np.concatenate(y_pred)
classif_rate = np.mean(y_pred == y_true) * 100
print "Classification rate: %f" % classif_rate


    
image_support_ = rfe.support_.reshape(digits.images[0].shape)
import pylab as pl
pl.matshow(image_support_)
pl.colorbar()
pl.title('Voxels selected with RFE')

image_support_ = rfe.ranking_.reshape(digits.images[0].shape)
import pylab as pl
pl.matshow(image_support_)
pl.colorbar()
pl.title('Ranking of voxels with RFE')


### Create the happy curve
classif_score = []
for i in np.unique(ranking_):
    clf.fit(Xtrain[:,ranking_<=i],ytrain)
    classif_score.append(np.mean(clf.predict(Xtest[:,ranking_<=i])
                                == ytest) * 100)

import pylab as pl
pl.plot(classif_score)









