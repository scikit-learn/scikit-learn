"""
Recursive feature elimination
=======================================================================

A recursive feature elimination is performed prior to SVM classification.
"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold, GridSearchCV
from scikits.learn import datasets
from scikits.learn.rfe import RFE, RFECV



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
rfe = RFE(estimator=svc, n_features=1, percentage=0.1)
clf = SVC(kernel="linear",C=1)
rfe.fit(X,y)

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



from scikits.learn.cross_val import StratifiedKFold
def loss_func(y1, y2):
    return np.sum(y1 != y2)
svc = SVC(kernel='linear')
rfecv = RFECV(estimator=svc, n_features=2,
                percentage=0.1,loss_func=loss_func)
rfecv.fit(X, y, cv=StratifiedKFold(y, 2))
print rfe.support_

pl.figure()
pl.plot(rfecv.cv_scores_)

image_support_ = rfecv.ranking_.reshape(digits.images[0].shape)
pl.figure()
import pylab as pl
pl.matshow(image_support_)
pl.colorbar()
pl.title('Ranking of voxels with RFE')



