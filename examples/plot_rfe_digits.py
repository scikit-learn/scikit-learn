"""
=============================
Recursive feature elimination
=============================

A recursive feature elimination is performed prior to SVM classification.
"""
print __doc__

from scikits.learn.svm import SVC
from scikits.learn import datasets
from scikits.learn.feature_selection import RFE

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
rfe.fit(X, y)

image_ranking_ = rfe.ranking_.reshape(digits.images[0].shape)

import pylab as pl
pl.matshow(image_ranking_)
pl.colorbar()
pl.title('Ranking of voxels with RFE')
pl.show()

