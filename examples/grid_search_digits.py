"""Run parameter estimation using grid search
in a nested cross-validation setting.
"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold, KFold
from scikits.learn import datasets
from scikits.learn.grid_search import GridSearch

# The Digits dataset
digits = datasets.load_digits()

# # The data that we are interesting in is made of 8x8 images of digits,
# # let's have a look at the first 3 images. We know which digit they
# # represent: it is given in the 'target' of the dataset.
# for index, (image, label) in enumerate(zip(digits.images, digits.target)[:4]):
#     pl.subplot(2, 4, index+1)
#     pl.imshow(image, cmap=pl.cm.gray_r)
#     pl.title('Training: %i' % label)

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

parameters = {'kernel':('linear', 'rbf'), 'C':[0.1, 1]}

def loss_func(y1, y2):
    return np.mean(y1 != y2)

def cv(n_samples):
    return KFold(n_samples, 2)

clf = GridSearch(SVC, parameters, cv, loss_func, n_jobs=2)

"""
Run crossvalidation (different than the nested crossvalidation that is used
to select the best classifier). The classifier is optimized by "nested"
crossvalidation
"""
n_samples, n_features = X.shape
y_pred = np.zeros_like(y)
for train, test in StratifiedKFold(y, 2):
    y_pred[test] = clf.fit(X[train], y[train]).predict(X[test]).astype(np.int)

classif_rate = np.mean(y_pred == y) * 100
print "Classification rate : %f" % classif_rate
