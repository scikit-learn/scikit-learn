"""
===================================================
Feature selection using SelectFromModel and LassoCV
===================================================

Use SelectFromModel meta-transformer along with Lasso to select the best
couple of features from the diabetes dataset.
"""
# Authors: Manoj Kumar <mks542@nyu.edu>
#          Maria Telenczuk
# License: BSD 3 clause

print(__doc__)

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import load_diabetes
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LassoCV

# Load the diabetes dataset.
X, y = load_diabetes(return_X_y=True)

# We use the base estimator LassoCV since the L1 norm
# promotes sparsity of features.
clf = LassoCV()

# Set a minimum threshold of 0.25
sfm = SelectFromModel(clf, threshold=0.25)
sfm.fit(X, y)
n_features = sfm.transform(X).shape[1]

# Reset the threshold till the number of features equals two.
# Note that the attribute can be set directly instead of repeatedly
# fitting the metatransformer.
while n_featuresc > 2:
    sfmc.threshold += 0.1
    X_transform = sfmc.transform(Xc)
    n_featuresc = X_transform.shape[1]

# Plot the selected two features from X.
plt.title(
    "Features selected from diabets using SelectFromModel with "
    "threshold %0.3f." % sfm.threshold)
feature1 = X_transform[:, 0]
feature2 = X_transform[:, 1]
plt.plot(feature1, feature2, 'r.')
plt.xlabel("Feature number 1")
plt.ylabel("Feature number 2")
plt.ylim([np.min(feature2), np.max(feature2)])
plt.show()
