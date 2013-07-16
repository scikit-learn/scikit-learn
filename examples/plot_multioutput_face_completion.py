"""
============================================
Face completion with multi-output estimators
============================================

This example shows the use of multi-output estimator to complete images.
The goal is to predict the lower half of a face given its upper half.

The first row of images shows true faces. The second half illustrates
how an estimator completes the lower half of those faces.

"""
print(__doc__)

import numpy as np
import pylab as pl

from sklearn.datasets import fetch_olivetti_faces
from sklearn.utils.validation import check_random_state

from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import RidgeCV

# Load the faces datasets
data = fetch_olivetti_faces()
targets = data.target

data = data.images.reshape((len(data.images), -1))
train = data[targets < 30]
test = data[targets >= 30]  # Test on independent people
n_pixels = data.shape[1]

X_train = train[:, :int(0.5 * n_pixels)]  # Upper half of the faces
y_train = train[:, int(0.5 * n_pixels):]  # Lower half of the faces
X_test = test[:, :int(0.5 * n_pixels)]
y_test = test[:, int(0.5 * n_pixels):]

n_faces = 5
rng = check_random_state(4)
face_ids = rng.randint(X_test.shape[0], size=(n_faces, ))
X_test = X_test[face_ids, :]
y_test = y_test[face_ids, :]

# Build a multi-output forest
ESTIMATORS = {
    "extra trees": ExtraTreesRegressor(n_estimators=10, max_features=32,
                                       random_state=0),
    "knn": KNeighborsRegressor(),
    "linear regression": LinearRegression(),
    "ridge": RidgeCV(),
}

y_test_predict = dict()
for name, estimator in ESTIMATORS.items():
    estimator.fit(X_train, y_train)
    y_test_predict[name] = estimator.predict(X_test)

# Plot the completed faces
image_shape = (64, 64)

n_cols = 1 + len(ESTIMATORS)
pl.figure(figsize=(2. * n_cols, 2.26 * n_faces))
pl.suptitle("Face completion with multi-output estimators", size=16)

for i in range(n_faces):
    true_face = np.hstack((X_test[i], y_test[i]))

    if i != 0:
        sub = pl.subplot(n_faces, n_cols, i * n_cols + 1)
    else:
        sub = pl.subplot(n_faces, n_cols, i * n_cols + 1,
                         title="true faces")
    sub.axis("off")
    sub.imshow(true_face.reshape(image_shape),
               cmap=pl.cm.gray,
               interpolation="nearest")

    for j, est in enumerate(sorted(ESTIMATORS)):
        completed_face = np.hstack((X_test[i], y_test_predict[est][i]))

        if i != 0:
            sub = pl.subplot(n_faces, n_cols, i * n_cols + 2 + j)
        else:
            sub = pl.subplot(n_faces, n_cols, i * n_cols + 2 + j,
                             title=est)
        sub.axis("off")
        sub.imshow(completed_face.reshape(image_shape),
                   cmap=pl.cm.gray,
                   interpolation="nearest")

pl.show()
