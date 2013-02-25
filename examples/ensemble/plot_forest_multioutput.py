"""
=========================================
Face completion with multi-output forests
=========================================

This example shows the use of multi-output forests to complete images.
The goal is to predict the lower half of a face given its upper half.

The first row of images shows true faces. The second half illustrates
how the forest completes the lower half of those faces.

"""
print(__doc__)

import numpy as np
import pylab as pl

from sklearn.datasets import fetch_olivetti_faces
from sklearn.ensemble import ExtraTreesRegressor


# Load the faces datasets
data = fetch_olivetti_faces()
targets = data.target

data = data.images.reshape((len(data.images), -1))
train = data[targets < 30]
test = data[targets >= 30]  # Test on independent people
n_pixels = data.shape[1]

X_train = train[:, :int(0.5 * n_pixels)]  # Upper half of the faces
Y_train = train[:, int(0.5 * n_pixels):]  # Lower half of the faces
X_test = test[:, :int(0.5 * n_pixels)]
Y_test = test[:, int(0.5 * n_pixels):]

# Build a multi-output forest
forest = ExtraTreesRegressor(n_estimators=10,
                             max_features=32,
                             random_state=0)

forest.fit(X_train, Y_train)
Y_test_predict = forest.predict(X_test)

# Plot the completed faces
n_faces = 5
image_shape = (64, 64)

pl.figure(figsize=(2. * n_faces, 2.26 * 2))
pl.suptitle("Face completion with multi-output forests", size=16)

for i in range(1, 1 + n_faces):
    face_id = np.random.randint(X_test.shape[0])

    true_face = np.hstack((X_test[face_id], Y_test[face_id]))
    completed_face = np.hstack((X_test[face_id], Y_test_predict[face_id]))

    pl.subplot(2, n_faces, i)
    pl.axis("off")
    pl.imshow(true_face.reshape(image_shape),
              cmap=pl.cm.gray,
              interpolation="nearest")

    pl.subplot(2, n_faces, n_faces + i)
    pl.axis("off")
    pl.imshow(completed_face.reshape(image_shape),
              cmap=pl.cm.gray,
              interpolation="nearest")

pl.show()
