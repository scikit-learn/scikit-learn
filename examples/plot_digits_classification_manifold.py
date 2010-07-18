"""
=============================================================
Recognizing hand-written digits with dimensionality reduction
=============================================================

An example showing how the scikit-learn can be used to recognize images of 
hand-written digits.

"""
# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Author: Matthieu Brucher <matthieu dot brucher at gmail dot com>
# License: Simplified BSD

# Standard scientific Python imports
import pylab as pl
import numpy as np

# The digits dataset
from scikits.learn import datasets
digits = datasets.load_digits()

fig1 = pl.figure()

# The data that we are interesting in is made of 8x8 images of digits,
# let's have a look at the first 3 images. We know which digit they
# represent: it is given in the 'target' of the dataset.
for index, (image, label) in enumerate(zip(digits.images, digits.target)[:4]):
    axes = fig1.add_subplot(2, 4, index+1)
    axes.imshow(image, cmap=pl.cm.gray_r)
    axes.set_title('Training: %i' % label)

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
data = digits.images.reshape((n_samples, -1))

from scikits.learn.manifold import Isomap

isomap = Isomap(embedded_opts={'nb_coords' : 5, 'neighbors' : 10},
    mapping_opts={'neighbors' : 9, 'tol' : 1e-9})
isomap.fit(data[:n_samples/2])

colors = np.array([(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1),
                   (0,0,0), (1,1,1), (.5, 0, 0), (0, .5, 0), (0, 0, .5)])

from mpl_toolkits.mplot3d import Axes3D

fig2 = pl.figure()
ax2 = Axes3D(fig2)
ax2.scatter(isomap.embedding_[:,0], isomap.embedding_[:,1],
    isomap.embedding_[:,2], c=colors[digits.target[:n_samples/2]])
fig2.suptitle('Training samples in embedded space')

reduced_test_data = isomap.predict(data[n_samples/2:])
fig3 = pl.figure()
ax3 = Axes3D(fig3)
ax3.scatter(reduced_test_data[:,0], reduced_test_data[:,1],
    reduced_test_data[:,2], c=colors[digits.target[n_samples/2:]])
fig3.suptitle('Predicted samples in embedded space')

# Import a classifier:
from scikits.learn import svm
classifier = svm.SVC()

# We learn the digits on the first half of the digits
classifier.fit(isomap.embedding_, digits.target[:n_samples/2])

# Now predict the value of the digit on the second half:
predicted = classifier.predict(reduced_test_data)

for index, (image, prediction) in enumerate(zip(
                                       digits.images[n_samples/2:], 
                                       predicted
                                    )[:4]):
    axes = fig1.add_subplot(2, 4, index+5)
    axes.imshow(image, cmap=pl.cm.gray_r)
    axes.set_title('Prediction: %i' % prediction)

from scikits.learn.metrics import confusion_matrix
cm = confusion_matrix(digits.target[n_samples/2:], predicted)

print "Correctly predicted samples: ", np.trace(cm), "/", np.sum(cm)

pl.matshow(cm)
pl.title('Confusion matrix')
pl.colorbar()

pl.show()
