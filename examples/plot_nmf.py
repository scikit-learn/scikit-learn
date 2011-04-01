"""
===================================================
NMF for faces feature extraction
===================================================

The dataset used in this example is a preprocessed excerpt of the
"Labeled Faces in the Wild", aka LFW_:

  http://vis-www.cs.umass.edu/lfw/lfw-funneled.tgz (233MB)

.. _LFW: http://vis-www.cs.umass.edu/lfw/


"""


print __doc__

from time import time
import logging
import numpy as np
import pylab as pl

from scikits.learn.pca import RandomizedPCA
from scikits.learn.nmf import NMF
from scikits.learn.datasets import fetch_lfw_people


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
                    
lfw_people = fetch_lfw_people(min_faces_per_person=70, resize=0.4)

# reshape the data using the traditional (n_samples, n_features) shape
faces = lfw_people.data
n_samples, h, w = faces.shape

X = faces.reshape((n_samples, h * w))
n_features = X.shape[1]

################################################################################
# Compute a PCA (eigenfaces) on the face dataset (treated as unlabeled
# dataset): unsupervised feature extraction / dimensionality reduction
n_components = 42

print "Extracting the top %d eigenfaces from %d faces" % (
    n_components, X.shape[0])
t0 = time()
pca = RandomizedPCA(n_components=n_components, whiten=True).fit(X)
print "done in %0.3fs" % (time() - t0)

eigenfaces = pca.components_.T.reshape((n_components, h, w))

#print "Projecting the data on the eigenfaces orthonormal basis"
#t0 = time()
#X_pca = pca.transform(X)
#print "done in %0.3fs" % (time() - t0)

################################################################################
# Compute the NMF on the same data

print "Extracting %d non-negative features from %d faces" % (
    n_components, X.shape[0])
t0 = time()
nmf = NMF(n_comp=n_components, initial='nndsvd', tol=1e-2,
          beta=30, eta=0.01, sparseness="components").fit(X)
print "done in %0.3fs" % (time() - t0)

nmfaces = nmf.components_.T.reshape((n_components, h, w))

################################################################################
# Plotting the results

n_row, n_col = 7, 6

f1 = pl.figure(figsize=(1.8 * n_col, 2.4 * n_row))
f1.text(.5, .95, 'Principal components', horizontalalignment='center') 
pl.subplots_adjust(bottom=0, left=.01, right=.99, top=.90, hspace=.35)
for i in range(n_row * n_col):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(eigenfaces[i].reshape((h, w)), cmap=pl.cm.gray)
    pl.xticks(())
    pl.yticks(())

f2 = pl.figure(figsize=(1.8 * n_col, 2.4 * n_row))
f2.text(.5, .95, 'Non-negative components', horizontalalignment='center') 
pl.subplots_adjust(bottom=0, left=.01, right=.99, top=.90, hspace=.35)
for i in range(n_row * n_col):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(nmfaces[i].reshape((h, w)), cmap=pl.cm.gray)
    pl.xticks(())
    pl.yticks(())  

pl.show()