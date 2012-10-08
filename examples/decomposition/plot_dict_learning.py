"""
"""
print __doc__

from time import time

import pylab as pl
import numpy as np

from scipy.misc import lena

from sklearn import datasets
from sklearn.cluster import MiniBatchKMeans

faces = datasets.fetch_lfw_people()
data = faces.data

###############################################################################
# Learn the dictionary of images

print 'Learning the dictionary... '
t0 = time()
kmeans = MiniBatchKMeans(n_clusters=100, batch_size=200, verbose=1,
                         max_no_improvement=30)
kmeans.fit(data)
dt = time() - t0
print 'done in %.2fs.' % dt

pl.figure(figsize=(4.2, 4))
for i, label in enumerate(np.unique(kmeans.labels_)):
    comp = np.mean(data[label == kmeans.labels_], axis=0)
    pl.subplot(10, 10, i + 1)
    pl.imshow(comp.reshape(faces.images[0].shape), cmap=pl.cm.gray,
              interpolation='nearest')
    pl.xticks(())
    pl.yticks(())

pl.suptitle('Dictionary of faces\n' +
            'Train time %.1fs on %d patches' % (dt, len(data)),
            fontsize=16)
pl.subplots_adjust(0.08, 0.02, 0.92, 0.85, 0.08, 0.23)


