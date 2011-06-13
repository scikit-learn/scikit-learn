"""
========================================
Sparse PCA for digits feature extraction
========================================

:ref:`SparsePCA` extracting the first 12 components out of the subset of
the digits dataset containing only the digit 3.

"""


print __doc__

import numpy as np
import matplotlib.pyplot as plt

from scikits.learn.decomposition import SparsePCA
from scikits.learn.datasets import load_digits

rows, cols = 4, 3
digits = load_digits()
threes = digits.data[digits.target == 3]
threes -= threes.mean(axis=0)  # todo: use preprocessors
# low tolerance for high speed
model = SparsePCA(n_components=rows * cols, tol=1e-3)
model.fit(threes)
span = np.max(np.abs(model.components_))

fig = plt.figure(figsize=(1.2 * cols, 1.4 * rows))
for i, comp in enumerate(model.components_):
    plt.subplot(rows, cols, i + 1)
    plt.imshow(np.reshape(comp, (8, 8)), interpolation='nearest',
               vmin=-span, vmax=span, cmap=plt.cm.PuOr)
    plt.xticks(())
    plt.yticks(())
plt.subplots_adjust(0.01, 0.15, 0.99, 0.99, 0.04, 0.)
cax = fig.add_axes([0.1, 0.06, 0.8, 0.04])
plt.colorbar(cax=cax, orientation='horizontal')
plt.show()
