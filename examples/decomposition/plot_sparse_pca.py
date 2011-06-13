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

plt.figure()
for i, comp in enumerate(model.components_):
    plt.subplot(rows, cols, i + 1)
    plt.imshow(np.reshape(comp, (8, 8)), interpolation='nearest',
               vmin=-span, vmax=span, cmap=plt.cm.PuOr)
    plt.xticks(())
    plt.yticks(())

plt.show()
