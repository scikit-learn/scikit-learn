"""
===========================================================
A demo of SelfOrganisingMap with colored neurons
===========================================================

XXX : add description of example.

"""
print __doc__
import pylab as pl
from matplotlib.colors import ListedColormap, NoNorm, rgb2hex
import numpy as np
from scikits.learn.cluster import SelfOrganizingMap


def plot(neurons):
    assert neurons.shape[-1] == 3
    h, w, d = neurons.shape
    hexmap = np.apply_along_axis(rgb2hex, 1, neurons.reshape(-1, 3) / 256)
    index = np.arange(h * w).reshape(h, w)
    pl.pcolor(index, cmap=ListedColormap(hexmap), norm=NoNorm())

train = np.array([[0, 0, 0],       # black
                  [255, 255, 255], # white
                  [255, 0, 0],     # red
                  [0, 255, 0],     # green
                  [0, 0, 255],     # blue
                  [255, 255, 0],   # yellow
                  [0, 255, 255],   # cyan
                  [255, 0, 255]])  # magenta

init = np.random.rand(16, 16, 3) * 255

pl.subplot(1, 2, 1, aspect='equal')
plot(init)
pl.title('Initial map')

som = SelfOrganizingMap(init, n_iterations=1024,
                        init='matrix', learning_rate=1)
som.fit(train)

pl.subplot(1, 2, 2, aspect='equal')
plot(som.neurons_)
pl.title('Organized Map')
F = pl.gcf()
F.set_size_inches((40, 20))
pl.show()
