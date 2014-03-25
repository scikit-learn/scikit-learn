"""
===========================================================
A demo of SelfOrganisingMap with colored neurons
===========================================================

Example for SOM clustering using 3 dimensionals vectors (RGB)
with 8 colors (black, white, red, green, blue, yellow, cyan, magenta)

"""

print __doc__

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, NoNorm, rgb2hex
from sklearn.cluster import SelfOrganizingMap


def plot(neurons):
    assert neurons.shape[-1] == 3
    n_centers, dim = neurons.shape
    neurons[np.where(neurons > 256)] = 256
    neurons[np.where(neurons < 0)] = 0
    hexmap = np.apply_along_axis(rgb2hex, 1, neurons / 256)
    index = np.arange(n_centers).reshape(np.floor(np.sqrt(n_centers)),
                                         n_centers/np.floor(np.sqrt(n_centers)))
    plt.pcolor(index, cmap=ListedColormap(hexmap), norm=NoNorm())

train = np.array([[0, 0, 0],        # black
                  [255, 255, 255],  # white
                  [255, 0, 0],      # red
                  [0, 255, 0],      # green
                  [0, 0, 255],      # blue
                  [255, 255, 0],    # yellow
                  [0, 255, 255],    # cyan
                  [255, 0, 255]])   # magenta

init = np.random.rand(256, 3) * 255

plt.subplot(1, 2, 1, aspect='equal')
plot(init)
plt.title('Initial map')

som = SelfOrganizingMap(adjacency=(16, 16), n_iterations=1024,
                        init=init, learning_rate=1)
som.fit(train)

plt.subplot(1, 2, 2, aspect='equal')
plot(som.centers_)
plt.title('Organized Map')
F = plt.gcf()
F.set_size_inches((40, 20))
plt.show()
