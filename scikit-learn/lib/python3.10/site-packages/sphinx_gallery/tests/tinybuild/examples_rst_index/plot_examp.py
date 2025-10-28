"""
Introductory example - Plotting sin
===================================

"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xlabel(r"$x$")
plt.ylabel(r"$\sin(x)$")
# To avoid matplotlib text output
plt.show()
