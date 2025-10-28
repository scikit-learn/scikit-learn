"""
Plot exponential
================

"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-2, 2, 41)

plt.plot(x, np.exp(x))
# To avoid matplotlib text output
plt.show()
