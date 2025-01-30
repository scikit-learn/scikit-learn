"""
Command line arguments support
==============================

Use command line arguments to control example script.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == "plot":
    fig_0, ax_0 = plt.subplots(figsize=(5, 1))

    x = np.arange(0, 10.0, 1)
    ax_0.plot(x, x**2)
    plt.show()
