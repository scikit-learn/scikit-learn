"""
==================
"SVG":-`graphics_`
==================

Make sure we can embed SVG graphics.
Use title that has punctuation marks.
"""

import matplotlib.pyplot as plt
import numpy as np
from local_module import N  # N = 1000

t = np.arange(N) / float(N)
win = np.hanning(N)
plt.figure()
plt.plot(t, win, color="r")
plt.text(0, 1, "svg", size=40, va="top")
