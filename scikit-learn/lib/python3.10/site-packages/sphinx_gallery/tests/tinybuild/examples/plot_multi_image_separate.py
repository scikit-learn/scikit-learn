"""
Force images to be displayed on separate lines
==============================================
This example demonstrates how the ``sphinx_gallery_multi_image`` option can be used to
force images to be displayed on separate lines, rather than the default behaviour of
displaying them side-by-side.
"""

# sphinx_gallery_multi_image = "single"

# %%

import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.pcolormesh(np.random.randn(100, 100))

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.pcolormesh(np.random.randn(100, 100))
