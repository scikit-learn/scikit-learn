"""
Force images to be displayed on separate lines per block
========================================================
This example demonstrates how the ``sphinx_gallery_multi_image_block`` option can be
used to force images to be displayed on separate lines for a specific block, rather than
the default behaviour of displaying them side-by-side.
"""

import matplotlib.pyplot as plt
import numpy as np

# %%

# Default behaviour
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.pcolormesh(np.random.randn(100, 100))

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.pcolormesh(np.random.randn(100, 100))


# %%

# sphinx_gallery_multi_image_block = "single"

# Non-default behaviour for just this cell
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.pcolormesh(np.random.randn(100, 100))

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.pcolormesh(np.random.randn(100, 100))
