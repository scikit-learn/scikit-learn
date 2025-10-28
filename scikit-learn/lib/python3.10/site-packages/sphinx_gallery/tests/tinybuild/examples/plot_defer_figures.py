"""
Test plot deferring
===================

This tests the ``sphinx_gallery_defer_figures`` flag.
"""

import matplotlib.pyplot as plt

# %%
# This code block should produce no plot.

plt.plot([0, 1])
plt.plot([1, 0])
# sphinx_gallery_defer_figures

# %%
# This code block should produce a plot with three lines.

plt.plot([2, 2])
plt.show()
