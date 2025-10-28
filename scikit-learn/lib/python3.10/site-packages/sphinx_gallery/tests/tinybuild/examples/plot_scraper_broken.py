"""
Error during scraping
=====================

The error is actually introduced by the resetter in ``conf.py``. It mocks
a "zero-size reduction" error in ``fig.savefig``.
"""

import matplotlib.pyplot as plt

fig = plt.figure()
