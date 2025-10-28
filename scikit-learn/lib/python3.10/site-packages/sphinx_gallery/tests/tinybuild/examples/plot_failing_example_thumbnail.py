"""
Failing example test with normal thumbnail behaviour
====================================================
Test files in expected_failing_examples run, but normal thumbnail behaviour is
retained when sphinx_gallery_failing_thumbnail = False.
"""

# sphinx_gallery_failing_thumbnail = False

import matplotlib.pyplot as plt
import numpy as np

plt.pcolormesh(np.random.randn(100, 100))

# %%
# will raise AssertionError

assert False
