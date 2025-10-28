"""
Failing example test
====================
Test failing thumbnail appears for files in expected_failing_examples.
"""

import matplotlib.pyplot as plt
import numpy as np

plt.pcolormesh(np.random.randn(100, 100))

# %%
# will raise AssertionError

assert False
