"""
Setting the Matplotlib backend
==============================
"""

# %%
# The Matplotlib backend should start as `agg`

import matplotlib

print(f"Matplotlib backend is {matplotlib.get_backend()}")
assert matplotlib.get_backend() == "agg"

# %%
# Changing the Matplotlib backend to `svg` should be possible

matplotlib.use("svg")
print(f"Matplotlib backend is {matplotlib.get_backend()}")
assert matplotlib.get_backend() == "svg"

# %%
# In a new code block, the Matplotlib backend should continue to be `svg`

print(f"Matplotlib backend is {matplotlib.get_backend()}")
assert matplotlib.get_backend() == "svg"
