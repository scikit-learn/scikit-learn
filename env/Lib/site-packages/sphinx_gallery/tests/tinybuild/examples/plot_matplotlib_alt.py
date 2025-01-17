"""
======================
Matplotlib alt text
======================

This example tests that the alt text is generated correctly for matplotlib
figures.
"""

import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 1, constrained_layout=True)
axs[0].plot([1, 2, 3])
axs[0].set_title("subplot 1")
axs[0].set_xlabel("x label")
axs[0].set_ylabel("y lab")
fig.suptitle("This is a\nsup title")

axs[1].plot([2, 3, 4])
axs[1].set_title("subplot 2")
axs[1].set_xlabel("x label")
axs[1].set_ylabel("y label")

plt.show()


# %%
# Several titles.

# sphinx_gallery_thumbnail_number = -1

plt.plot(range(10))

plt.title("Center Title")
plt.title("Left Title", loc="left")
plt.title("Right Title", loc="right")

plt.show()
