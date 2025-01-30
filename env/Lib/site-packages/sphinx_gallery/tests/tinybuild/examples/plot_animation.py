"""
Animation support
=================

Show an animation, which should end up nicely embedded below.
"""

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

# Adapted from
# https://matplotlib.org/gallery/animation/basic_example.html


def _update_line(num):
    line.set_data(data[..., :num])
    return (line,)


fig_0, ax_0 = plt.subplots(figsize=(5, 1))

# sphinx_gallery_thumbnail_number = 2
fig_1, ax_1 = plt.subplots(figsize=(5, 5))
data = np.random.RandomState(0).rand(2, 25)
(line,) = ax_1.plot([], [], "r-")
ax_1.set(xlim=(0, 1), ylim=(0, 1))
ani = animation.FuncAnimation(fig_1, _update_line, 25, interval=100, blit=True)

fig_2, ax_2 = plt.subplots(figsize=(5, 5))
