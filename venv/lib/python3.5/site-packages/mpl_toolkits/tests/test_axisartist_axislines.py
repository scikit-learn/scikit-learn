from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison

from mpl_toolkits.axisartist.axislines import SubplotZero, Subplot


@image_comparison(baseline_images=['SubplotZero'],
                  extensions=['png'], style='default')
def test_SubplotZero():
    fig = plt.figure()

    ax = SubplotZero(fig, 1, 1, 1)
    fig.add_subplot(ax)

    ax.axis["xzero"].set_visible(True)
    ax.axis["xzero"].label.set_text("Axis Zero")

    for n in ["top", "right"]:
        ax.axis[n].set_visible(False)

    xx = np.arange(0, 2 * np.pi, 0.01)
    ax.plot(xx, np.sin(xx))
    ax.set_ylabel("Test")


@image_comparison(baseline_images=['Subplot'],
                  extensions=['png'], style='default')
def test_Subplot():
    fig = plt.figure()

    ax = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax)

    xx = np.arange(0, 2 * np.pi, 0.01)
    ax.plot(xx, np.sin(xx))
    ax.set_ylabel("Test")

    ax.axis["top"].major_ticks.set_tick_out(True)
    ax.axis["bottom"].major_ticks.set_tick_out(True)

    ax.axis["bottom"].set_label("Tk0")
