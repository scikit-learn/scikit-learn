from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison
from matplotlib.transforms import Bbox

from mpl_toolkits.axisartist.clip_path import clip_line_to_rect


@image_comparison(baseline_images=['clip_path'],
                  extensions=['png'], style='default')
def test_clip_path():
    x = np.array([-3, -2, -1, 0., 1, 2, 3, 2, 1, 0, -1, -2, -3, 5])
    y = np.arange(len(x))

    fig, ax = plt.subplots()
    ax.plot(x, y, lw=1)

    bbox = Bbox.from_extents(-2, 3, 2, 12.5)
    rect = plt.Rectangle(bbox.p0, bbox.width, bbox.height,
                         facecolor='none', edgecolor='k', ls='--')
    ax.add_patch(rect)

    clipped_lines, ticks = clip_line_to_rect(x, y, bbox)
    for lx, ly in clipped_lines:
        ax.plot(lx, ly, lw=1, color='C1')
        for px, py in zip(lx, ly):
            assert bbox.contains(px, py)

    ccc = iter(['C3o', 'C2x', 'C3o', 'C2x'])
    for ttt in ticks:
        cc = six.next(ccc)
        for (xx, yy), aa in ttt:
            ax.plot([xx], [yy], cc)
