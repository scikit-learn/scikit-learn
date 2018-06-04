from __future__ import absolute_import, division, print_function

import six

import matplotlib
from matplotlib.font_manager import FontProperties
from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
import os.path


@image_comparison(baseline_images=["truetype-conversion"],
                  extensions=["pdf"])
def test_truetype_conversion():
    fontname = os.path.join(os.path.dirname(__file__), 'mpltest.ttf')
    fontname = os.path.abspath(fontname)
    fontprop = FontProperties(fname=fontname, size=80)
    matplotlib.rcParams['pdf.fonttype'] = 3
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.text(0, 0, "ABCDE", fontproperties=fontprop)
    ax.set_xticks([])
    ax.set_yticks([])
