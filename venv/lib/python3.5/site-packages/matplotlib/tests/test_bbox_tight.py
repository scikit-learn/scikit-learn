from __future__ import absolute_import, division, print_function

import numpy as np

from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter


@image_comparison(baseline_images=['bbox_inches_tight'], remove_text=True,
                  savefig_kwarg=dict(bbox_inches='tight'))
def test_bbox_inches_tight():
    #: Test that a figure saved using bbox_inches='tight' is clipped correctly
    data = [[66386, 174296, 75131, 577908, 32015],
            [58230, 381139, 78045, 99308, 160454],
            [89135, 80552, 152558, 497981, 603535],
            [78415, 81858, 150656, 193263, 69638],
            [139361, 331509, 343164, 781380, 52269]]

    colLabels = rowLabels = [''] * 5

    rows = len(data)
    ind = np.arange(len(colLabels)) + 0.3  # the x locations for the groups
    cellText = []
    width = 0.4     # the width of the bars
    yoff = np.zeros(len(colLabels))
    # the bottom values for stacked bar chart
    fig, ax = plt.subplots(1, 1)
    for row in range(rows):
        ax.bar(ind, data[row], width, bottom=yoff, color='b')
        yoff = yoff + data[row]
        cellText.append([''])
    plt.xticks([])
    plt.legend([''] * 5, loc=(1.2, 0.2))
    # Add a table at the bottom of the axes
    cellText.reverse()
    the_table = plt.table(cellText=cellText,
                          rowLabels=rowLabels,
                          colLabels=colLabels, loc='bottom')


@image_comparison(baseline_images=['bbox_inches_tight_suptile_legend'],
                  remove_text=False, savefig_kwarg={'bbox_inches': 'tight'})
def test_bbox_inches_tight_suptile_legend():
    plt.plot(np.arange(10), label='a straight line')
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, )
    plt.title('Axis title')
    plt.suptitle('Figure title')

    # put an extra long y tick on to see that the bbox is accounted for
    def y_formatter(y, pos):
        if int(y) == 4:
            return 'The number 4'
        else:
            return str(y)
    plt.gca().yaxis.set_major_formatter(FuncFormatter(y_formatter))

    plt.xlabel('X axis')


@image_comparison(baseline_images=['bbox_inches_tight_clipping'],
                  remove_text=True, savefig_kwarg={'bbox_inches': 'tight'})
def test_bbox_inches_tight_clipping():
    # tests bbox clipping on scatter points, and path clipping on a patch
    # to generate an appropriately tight bbox
    plt.scatter(np.arange(10), np.arange(10))
    ax = plt.gca()
    ax.set_xlim([0, 5])
    ax.set_ylim([0, 5])

    # make a massive rectangle and clip it with a path
    patch = mpatches.Rectangle([-50, -50], 100, 100,
                               transform=ax.transData,
                               facecolor='blue', alpha=0.5)

    path = mpath.Path.unit_regular_star(5).deepcopy()
    path.vertices *= 0.25
    patch.set_clip_path(path, transform=ax.transAxes)
    plt.gcf().artists.append(patch)


@image_comparison(baseline_images=['bbox_inches_tight_raster'],
                  remove_text=True, savefig_kwarg={'bbox_inches': 'tight'})
def test_bbox_inches_tight_raster():
    """Test rasterization with tight_layout"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([1.0, 2.0], rasterized=True)
