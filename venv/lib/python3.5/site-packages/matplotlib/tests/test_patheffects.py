from __future__ import absolute_import, division, print_function

import numpy as np
import pytest

from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


@image_comparison(baseline_images=['patheffect1'], remove_text=True)
def test_patheffect1():
    ax1 = plt.subplot(111)
    ax1.imshow([[1, 2], [2, 3]])
    txt = ax1.annotate("test", (1., 1.), (0., 0),
                       arrowprops=dict(arrowstyle="->",
                                       connectionstyle="angle3", lw=2),
                       size=20, ha="center",
                       path_effects=[path_effects.withStroke(linewidth=3,
                                                             foreground="w")])
    txt.arrow_patch.set_path_effects([path_effects.Stroke(linewidth=5,
                                                          foreground="w"),
                                      path_effects.Normal()])

    pe = [path_effects.withStroke(linewidth=3, foreground="w")]
    ax1.grid(True, linestyle="-", path_effects=pe)


@image_comparison(baseline_images=['patheffect2'], remove_text=True,
                  style='mpl20')
def test_patheffect2():

    ax2 = plt.subplot(111)
    arr = np.arange(25).reshape((5, 5))
    ax2.imshow(arr)
    cntr = ax2.contour(arr, colors="k")

    plt.setp(cntr.collections,
             path_effects=[path_effects.withStroke(linewidth=3,
                                                   foreground="w")])

    clbls = ax2.clabel(cntr, fmt="%2.0f", use_clabeltext=True)
    plt.setp(clbls,
             path_effects=[path_effects.withStroke(linewidth=3,
                                                   foreground="w")])


@image_comparison(baseline_images=['patheffect3'])
def test_patheffect3():
    p1, = plt.plot([1, 3, 5, 4, 3], 'o-b', lw=4)
    p1.set_path_effects([path_effects.SimpleLineShadow(),
                         path_effects.Normal()])
    plt.title(r'testing$^{123}$',
        path_effects=[path_effects.withStroke(linewidth=1, foreground="r")])
    leg = plt.legend([p1], [r'Line 1$^2$'], fancybox=True, loc=2)
    leg.legendPatch.set_path_effects([path_effects.withSimplePatchShadow()])

    text = plt.text(2, 3, 'Drop test', color='white',
                    bbox={'boxstyle': 'circle,pad=0.1', 'color': 'red'})
    pe = [path_effects.Stroke(linewidth=3.75, foreground='k'),
          path_effects.withSimplePatchShadow((6, -3), shadow_rgbFace='blue')]
    text.set_path_effects(pe)
    text.get_bbox_patch().set_path_effects(pe)

    pe = [path_effects.PathPatchEffect(offset=(4, -4), hatch='xxxx',
                                       facecolor='gray'),
          path_effects.PathPatchEffect(edgecolor='white', facecolor='black',
                                       lw=1.1)]

    t = plt.gcf().text(0.02, 0.1, 'Hatch shadow', fontsize=75, weight=1000,
                       va='center')
    t.set_path_effects(pe)


@image_comparison(baseline_images=['stroked_text'], extensions=['png'])
def test_patheffects_stroked_text():
    text_chunks = [
        'A B C D E F G H I J K L',
        'M N O P Q R S T U V W',
        'X Y Z a b c d e f g h i j',
        'k l m n o p q r s t u v',
        'w x y z 0123456789',
        r"!@#$%^&*()-=_+[]\;'",
        ',./{}|:"<>?'
    ]
    font_size = 50

    ax = plt.axes([0, 0, 1, 1])
    for i, chunk in enumerate(text_chunks):
        text = ax.text(x=0.01, y=(0.9 - i * 0.13), s=chunk,
                       fontdict={'ha': 'left', 'va': 'center',
                                 'size': font_size, 'color': 'white'})

        text.set_path_effects([path_effects.Stroke(linewidth=font_size / 10,
                                                   foreground='black'),
                               path_effects.Normal()])

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')


@pytest.mark.xfail
def test_PathEffect_points_to_pixels():
    fig = plt.figure(dpi=150)
    p1, = plt.plot(range(10))
    p1.set_path_effects([path_effects.SimpleLineShadow(),
                         path_effects.Normal()])

    renderer = fig.canvas.get_renderer()
    pe_renderer = path_effects.SimpleLineShadow().get_proxy_renderer(renderer)

    assert isinstance(pe_renderer, path_effects.PathEffectRenderer)
    # Confirm that using a path effects renderer maintains point sizes
    # appropriately. Otherwise rendered font would be the wrong size.
    assert renderer.points_to_pixels(15) == pe_renderer.points_to_pixels(15)


def test_SimplePatchShadow_offset():
    pe = path_effects.SimplePatchShadow(offset=(4, 5))
    assert pe._offset == (4, 5)


@image_comparison(baseline_images=['collection'], tol=0.02)
def test_collection():
    x, y = np.meshgrid(np.linspace(0, 10, 150), np.linspace(-5, 5, 100))
    data = np.sin(x) + np.cos(y)
    cs = plt.contour(data)
    pe = [path_effects.PathPatchEffect(edgecolor='black', facecolor='none',
                                       linewidth=12),
          path_effects.Stroke(linewidth=5)]

    for collection in cs.collections:
        collection.set_path_effects(pe)

    for text in plt.clabel(cs, colors='white'):
        text.set_path_effects([path_effects.withStroke(foreground='k',
                                                       linewidth=3)])
        text.set_bbox({'boxstyle': 'sawtooth', 'facecolor': 'none',
                       'edgecolor': 'blue'})
