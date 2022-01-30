import numpy as np
import matplotlib as mpl
from matplotlib.colors import to_rgb, to_rgba
from numpy.testing import assert_array_equal


LINE_PROPS = [
    "alpha",
    "color",
    "linewidth",
    "linestyle",
    "xydata",
    "zorder",
]

COLLECTION_PROPS = [
    "alpha",
    "edgecolor",
    "facecolor",
    "fill",
    "hatch",
    "linestyle",
    "linewidth",
    "paths",
    "zorder",
]

BAR_PROPS = [
    "alpha",
    "edgecolor",
    "facecolor",
    "fill",
    "hatch",
    "height",
    "linestyle",
    "linewidth",
    "xy",
    "zorder",
]


def assert_colors_equal(a, b, check_alpha=True):

    def handle_array(x):

        if isinstance(x, np.ndarray):
            if x.ndim > 1:
                x = np.unique(x, axis=0).squeeze()
            if x.ndim > 1:
                raise ValueError("Color arrays must be 1 dimensional")
        return x

    a = handle_array(a)
    b = handle_array(b)

    f = to_rgba if check_alpha else to_rgb
    assert f(a) == f(b)


def assert_artists_equal(list1, list2, properties):

    assert len(list1) == len(list2)
    for a1, a2 in zip(list1, list2):
        prop1 = a1.properties()
        prop2 = a2.properties()
        for key in properties:
            v1 = prop1[key]
            v2 = prop2[key]
            if key == "paths":
                for p1, p2 in zip(v1, v2):
                    assert_array_equal(p1.vertices, p2.vertices)
                    assert_array_equal(p1.codes, p2.codes)
            elif isinstance(v1, np.ndarray):
                assert_array_equal(v1, v2)
            elif key == "color":
                v1 = mpl.colors.to_rgba(v1)
                v2 = mpl.colors.to_rgba(v2)
                assert v1 == v2
            else:
                assert v1 == v2


def assert_legends_equal(leg1, leg2):

    assert leg1.get_title().get_text() == leg2.get_title().get_text()
    for t1, t2 in zip(leg1.get_texts(), leg2.get_texts()):
        assert t1.get_text() == t2.get_text()

    assert_artists_equal(
        leg1.get_patches(), leg2.get_patches(), BAR_PROPS,
    )
    assert_artists_equal(
        leg1.get_lines(), leg2.get_lines(), LINE_PROPS,
    )


def assert_plots_equal(ax1, ax2, labels=True):

    assert_artists_equal(ax1.patches, ax2.patches, BAR_PROPS)
    assert_artists_equal(ax1.lines, ax2.lines, LINE_PROPS)

    poly1 = ax1.findobj(mpl.collections.PolyCollection)
    poly2 = ax2.findobj(mpl.collections.PolyCollection)
    assert_artists_equal(poly1, poly2, COLLECTION_PROPS)

    if labels:
        assert ax1.get_xlabel() == ax2.get_xlabel()
        assert ax1.get_ylabel() == ax2.get_ylabel()
