from __future__ import absolute_import, division, print_function

import six

import io
import warnings

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison


needs_usetex = pytest.mark.xfail(
    not matplotlib.checkdep_usetex(True),
    reason="This test needs a TeX installation")


@image_comparison(baseline_images=['font_styles'])
def test_font_styles():
    from matplotlib import _get_data_path
    data_path = _get_data_path()

    def find_matplotlib_font(**kw):
        prop = FontProperties(**kw)
        path = findfont(prop, directory=data_path)
        return FontProperties(fname=path)

    from matplotlib.font_manager import FontProperties, findfont
    warnings.filterwarnings(
        'ignore',
        r"findfont: Font family \[u?'Foo'\] not found. Falling back to .",
        UserWarning,
        module='matplotlib.font_manager')

    plt.figure()
    ax = plt.subplot(1, 1, 1)

    normalFont = find_matplotlib_font(
        family="sans-serif",
        style="normal",
        variant="normal",
        size=14)
    ax.annotate(
        "Normal Font",
        (0.1, 0.1),
        xycoords='axes fraction',
        fontproperties=normalFont)

    boldFont = find_matplotlib_font(
        family="Foo",
        style="normal",
        variant="normal",
        weight="bold",
        stretch=500,
        size=14)
    ax.annotate(
        "Bold Font",
        (0.1, 0.2),
        xycoords='axes fraction',
        fontproperties=boldFont)

    boldItemFont = find_matplotlib_font(
        family="sans serif",
        style="italic",
        variant="normal",
        weight=750,
        stretch=500,
        size=14)
    ax.annotate(
        "Bold Italic Font",
        (0.1, 0.3),
        xycoords='axes fraction',
        fontproperties=boldItemFont)

    lightFont = find_matplotlib_font(
        family="sans-serif",
        style="normal",
        variant="normal",
        weight=200,
        stretch=500,
        size=14)
    ax.annotate(
        "Light Font",
        (0.1, 0.4),
        xycoords='axes fraction',
        fontproperties=lightFont)

    condensedFont = find_matplotlib_font(
        family="sans-serif",
        style="normal",
        variant="normal",
        weight=500,
        stretch=100,
        size=14)
    ax.annotate(
        "Condensed Font",
        (0.1, 0.5),
        xycoords='axes fraction',
        fontproperties=condensedFont)

    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(baseline_images=['multiline'])
def test_multiline():
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.set_title("multiline\ntext alignment")

    plt.text(
        0.2, 0.5, "TpTpTp\n$M$\nTpTpTp", size=20, ha="center", va="top")

    plt.text(
        0.5, 0.5, "TpTpTp\n$M^{M^{M^{M}}}$\nTpTpTp", size=20,
        ha="center", va="top")

    plt.text(
        0.8, 0.5, "TpTpTp\n$M_{q_{q_{q}}}$\nTpTpTp", size=20,
        ha="center", va="top")

    plt.xlim(0, 1)
    plt.ylim(0, 0.8)

    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(baseline_images=['antialiased'], extensions=['png'])
def test_antialiasing():
    matplotlib.rcParams['text.antialiased'] = True

    fig = plt.figure(figsize=(5.25, 0.75))
    fig.text(0.5, 0.75, "antialiased", horizontalalignment='center',
             verticalalignment='center')
    fig.text(0.5, 0.25, r"$\sqrt{x}$", horizontalalignment='center',
             verticalalignment='center')
    # NOTE: We don't need to restore the rcParams here, because the
    # test cleanup will do it for us.  In fact, if we do it here, it
    # will turn antialiasing back off before the images are actually
    # rendered.


def test_afm_kerning():
    from matplotlib.afm import AFM
    from matplotlib.font_manager import findfont

    fn = findfont("Helvetica", fontext="afm")
    with open(fn, 'rb') as fh:
        afm = AFM(fh)
    assert afm.string_width_height('VAVAVAVAVAVA') == (7174.0, 718)


@image_comparison(baseline_images=['text_contains'], extensions=['png'])
def test_contains():
    import matplotlib.backend_bases as mbackend

    fig = plt.figure()
    ax = plt.axes()

    mevent = mbackend.MouseEvent(
        'button_press_event', fig.canvas, 0.5, 0.5, 1, None)

    xs = np.linspace(0.25, 0.75, 30)
    ys = np.linspace(0.25, 0.75, 30)
    xs, ys = np.meshgrid(xs, ys)

    txt = plt.text(
        0.48, 0.52, 'hello world', ha='center', fontsize=30, rotation=30)
    # uncomment to draw the text's bounding box
    # txt.set_bbox(dict(edgecolor='black', facecolor='none'))

    # draw the text. This is important, as the contains method can only work
    # when a renderer exists.
    fig.canvas.draw()

    for x, y in zip(xs.flat, ys.flat):
        mevent.x, mevent.y = plt.gca().transAxes.transform_point([x, y])
        contains, _ = txt.contains(mevent)
        color = 'yellow' if contains else 'red'

        # capture the viewLim, plot a point, and reset the viewLim
        vl = ax.viewLim.frozen()
        ax.plot(x, y, 'o', color=color)
        ax.viewLim.set(vl)


@image_comparison(baseline_images=['titles'])
def test_titles():
    # left and right side titles
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.set_title("left title", loc="left")
    ax.set_title("right title", loc="right")
    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(baseline_images=['text_alignment'])
def test_alignment():
    plt.figure()
    ax = plt.subplot(1, 1, 1)

    x = 0.1
    for rotation in (0, 30):
        for alignment in ('top', 'bottom', 'baseline', 'center'):
            ax.text(
                x, 0.5, alignment + " Tj", va=alignment, rotation=rotation,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            ax.text(
                x, 1.0, r'$\sum_{i=0}^{j}$', va=alignment, rotation=rotation)
            x += 0.1

    ax.plot([0, 1], [0.5, 0.5])
    ax.plot([0, 1], [1.0, 1.0])

    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.5])
    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(baseline_images=['axes_titles'], extensions=['png'])
def test_axes_titles():
    # Related to issue #3327
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.set_title('center', loc='center', fontsize=20, fontweight=700)
    ax.set_title('left', loc='left', fontsize=12, fontweight=400)
    ax.set_title('right', loc='right', fontsize=12, fontweight=400)


def test_set_position():
    fig, ax = plt.subplots()

    # test set_position
    ann = ax.annotate(
        'test', (0, 0), xytext=(0, 0), textcoords='figure pixels')
    fig.canvas.draw()

    init_pos = ann.get_window_extent(fig.canvas.renderer)
    shift_val = 15
    ann.set_position((shift_val, shift_val))
    fig.canvas.draw()
    post_pos = ann.get_window_extent(fig.canvas.renderer)

    for a, b in zip(init_pos.min, post_pos.min):
        assert a + shift_val == b

    # test xyann
    ann = ax.annotate(
        'test', (0, 0), xytext=(0, 0), textcoords='figure pixels')
    fig.canvas.draw()

    init_pos = ann.get_window_extent(fig.canvas.renderer)
    shift_val = 15
    ann.xyann = (shift_val, shift_val)
    fig.canvas.draw()
    post_pos = ann.get_window_extent(fig.canvas.renderer)

    for a, b in zip(init_pos.min, post_pos.min):
        assert a + shift_val == b


def test_get_rotation_string():
    from matplotlib import text
    assert text.get_rotation('horizontal') == 0.
    assert text.get_rotation('vertical') == 90.
    assert text.get_rotation('15.') == 15.


def test_get_rotation_float():
    from matplotlib import text
    for i in [15., 16.70, 77.4]:
        assert text.get_rotation(i) == i


def test_get_rotation_int():
    from matplotlib import text
    for i in [67, 16, 41]:
        assert text.get_rotation(i) == float(i)


def test_get_rotation_raises():
    from matplotlib import text
    with pytest.raises(ValueError):
        text.get_rotation('hozirontal')


def test_get_rotation_none():
    from matplotlib import text
    assert text.get_rotation(None) == 0.0


def test_get_rotation_mod360():
    from matplotlib import text
    for i, j in zip([360., 377., 720+177.2], [0., 17., 177.2]):
        assert_almost_equal(text.get_rotation(i), j)


@image_comparison(baseline_images=['text_bboxclip'])
def test_bbox_clipping():
    plt.text(0.9, 0.2, 'Is bbox clipped?', backgroundcolor='r', clip_on=True)
    t = plt.text(0.9, 0.5, 'Is fancy bbox clipped?', clip_on=True)
    t.set_bbox({"boxstyle": "round, pad=0.1"})


@image_comparison(baseline_images=['annotation_negative_ax_coords'],
                  extensions=['png'])
def test_annotation_negative_ax_coords():
    fig, ax = plt.subplots()

    ax.annotate('+ pts',
                xytext=[30, 20], textcoords='axes points',
                xy=[30, 20], xycoords='axes points', fontsize=32)
    ax.annotate('- pts',
                xytext=[30, -20], textcoords='axes points',
                xy=[30, -20], xycoords='axes points', fontsize=32,
                va='top')
    ax.annotate('+ frac',
                xytext=[0.75, 0.05], textcoords='axes fraction',
                xy=[0.75, 0.05], xycoords='axes fraction', fontsize=32)
    ax.annotate('- frac',
                xytext=[0.75, -0.05], textcoords='axes fraction',
                xy=[0.75, -0.05], xycoords='axes fraction', fontsize=32,
                va='top')

    ax.annotate('+ pixels',
                xytext=[160, 25], textcoords='axes pixels',
                xy=[160, 25], xycoords='axes pixels', fontsize=32)
    ax.annotate('- pixels',
                xytext=[160, -25], textcoords='axes pixels',
                xy=[160, -25], xycoords='axes pixels', fontsize=32,
                va='top')


@image_comparison(baseline_images=['annotation_negative_fig_coords'],
                  extensions=['png'])
def test_annotation_negative_fig_coords():
    fig, ax = plt.subplots()

    ax.annotate('+ pts',
                xytext=[10, 120], textcoords='figure points',
                xy=[10, 120], xycoords='figure points', fontsize=32)
    ax.annotate('- pts',
                xytext=[-10, 180], textcoords='figure points',
                xy=[-10, 180], xycoords='figure points', fontsize=32,
                va='top')
    ax.annotate('+ frac',
                xytext=[0.05, 0.55], textcoords='figure fraction',
                xy=[0.05, 0.55], xycoords='figure fraction', fontsize=32)
    ax.annotate('- frac',
                xytext=[-0.05, 0.5], textcoords='figure fraction',
                xy=[-0.05, 0.5], xycoords='figure fraction', fontsize=32,
                va='top')

    ax.annotate('+ pixels',
                xytext=[50, 50], textcoords='figure pixels',
                xy=[50, 50], xycoords='figure pixels', fontsize=32)
    ax.annotate('- pixels',
                xytext=[-50, 100], textcoords='figure pixels',
                xy=[-50, 100], xycoords='figure pixels', fontsize=32,
                va='top')


def test_text_stale():
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.draw_all()
    assert not ax1.stale
    assert not ax2.stale
    assert not fig.stale

    txt1 = ax1.text(.5, .5, 'aardvark')
    assert ax1.stale
    assert txt1.stale
    assert fig.stale

    ann1 = ax2.annotate('aardvark', xy=[.5, .5])
    assert ax2.stale
    assert ann1.stale
    assert fig.stale

    plt.draw_all()
    assert not ax1.stale
    assert not ax2.stale
    assert not fig.stale


@image_comparison(baseline_images=['agg_text_clip'],
                  extensions=['png'])
def test_agg_text_clip():
    np.random.seed(1)
    fig, (ax1, ax2) = plt.subplots(2)
    for x, y in np.random.rand(10, 2):
        ax1.text(x, y, "foo", clip_on=True)
        ax2.text(x, y, "foo")


def test_text_size_binding():
    from matplotlib.font_manager import FontProperties

    matplotlib.rcParams['font.size'] = 10
    fp = FontProperties(size='large')
    sz1 = fp.get_size_in_points()
    matplotlib.rcParams['font.size'] = 100

    assert sz1 == fp.get_size_in_points()


@image_comparison(baseline_images=['font_scaling'],
                  extensions=['pdf'])
def test_font_scaling():
    matplotlib.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(6.4, 12.4))
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.set_ylim(-10, 600)

    for i, fs in enumerate(range(4, 43, 2)):
        ax.text(0.1, i*30, "{fs} pt font size".format(fs=fs), fontsize=fs)


@pytest.mark.parametrize('spacing1, spacing2', [(0.4, 2), (2, 0.4), (2, 2)])
def test_two_2line_texts(spacing1, spacing2):
    text_string = 'line1\nline2'
    fig = plt.figure()
    renderer = fig.canvas.get_renderer()

    text1 = plt.text(0.25, 0.5, text_string, linespacing=spacing1)
    text2 = plt.text(0.25, 0.5, text_string, linespacing=spacing2)
    fig.canvas.draw()

    box1 = text1.get_window_extent(renderer=renderer)
    box2 = text2.get_window_extent(renderer=renderer)

    # line spacing only affects height
    assert box1.width == box2.width
    if (spacing1 == spacing2):
        assert box1.height == box2.height
    else:
        assert box1.height != box2.height


def test_nonfinite_pos():
    fig, ax = plt.subplots()
    ax.text(0, np.nan, 'nan')
    ax.text(np.inf, 0, 'inf')
    fig.canvas.draw()


def test_hinting_factor_backends():
    plt.rcParams['text.hinting_factor'] = 1
    fig = plt.figure()
    t = fig.text(0.5, 0.5, 'some text')

    fig.savefig(io.BytesIO(), format='svg')
    expected = t.get_window_extent().intervalx

    fig.savefig(io.BytesIO(), format='png')
    # Backends should apply hinting_factor consistently (within 10%).
    np.testing.assert_allclose(t.get_window_extent().intervalx, expected,
                               rtol=0.1)


@needs_usetex
def test_single_artist_usetex():
    # Check that a single artist marked with usetex does not get passed through
    # the mathtext parser at all (for the Agg backend) (the mathtext parser
    # currently fails to parse \frac12, requiring \frac{1}{2} instead).
    fig, ax = plt.subplots()
    ax.text(.5, .5, r"$\frac12$", usetex=True)
    fig.canvas.draw()
