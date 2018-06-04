# -*- encoding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import os
import shutil

import numpy as np
import pytest

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.compat import subprocess
from matplotlib.testing.compare import compare_images, ImageComparisonFailure
from matplotlib.testing.decorators import image_comparison, _image_directories

baseline_dir, result_dir = _image_directories(lambda: 'dummy func')


def check_for(texsystem):
    header = """
    \\documentclass{minimal}
    \\usepackage{pgf}
    \\begin{document}
    \\typeout{pgfversion=\\pgfversion}
    \\makeatletter
    \\@@end
    """
    try:
        latex = subprocess.Popen([str(texsystem), "-halt-on-error"],
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
        stdout, stderr = latex.communicate(header.encode("utf8"))
    except OSError:
        return False

    return latex.returncode == 0


needs_xelatex = pytest.mark.skipif(not check_for('xelatex'),
                                   reason='xelatex + pgf is required')
needs_pdflatex = pytest.mark.skipif(not check_for('pdflatex'),
                                    reason='pdflatex + pgf is required')


def compare_figure(fname, savefig_kwargs={}, tol=0):
    actual = os.path.join(result_dir, fname)
    plt.savefig(actual, **savefig_kwargs)

    expected = os.path.join(result_dir, "expected_%s" % fname)
    shutil.copyfile(os.path.join(baseline_dir, fname), expected)
    err = compare_images(expected, actual, tol=tol)
    if err:
        raise ImageComparisonFailure(err)


def create_figure():
    plt.figure()
    x = np.linspace(0, 1, 15)

    # line plot
    plt.plot(x, x ** 2, "b-")

    # marker
    plt.plot(x, 1 - x**2, "g>")

    # filled paths and patterns
    plt.fill_between([0., .4], [.4, 0.], hatch='//', facecolor="lightgray",
                     edgecolor="red")
    plt.fill([3, 3, .8, .8, 3], [2, -2, -2, 0, 2], "b")

    # text and typesetting
    plt.plot([0.9], [0.5], "ro", markersize=3)
    plt.text(0.9, 0.5, u'unicode (ü, °, µ) and math ($\\mu_i = x_i^2$)',
             ha='right', fontsize=20)
    plt.ylabel('sans-serif, blue, $\\frac{\\sqrt{x}}{y^2}$..',
               family='sans-serif', color='blue')

    plt.xlim(0, 1)
    plt.ylim(0, 1)


# test compiling a figure to pdf with xelatex
@needs_xelatex
@pytest.mark.backend('pgf')
@image_comparison(baseline_images=['pgf_xelatex'], extensions=['pdf'],
                  style='default')
def test_xelatex():
    rc_xelatex = {'font.family': 'serif',
                  'pgf.rcfonts': False}
    mpl.rcParams.update(rc_xelatex)
    create_figure()


# test compiling a figure to pdf with pdflatex
@needs_pdflatex
@pytest.mark.backend('pgf')
@image_comparison(baseline_images=['pgf_pdflatex'], extensions=['pdf'],
                  style='default')
def test_pdflatex():
    import os
    if os.environ.get('APPVEYOR', False):
        pytest.xfail("pdflatex test does not work on appveyor due to missing "
                     "LaTeX fonts")

    rc_pdflatex = {'font.family': 'serif',
                   'pgf.rcfonts': False,
                   'pgf.texsystem': 'pdflatex',
                   'pgf.preamble': ['\\usepackage[utf8x]{inputenc}',
                                    '\\usepackage[T1]{fontenc}']}
    mpl.rcParams.update(rc_pdflatex)
    create_figure()


# test updating the rc parameters for each figure
@needs_xelatex
@needs_pdflatex
@pytest.mark.style('default')
@pytest.mark.backend('pgf')
def test_rcupdate():
    rc_sets = []
    rc_sets.append({'font.family': 'sans-serif',
                    'font.size': 30,
                    'figure.subplot.left': .2,
                    'lines.markersize': 10,
                    'pgf.rcfonts': False,
                    'pgf.texsystem': 'xelatex'})
    rc_sets.append({'font.family': 'monospace',
                    'font.size': 10,
                    'figure.subplot.left': .1,
                    'lines.markersize': 20,
                    'pgf.rcfonts': False,
                    'pgf.texsystem': 'pdflatex',
                    'pgf.preamble': ['\\usepackage[utf8x]{inputenc}',
                                     '\\usepackage[T1]{fontenc}',
                                     '\\usepackage{sfmath}']})
    tol = (6, 0)
    original_params = mpl.rcParams.copy()
    for i, rc_set in enumerate(rc_sets):
        mpl.rcParams.clear()
        mpl.rcParams.update(original_params)
        mpl.rcParams.update(rc_set)
        create_figure()
        compare_figure('pgf_rcupdate%d.pdf' % (i + 1), tol=tol[i])


# test backend-side clipping, since large numbers are not supported by TeX
@needs_xelatex
@pytest.mark.style('default')
@pytest.mark.backend('pgf')
def test_pathclip():
    rc_xelatex = {'font.family': 'serif',
                  'pgf.rcfonts': False}
    mpl.rcParams.update(rc_xelatex)

    plt.figure()
    plt.plot([0., 1e100], [0., 1e100])
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    # this test passes if compiling/saving to pdf works (no image comparison)
    plt.savefig(os.path.join(result_dir, "pgf_pathclip.pdf"))


# test mixed mode rendering
@needs_xelatex
@pytest.mark.backend('pgf')
@image_comparison(baseline_images=['pgf_mixedmode'], extensions=['pdf'],
                  style='default')
def test_mixedmode():
    rc_xelatex = {'font.family': 'serif',
                  'pgf.rcfonts': False}
    mpl.rcParams.update(rc_xelatex)

    Y, X = np.ogrid[-1:1:40j, -1:1:40j]
    plt.figure()
    plt.pcolor(X**2 + Y**2).set_rasterized(True)


# test bbox_inches clipping
@needs_xelatex
@pytest.mark.style('default')
@pytest.mark.backend('pgf')
def test_bbox_inches():
    rc_xelatex = {'font.family': 'serif',
                  'pgf.rcfonts': False}
    mpl.rcParams.update(rc_xelatex)

    Y, X = np.ogrid[-1:1:40j, -1:1:40j]
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(range(5))
    ax2 = fig.add_subplot(122)
    ax2.plot(range(5))
    plt.tight_layout()

    bbox = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    compare_figure('pgf_bbox_inches.pdf', savefig_kwargs={'bbox_inches': bbox},
                   tol=0)
