from __future__ import absolute_import, division, print_function

import six

import numpy as np
from io import BytesIO
import os
import tempfile
import xml.parsers.expat

import pytest

import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison
import matplotlib
from matplotlib import dviread


needs_usetex = pytest.mark.xfail(
    not matplotlib.checkdep_usetex(True),
    reason="This test needs a TeX installation")


def test_visibility():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    x = np.linspace(0, 4 * np.pi, 50)
    y = np.sin(x)
    yerr = np.ones_like(y)

    a, b, c = ax.errorbar(x, y, yerr=yerr, fmt='ko')
    for artist in b:
        artist.set_visible(False)

    fd = BytesIO()
    fig.savefig(fd, format='svg')

    fd.seek(0)
    buf = fd.read()
    fd.close()

    parser = xml.parsers.expat.ParserCreate()
    parser.Parse(buf)  # this will raise ExpatError if the svg is invalid


@image_comparison(baseline_images=['fill_black_with_alpha'], remove_text=True,
                  extensions=['svg'])
def test_fill_black_with_alpha():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(x=[0, 0.1, 1], y=[0, 0, 0], c='k', alpha=0.1, s=10000)


@image_comparison(baseline_images=['noscale'], remove_text=True)
def test_noscale():
    X, Y = np.meshgrid(np.arange(-5, 5, 1), np.arange(-5, 5, 1))
    Z = np.sin(Y ** 2)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(Z, cmap='gray', interpolation='none')


def test_text_urls():
    fig = plt.figure()

    test_url = "http://test_text_urls.matplotlib.org"
    fig.suptitle("test_text_urls", url=test_url)

    fd = BytesIO()
    fig.savefig(fd, format='svg')
    fd.seek(0)
    buf = fd.read().decode()
    fd.close()

    expected = '<a xlink:href="{0}">'.format(test_url)
    assert expected in buf


@image_comparison(baseline_images=['bold_font_output'], extensions=['svg'])
def test_bold_font_output():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(np.arange(10), np.arange(10))
    ax.set_xlabel('nonbold-xlabel')
    ax.set_ylabel('bold-ylabel', fontweight='bold')
    ax.set_title('bold-title', fontweight='bold')


@image_comparison(baseline_images=['bold_font_output_with_none_fonttype'],
                  extensions=['svg'])
def test_bold_font_output_with_none_fonttype():
    plt.rcParams['svg.fonttype'] = 'none'
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(np.arange(10), np.arange(10))
    ax.set_xlabel('nonbold-xlabel')
    ax.set_ylabel('bold-ylabel', fontweight='bold')
    ax.set_title('bold-title', fontweight='bold')


def _test_determinism_save(filename, usetex):
    # This function is mostly copy&paste from "def test_visibility"
    # To require no GUI, we use Figure and FigureCanvasSVG
    # instead of plt.figure and fig.savefig
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_svg import FigureCanvasSVG
    from matplotlib import rc
    rc('svg', hashsalt='asdf')
    rc('text', usetex=usetex)

    fig = Figure()
    ax = fig.add_subplot(111)

    x = np.linspace(0, 4 * np.pi, 50)
    y = np.sin(x)
    yerr = np.ones_like(y)

    a, b, c = ax.errorbar(x, y, yerr=yerr, fmt='ko')
    for artist in b:
        artist.set_visible(False)
    ax.set_title('A string $1+2+\\sigma$')
    ax.set_xlabel('A string $1+2+\\sigma$')
    ax.set_ylabel('A string $1+2+\\sigma$')

    FigureCanvasSVG(fig).print_svg(filename)


@pytest.mark.parametrize(
    "filename, usetex",
    # unique filenames to allow for parallel testing
    [("determinism_notex.svg", False),
     needs_usetex(("determinism_tex.svg", True))])
def test_determinism(filename, usetex):
    import sys
    from subprocess import check_output, STDOUT, CalledProcessError
    plots = []
    for i in range(3):
        # Using check_output and setting stderr to STDOUT will capture the real
        # problem in the output property of the exception
        try:
            check_output(
                [sys.executable, '-R', '-c',
                 'import matplotlib; '
                 'matplotlib._called_from_pytest = True; '
                 'matplotlib.use("svg"); '
                 'from matplotlib.tests.test_backend_svg '
                 'import _test_determinism_save;'
                 '_test_determinism_save(%r, %r)' % (filename, usetex)],
                stderr=STDOUT)
        except CalledProcessError as e:
            # it's easier to use utf8 and ask for forgiveness than try
            # to figure out what the current console has as an
            # encoding :-/
            print(e.output.decode(encoding="utf-8", errors="ignore"))
            raise e
        else:
            with open(filename, 'rb') as fd:
                plots.append(fd.read())
        finally:
            os.unlink(filename)
    for p in plots[1:]:
        assert p == plots[0]


@needs_usetex
def test_missing_psfont(monkeypatch):
    """An error is raised if a TeX font lacks a Type-1 equivalent"""
    from matplotlib import rc

    def psfont(*args, **kwargs):
        return dviread.PsFont(texname='texfont', psname='Some Font',
                              effects=None, encoding=None, filename=None)

    monkeypatch.setattr(dviread.PsfontsMap, '__getitem__', psfont)
    rc('text', usetex=True)
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, 'hello')
    with tempfile.TemporaryFile() as tmpfile, pytest.raises(ValueError):
        fig.savefig(tmpfile, format='svg')
