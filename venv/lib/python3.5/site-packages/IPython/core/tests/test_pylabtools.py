"""Tests for pylab tools module.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.


from io import UnsupportedOperation, BytesIO

import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure

from nose import SkipTest
import nose.tools as nt

from matplotlib import pyplot as plt
import numpy as np

from IPython.core.getipython import get_ipython
from IPython.core.interactiveshell import InteractiveShell
from IPython.core.display import _PNG, _JPEG
from .. import pylabtools as pt

from IPython.testing import decorators as dec


def test_figure_to_svg():
    # simple empty-figure test
    fig = plt.figure()
    nt.assert_equal(pt.print_figure(fig, 'svg'), None)

    plt.close('all')

    # simple check for at least svg-looking output
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot([1,2,3])
    plt.draw()
    svg = pt.print_figure(fig, 'svg')[:100].lower()
    nt.assert_in(u'doctype svg', svg)

def _check_pil_jpeg_bytes():
    """Skip if PIL can't write JPEGs to BytesIO objects"""
    # PIL's JPEG plugin can't write to BytesIO objects
    # Pillow fixes this
    from PIL import Image
    buf = BytesIO()
    img = Image.new("RGB", (4,4))
    try:
        img.save(buf, 'jpeg')
    except Exception as e:
        ename = e.__class__.__name__
        raise SkipTest("PIL can't write JPEG to BytesIO: %s: %s" % (ename, e))

@dec.skip_without("PIL.Image")
def test_figure_to_jpeg():
    _check_pil_jpeg_bytes()
    # simple check for at least jpeg-looking output
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot([1,2,3])
    plt.draw()
    jpeg = pt.print_figure(fig, 'jpeg', quality=50)[:100].lower()
    assert jpeg.startswith(_JPEG)

def test_retina_figure():
    # simple empty-figure test
    fig = plt.figure()
    nt.assert_equal(pt.retina_figure(fig), None)
    plt.close('all')

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot([1,2,3])
    plt.draw()
    png, md = pt.retina_figure(fig)
    assert png.startswith(_PNG)
    nt.assert_in('width', md)
    nt.assert_in('height', md)

_fmt_mime_map = {
    'png': 'image/png',
    'jpeg': 'image/jpeg',
    'pdf': 'application/pdf',
    'retina': 'image/png',
    'svg': 'image/svg+xml',
}

def test_select_figure_formats_str():
    ip = get_ipython()
    for fmt, active_mime in _fmt_mime_map.items():
        pt.select_figure_formats(ip, fmt)
        for mime, f in ip.display_formatter.formatters.items():
            if mime == active_mime:
                nt.assert_in(Figure, f)
            else:
                nt.assert_not_in(Figure, f)

def test_select_figure_formats_kwargs():
    ip = get_ipython()
    kwargs = dict(quality=10, bbox_inches='tight')
    pt.select_figure_formats(ip, 'png', **kwargs)
    formatter = ip.display_formatter.formatters['image/png']
    f = formatter.lookup_by_type(Figure)
    cell = f.__closure__[0].cell_contents
    nt.assert_equal(cell, kwargs)
    
    # check that the formatter doesn't raise
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot([1,2,3])
    plt.draw()
    formatter.enabled = True
    png = formatter(fig)
    assert png.startswith(_PNG)

def test_select_figure_formats_set():
    ip = get_ipython()
    for fmts in [
        {'png', 'svg'},
        ['png'],
        ('jpeg', 'pdf', 'retina'),
        {'svg'},
    ]:
        active_mimes = {_fmt_mime_map[fmt] for fmt in fmts}
        pt.select_figure_formats(ip, fmts)
        for mime, f in ip.display_formatter.formatters.items():
            if mime in active_mimes:
                nt.assert_in(Figure, f)
            else:
                nt.assert_not_in(Figure, f)

def test_select_figure_formats_bad():
    ip = get_ipython()
    with nt.assert_raises(ValueError):
        pt.select_figure_formats(ip, 'foo')
    with nt.assert_raises(ValueError):
        pt.select_figure_formats(ip, {'png', 'foo'})
    with nt.assert_raises(ValueError):
        pt.select_figure_formats(ip, ['retina', 'pdf', 'bar', 'bad'])

def test_import_pylab():
    ns = {}
    pt.import_pylab(ns, import_all=False)
    nt.assert_true('plt' in ns)
    nt.assert_equal(ns['np'], np)

class TestPylabSwitch(object):
    class Shell(InteractiveShell):
        def enable_gui(self, gui):
            pass
    
    def setup(self):
        import matplotlib
        def act_mpl(backend):
            matplotlib.rcParams['backend'] = backend

        # Save rcParams since they get modified
        self._saved_rcParams = matplotlib.rcParams
        self._saved_rcParamsOrig = matplotlib.rcParamsOrig
        matplotlib.rcParams = dict(backend='Qt4Agg')
        matplotlib.rcParamsOrig = dict(backend='Qt4Agg')

        # Mock out functions
        self._save_am = pt.activate_matplotlib
        pt.activate_matplotlib = act_mpl
        self._save_ip = pt.import_pylab
        pt.import_pylab = lambda *a,**kw:None
        self._save_cis = pt.configure_inline_support
        pt.configure_inline_support = lambda *a,**kw:None

    def teardown(self):
        pt.activate_matplotlib = self._save_am
        pt.import_pylab = self._save_ip
        pt.configure_inline_support = self._save_cis
        import matplotlib
        matplotlib.rcParams = self._saved_rcParams
        matplotlib.rcParamsOrig = self._saved_rcParamsOrig

    def test_qt(self):
        s = self.Shell()
        gui, backend = s.enable_matplotlib(None)
        nt.assert_equal(gui, 'qt')
        nt.assert_equal(s.pylab_gui_select, 'qt')

        gui, backend = s.enable_matplotlib('inline')
        nt.assert_equal(gui, 'inline')
        nt.assert_equal(s.pylab_gui_select, 'qt')

        gui, backend = s.enable_matplotlib('qt')
        nt.assert_equal(gui, 'qt')
        nt.assert_equal(s.pylab_gui_select, 'qt')

        gui, backend = s.enable_matplotlib('inline')
        nt.assert_equal(gui, 'inline')
        nt.assert_equal(s.pylab_gui_select, 'qt')

        gui, backend = s.enable_matplotlib()
        nt.assert_equal(gui, 'qt')
        nt.assert_equal(s.pylab_gui_select, 'qt')

    def test_inline(self):
        s = self.Shell()
        gui, backend = s.enable_matplotlib('inline')
        nt.assert_equal(gui, 'inline')
        nt.assert_equal(s.pylab_gui_select, None)

        gui, backend = s.enable_matplotlib('inline')
        nt.assert_equal(gui, 'inline')
        nt.assert_equal(s.pylab_gui_select, None)

        gui, backend = s.enable_matplotlib('qt')
        nt.assert_equal(gui, 'qt')
        nt.assert_equal(s.pylab_gui_select, 'qt')

    def test_inline_twice(self):
        "Using '%matplotlib inline' twice should not reset formatters"

        ip = self.Shell()
        gui, backend = ip.enable_matplotlib('inline')
        nt.assert_equal(gui, 'inline')

        fmts =  {'png'}
        active_mimes = {_fmt_mime_map[fmt] for fmt in fmts}
        pt.select_figure_formats(ip, fmts)

        gui, backend = ip.enable_matplotlib('inline')
        nt.assert_equal(gui, 'inline')

        for mime, f in ip.display_formatter.formatters.items():
            if mime in active_mimes:
                nt.assert_in(Figure, f)
            else:
                nt.assert_not_in(Figure, f)

    def test_qt_gtk(self):
        s = self.Shell()
        gui, backend = s.enable_matplotlib('qt')
        nt.assert_equal(gui, 'qt')
        nt.assert_equal(s.pylab_gui_select, 'qt')

        gui, backend = s.enable_matplotlib('gtk')
        nt.assert_equal(gui, 'qt')
        nt.assert_equal(s.pylab_gui_select, 'qt')

