# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import json
import os
import warnings

from unittest import mock

import nose.tools as nt

from IPython.core import display
from IPython.core.getipython import get_ipython
from IPython.utils.io import capture_output
from IPython.utils.tempdir import NamedFileInTemporaryDirectory
from IPython import paths as ipath
from IPython.testing.tools import AssertPrints, AssertNotPrints

import IPython.testing.decorators as dec

def test_image_size():
    """Simple test for display.Image(args, width=x,height=y)"""
    thisurl = 'http://www.google.fr/images/srpr/logo3w.png'
    img = display.Image(url=thisurl, width=200, height=200)
    nt.assert_equal(u'<img src="%s" width="200" height="200"/>' % (thisurl), img._repr_html_())
    img = display.Image(url=thisurl, metadata={'width':200, 'height':200})
    nt.assert_equal(u'<img src="%s" width="200" height="200"/>' % (thisurl), img._repr_html_())
    img = display.Image(url=thisurl, width=200)
    nt.assert_equal(u'<img src="%s" width="200"/>' % (thisurl), img._repr_html_())
    img = display.Image(url=thisurl)
    nt.assert_equal(u'<img src="%s"/>' % (thisurl), img._repr_html_())
    img = display.Image(url=thisurl, unconfined=True)
    nt.assert_equal(u'<img src="%s" class="unconfined"/>' % (thisurl), img._repr_html_())


def test_image_mimes():
    fmt = get_ipython().display_formatter.format
    for format in display.Image._ACCEPTABLE_EMBEDDINGS:
        mime = display.Image._MIMETYPES[format]
        img = display.Image(b'garbage', format=format)
        data, metadata = fmt(img)
        nt.assert_equal(sorted(data), sorted([mime, 'text/plain']))


def test_geojson():

    gj = display.GeoJSON(data={
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [-81.327, 296.038]
            },
            "properties": {
                "name": "Inca City"
            }
        },
        url_template="http://s3-eu-west-1.amazonaws.com/whereonmars.cartodb.net/{basemap_id}/{z}/{x}/{y}.png",
        layer_options={
            "basemap_id": "celestia_mars-shaded-16k_global",
            "attribution": "Celestia/praesepe",
            "minZoom": 0,
            "maxZoom": 18,
        })
    nt.assert_equal(u'<IPython.core.display.GeoJSON object>', str(gj))

def test_retina_png():
    here = os.path.dirname(__file__)
    img = display.Image(os.path.join(here, "2x2.png"), retina=True)
    nt.assert_equal(img.height, 1)
    nt.assert_equal(img.width, 1)
    data, md = img._repr_png_()
    nt.assert_equal(md['width'], 1)
    nt.assert_equal(md['height'], 1)

def test_retina_jpeg():
    here = os.path.dirname(__file__)
    img = display.Image(os.path.join(here, "2x2.jpg"), retina=True)
    nt.assert_equal(img.height, 1)
    nt.assert_equal(img.width, 1)
    data, md = img._repr_jpeg_()
    nt.assert_equal(md['width'], 1)
    nt.assert_equal(md['height'], 1)

def test_base64image():
    display.Image("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABAQMAAAAl21bKAAAAA1BMVEUAAACnej3aAAAAAWJLR0QAiAUdSAAAAAlwSFlzAAALEwAACxMBAJqcGAAAAAd0SU1FB94BCRQnOqNu0b4AAAAKSURBVAjXY2AAAAACAAHiIbwzAAAAAElFTkSuQmCC")

def test_image_filename_defaults():
    '''test format constraint, and validity of jpeg and png'''
    tpath = ipath.get_ipython_package_dir()
    nt.assert_raises(ValueError, display.Image, filename=os.path.join(tpath, 'testing/tests/badformat.zip'),
                     embed=True)
    nt.assert_raises(ValueError, display.Image)
    nt.assert_raises(ValueError, display.Image, data='this is not an image', format='badformat', embed=True)
    # check boths paths to allow packages to test at build and install time
    imgfile = os.path.join(tpath, 'core/tests/2x2.png')
    img = display.Image(filename=imgfile)
    nt.assert_equal('png', img.format)
    nt.assert_is_not_none(img._repr_png_())
    img = display.Image(filename=os.path.join(tpath, 'testing/tests/logo.jpg'), embed=False)
    nt.assert_equal('jpeg', img.format)
    nt.assert_is_none(img._repr_jpeg_())

def _get_inline_config():
    from ipykernel.pylab.config import InlineBackend
    return InlineBackend.instance()

@dec.skip_without('matplotlib')
def test_set_matplotlib_close():
    cfg = _get_inline_config()
    cfg.close_figures = False
    display.set_matplotlib_close()
    assert cfg.close_figures
    display.set_matplotlib_close(False)
    assert not cfg.close_figures

_fmt_mime_map = {
    'png': 'image/png',
    'jpeg': 'image/jpeg',
    'pdf': 'application/pdf',
    'retina': 'image/png',
    'svg': 'image/svg+xml',
}

@dec.skip_without('matplotlib')
def test_set_matplotlib_formats():
    from matplotlib.figure import Figure
    formatters = get_ipython().display_formatter.formatters
    for formats in [
        ('png',),
        ('pdf', 'svg'),
        ('jpeg', 'retina', 'png'),
        (),
    ]:
        active_mimes = {_fmt_mime_map[fmt] for fmt in formats}
        display.set_matplotlib_formats(*formats)
        for mime, f in formatters.items():
            if mime in active_mimes:
                nt.assert_in(Figure, f)
            else:
                nt.assert_not_in(Figure, f)

@dec.skip_without('matplotlib')
def test_set_matplotlib_formats_kwargs():
    from matplotlib.figure import Figure
    ip = get_ipython()
    cfg = _get_inline_config()
    cfg.print_figure_kwargs.update(dict(foo='bar'))
    kwargs = dict(quality=10)
    display.set_matplotlib_formats('png', **kwargs)
    formatter = ip.display_formatter.formatters['image/png']
    f = formatter.lookup_by_type(Figure)
    cell = f.__closure__[0].cell_contents
    expected = kwargs
    expected.update(cfg.print_figure_kwargs)
    nt.assert_equal(cell, expected)

def test_display_available():
    """
    Test that display is available without import

    We don't really care if it's in builtin or anything else, but it should
    always be available.
    """
    ip = get_ipython()
    with AssertNotPrints('NameError'):
        ip.run_cell('display')
    try:
        ip.run_cell('del display')
    except NameError:
        pass # it's ok, it might be in builtins
    # even if deleted it should be back
    with AssertNotPrints('NameError'):
        ip.run_cell('display')

def test_textdisplayobj_pretty_repr():
     p = display.Pretty("This is a simple test")
     nt.assert_equal(repr(p), '<IPython.core.display.Pretty object>')
     nt.assert_equal(p.data, 'This is a simple test')

     p._show_mem_addr = True
     nt.assert_equal(repr(p), object.__repr__(p))

def test_displayobject_repr():
    h = display.HTML('<br />')
    nt.assert_equal(repr(h), '<IPython.core.display.HTML object>')
    h._show_mem_addr = True
    nt.assert_equal(repr(h), object.__repr__(h))
    h._show_mem_addr = False
    nt.assert_equal(repr(h), '<IPython.core.display.HTML object>')

    j = display.Javascript('')
    nt.assert_equal(repr(j), '<IPython.core.display.Javascript object>')
    j._show_mem_addr = True
    nt.assert_equal(repr(j), object.__repr__(j))
    j._show_mem_addr = False
    nt.assert_equal(repr(j), '<IPython.core.display.Javascript object>')

def test_progress():
    p = display.ProgressBar(10)
    nt.assert_in('0/10',repr(p))
    p.html_width = '100%'
    p.progress = 5
    nt.assert_equal(p._repr_html_(), "<progress style='width:100%' max='10' value='5'></progress>")

def test_progress_iter():
    with capture_output(display=False) as captured:
        for i in display.ProgressBar(5):
            out = captured.stdout
            nt.assert_in('{0}/5'.format(i), out)
    out = captured.stdout
    nt.assert_in('5/5', out)

def test_json():
    d = {'a': 5}
    lis = [d]
    md = {'expanded': False}
    md2 = {'expanded': True}
    j = display.JSON(d)
    j2 = display.JSON(d, expanded=True)
    nt.assert_equal(j._repr_json_(), (d, md))
    nt.assert_equal(j2._repr_json_(), (d, md2))

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        j = display.JSON(json.dumps(d))
        nt.assert_equal(len(w), 1)
        nt.assert_equal(j._repr_json_(), (d, md))
        nt.assert_equal(j2._repr_json_(), (d, md2))

    j = display.JSON(lis)
    j2 = display.JSON(lis, expanded=True)
    nt.assert_equal(j._repr_json_(), (lis, md))
    nt.assert_equal(j2._repr_json_(), (lis, md2))

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        j = display.JSON(json.dumps(lis))
        nt.assert_equal(len(w), 1)
        nt.assert_equal(j._repr_json_(), (lis, md))
        nt.assert_equal(j2._repr_json_(), (lis, md2))

def test_video_embedding():
    """use a tempfile, with dummy-data, to ensure that video embedding doesn't crash"""
    v = display.Video("http://ignored")
    assert not v.embed
    html = v._repr_html_()
    nt.assert_not_in('src="data:', html)
    nt.assert_in('src="http://ignored"', html)

    with nt.assert_raises(ValueError):
        v = display.Video(b'abc')

    with NamedFileInTemporaryDirectory('test.mp4') as f:
        f.write(b'abc')
        f.close()

        v = display.Video(f.name)
        assert not v.embed
        html = v._repr_html_()
        nt.assert_not_in('src="data:', html)

        v = display.Video(f.name, embed=True)
        html = v._repr_html_()
        nt.assert_in('src="data:video/mp4;base64,YWJj"',html)

        v = display.Video(f.name, embed=True, mimetype='video/other')
        html = v._repr_html_()
        nt.assert_in('src="data:video/other;base64,YWJj"',html)

        v = display.Video(b'abc', embed=True, mimetype='video/mp4')
        html = v._repr_html_()
        nt.assert_in('src="data:video/mp4;base64,YWJj"',html)

        v = display.Video(u'YWJj', embed=True, mimetype='video/xyz')
        html = v._repr_html_()
        nt.assert_in('src="data:video/xyz;base64,YWJj"',html)

def test_html_metadata():
    s = "<h1>Test</h1>"
    h = display.HTML(s, metadata={"isolated": True})
    nt.assert_equal(h._repr_html_(), (s, {"isolated": True}))

def test_display_id():
    ip = get_ipython()
    with mock.patch.object(ip.display_pub, 'publish') as pub:
        handle = display.display('x')
        nt.assert_is(handle, None)
        handle = display.display('y', display_id='secret')
        nt.assert_is_instance(handle, display.DisplayHandle)
        handle2 = display.display('z', display_id=True)
        nt.assert_is_instance(handle2, display.DisplayHandle)
    nt.assert_not_equal(handle.display_id, handle2.display_id)

    nt.assert_equal(pub.call_count, 3)
    args, kwargs = pub.call_args_list[0]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('x')
        },
        'metadata': {},
    })
    args, kwargs = pub.call_args_list[1]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('y')
        },
        'metadata': {},
        'transient': {
            'display_id': handle.display_id,
        },
    })
    args, kwargs = pub.call_args_list[2]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('z')
        },
        'metadata': {},
        'transient': {
            'display_id': handle2.display_id,
        },
    })


def test_update_display():
    ip = get_ipython()
    with mock.patch.object(ip.display_pub, 'publish') as pub:
        with nt.assert_raises(TypeError):
            display.update_display('x')
        display.update_display('x', display_id='1')
        display.update_display('y', display_id='2')
    args, kwargs = pub.call_args_list[0]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('x')
        },
        'metadata': {},
        'transient': {
            'display_id': '1',
        },
        'update': True,
    })
    args, kwargs = pub.call_args_list[1]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('y')
        },
        'metadata': {},
        'transient': {
            'display_id': '2',
        },
        'update': True,
    })


def test_display_handle():
    ip = get_ipython()
    handle = display.DisplayHandle()
    nt.assert_is_instance(handle.display_id, str)
    handle = display.DisplayHandle('my-id')
    nt.assert_equal(handle.display_id, 'my-id')
    with mock.patch.object(ip.display_pub, 'publish') as pub:
        handle.display('x')
        handle.update('y')

    args, kwargs = pub.call_args_list[0]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('x')
        },
        'metadata': {},
        'transient': {
            'display_id': handle.display_id,
        }
    })
    args, kwargs = pub.call_args_list[1]
    nt.assert_equal(args, ())
    nt.assert_equal(kwargs, {
        'data': {
            'text/plain': repr('y')
        },
        'metadata': {},
        'transient': {
            'display_id': handle.display_id,
        },
        'update': True,
    })

