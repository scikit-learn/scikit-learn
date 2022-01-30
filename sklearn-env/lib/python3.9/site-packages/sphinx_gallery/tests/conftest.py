# -*- coding: utf-8 -*-
"""
Pytest fixtures
"""
from __future__ import division, absolute_import, print_function

import collections
from contextlib import contextmanager
from io import StringIO
import os
import shutil

import pytest

import sphinx
from sphinx.application import Sphinx
from sphinx.errors import ExtensionError
from sphinx.util.docutils import docutils_namespace
from sphinx_gallery import (docs_resolv, gen_gallery, gen_rst, utils,
                            sphinx_compatibility, py_source_parser)
from sphinx_gallery.scrapers import _import_matplotlib
from sphinx_gallery.utils import _get_image


def pytest_report_header(config, startdir):
    """Add information to the pytest run header."""
    return 'Sphinx:  %s (%s)' % (sphinx.__version__, sphinx.__file__)


Params = collections.namedtuple('Params', 'args kwargs')


class FakeSphinxApp:
    def __init__(self):
        self.calls = collections.defaultdict(list)

    def status_iterator(self, *args, **kwargs):
        self.calls['status_iterator'].append(Params(args, kwargs))
        for it in args[0]:
            yield it

    def warning(self, *args, **kwargs):
        self.calls['warning'].append(Params(args, kwargs))

    def warn(self, *args, **kwargs):
        self.calls['warn'].append(Params(args, kwargs))

    def info(self, *args, **kwargs):
        self.calls['info'].append(Params(args, kwargs))

    def verbose(self, *args, **kwargs):
        self.calls['verbose'].append(Params(args, kwargs))

    def debug(self, *args, **kwargs):
        self.calls['debug'].append(Params(args, kwargs))


@pytest.fixture
def gallery_conf(tmpdir):
    """Set up a test sphinx-gallery configuration."""
    app = utils.Bunch()
    app.add_css_file = lambda x: None
    app.config = dict(source_suffix={'.rst': None})
    gallery_conf = gen_gallery._complete_gallery_conf(
        {}, str(tmpdir), True, False, app=app)
    gallery_conf.update(examples_dir=str(tmpdir), gallery_dir=str(tmpdir))
    return gallery_conf


@pytest.fixture
def fakesphinxapp():
    orig_app = sphinx_compatibility._app
    sphinx_compatibility._app = app = FakeSphinxApp()
    try:
        yield app
    finally:
        sphinx_compatibility._app = orig_app


@pytest.fixture
def log_collector():
    orig_dr_logger = docs_resolv.logger
    orig_gg_logger = gen_gallery.logger
    orig_gr_logger = gen_rst.logger
    orig_ps_logger = py_source_parser.logger
    app = FakeSphinxApp()
    docs_resolv.logger = app
    gen_gallery.logger = app
    py_source_parser.logger = app
    gen_rst.logger = app
    try:
        yield app
    finally:
        docs_resolv.logger = orig_dr_logger
        gen_gallery.logger = orig_gg_logger
        gen_rst.logger = orig_gr_logger
        py_source_parser.logger = orig_ps_logger


@pytest.fixture
def unicode_sample(tmpdir):
    """Return temporary python source file with Unicode in various places"""
    code_str = b"""# -*- coding: utf-8 -*-
'''
\xc3\x9anicode in header
=================

U\xc3\xb1icode in description
'''

# Code source: \xc3\x93scar N\xc3\xa1jera
# License: BSD 3 clause

import os
path = os.path.join('a','b')

a = 'hei\xc3\x9f'  # Unicode string

import sphinx_gallery.back_references as br
br.identify_names

from sphinx_gallery.back_references import identify_names
identify_names

"""

    fname = tmpdir.join("unicode_sample.py")
    fname.write(code_str, 'wb')
    return fname.strpath


@pytest.fixture
def req_mpl_jpg(tmpdir, req_mpl, scope='session'):
    """Raise SkipTest if JPEG support is not available."""
    # mostly this is needed because of
    # https://github.com/matplotlib/matplotlib/issues/16083
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(range(10))
    try:
        plt.savefig(str(tmpdir.join('testplot.jpg')))
    except Exception as exp:
        pytest.skip('Matplotlib jpeg saving failed: %s' % (exp,))
    finally:
        plt.close(fig)


@pytest.fixture(scope='session')
def req_mpl():
    try:
        _import_matplotlib()
    except (ImportError, ValueError):
        pytest.skip('Test requires matplotlib')


@pytest.fixture(scope='session')
def req_pil():
    try:
        _get_image()
    except ExtensionError:
        pytest.skip('Test requires pillow')


@pytest.fixture
def conf_file(request):
    try:
        env = request.node.get_closest_marker('conf_file')
    except AttributeError:  # old pytest
        env = request.node.get_marker('conf_file')
    kwargs = env.kwargs if env else {}
    result = {
        'content': "",
        'extensions': ['sphinx_gallery.gen_gallery'],
    }
    result.update(kwargs)

    return result


class SphinxAppWrapper(object):
    """Wrapper for sphinx.application.Application.

    This allows control over when the sphinx application is initialized, since
    part of the sphinx-gallery build is done in
    sphinx.application.Application.__init__ and the remainder is done in
    sphinx.application.Application.build.

    """

    def __init__(self, srcdir, confdir, outdir, doctreedir, buildername,
                 **kwargs):
        self.srcdir = srcdir
        self.confdir = confdir
        self.outdir = outdir
        self.doctreedir = doctreedir
        self.buildername = buildername
        self.kwargs = kwargs

    def create_sphinx_app(self):
        # Avoid warnings about re-registration, see:
        # https://github.com/sphinx-doc/sphinx/issues/5038
        with self.create_sphinx_app_context() as app:
            pass
        return app

    @contextmanager
    def create_sphinx_app_context(self):
        with docutils_namespace():
            app = Sphinx(self.srcdir, self.confdir, self.outdir,
                         self.doctreedir, self.buildername, **self.kwargs)
            sphinx_compatibility._app = app
            yield app

    def build_sphinx_app(self, *args, **kwargs):
        with self.create_sphinx_app_context() as app:
            # building should be done in the same docutils_namespace context
            app.build(*args, **kwargs)
        return app


@pytest.fixture
def sphinx_app_wrapper(tmpdir, conf_file, req_mpl, req_pil):
    _fixturedir = os.path.join(os.path.dirname(__file__), 'testconfs')
    srcdir = os.path.join(str(tmpdir), "config_test")
    shutil.copytree(_fixturedir, srcdir)
    shutil.copytree(os.path.join(_fixturedir, "src"),
                    os.path.join(str(tmpdir), "examples"))

    base_config = """
import os
import sphinx_gallery
extensions = %r
exclude_patterns = ['_build']
source_suffix = '.rst'
master_doc = 'index'
# General information about the project.
project = u'Sphinx-Gallery <Tests>'\n\n
""" % (conf_file['extensions'],)
    with open(os.path.join(srcdir, "conf.py"), "w") as conffile:
        conffile.write(base_config + conf_file['content'])

    return SphinxAppWrapper(
        srcdir, srcdir, os.path.join(srcdir, "_build"),
        os.path.join(srcdir, "_build", "toctree"), "html", warning=StringIO(),
        status=StringIO())
