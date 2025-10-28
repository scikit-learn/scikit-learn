"""Pytest fixtures."""

import shutil
from contextlib import contextmanager
from io import StringIO
from pathlib import Path
from unittest.mock import Mock

import pytest
import sphinx
from sphinx.application import Sphinx
from sphinx.errors import ExtensionError
from sphinx.util.docutils import docutils_namespace

from sphinx_gallery import docs_resolv, gen_gallery, gen_rst, py_source_parser
from sphinx_gallery.scrapers import _import_matplotlib
from sphinx_gallery.utils import _get_image

INDEX_RST = """
=============
Own index.rst
=============

Own index.rst file.

.. toctree::

    plot_1
    plot_2
    plot_3
"""

NESTED_PY = """\"\"\"
Header
======

Text.
\"\"\"

a = 1
"""

GALLERY_HEADER = """
Gallery header
==============

Some text.
"""


def pytest_report_header(config, startdir=None):
    """Add information to the pytest run header."""
    return f"Sphinx:  {sphinx.__version__} ({sphinx.__file__})"


@pytest.fixture
def gallery_conf(tmpdir):
    """Set up a test sphinx-gallery configuration."""
    app = Mock(
        spec=Sphinx,
        config=dict(source_suffix={".rst": None}, default_role=None),
        extensions=[],
    )
    gallery_conf = gen_gallery._fill_gallery_conf_defaults({}, app=app)
    gen_gallery._update_gallery_conf_builder_inited(gallery_conf, str(tmpdir))
    gallery_conf.update(examples_dir=str(tmpdir), gallery_dir=str(tmpdir))
    return gallery_conf


@pytest.fixture
def log_collector(monkeypatch):
    app = Mock(spec=Sphinx, name="FakeSphinxApp")()
    monkeypatch.setattr(docs_resolv, "logger", app)
    monkeypatch.setattr(gen_gallery, "logger", app)
    monkeypatch.setattr(py_source_parser, "logger", app)
    monkeypatch.setattr(gen_rst, "logger", app)
    yield app


@pytest.fixture
def unicode_sample(tmpdir):
    """Return temporary python source file with Unicode in various places."""
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

from sphinx_gallery._dummy import DummyClass
DummyClass().prop

from sphinx_gallery._dummy.nested import NestedDummyClass
NestedDummyClass().prop

import matplotlib.pyplot as plt
_ = plt.figure()

"""

    fname = tmpdir.join("unicode_sample.py")
    fname.write(code_str, "wb")
    return fname.strpath


@pytest.fixture
def req_mpl_jpg(tmpdir, req_mpl, scope="session"):
    """Raise SkipTest if JPEG support is not available."""
    # mostly this is needed because of
    # https://github.com/matplotlib/matplotlib/issues/16083
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(range(10))
    try:
        plt.savefig(str(tmpdir.join("testplot.jpg")))
    except Exception as exp:
        pytest.skip(f"Matplotlib jpeg saving failed: {exp}")
    finally:
        plt.close(fig)


@pytest.fixture(scope="session")
def req_mpl():
    try:
        _import_matplotlib()
    except (ImportError, ValueError):
        pytest.skip("Test requires matplotlib")


@pytest.fixture(scope="session")
def req_pil():
    try:
        _get_image()
    except ExtensionError:
        pytest.skip("Test requires pillow")


@pytest.fixture
def conf_file(request):
    try:
        env = request.node.get_closest_marker("add_conf")
    except AttributeError:  # old pytest
        env = request.node.get_marker("add_conf")
    kwargs = env.kwargs if env else {}
    result = {
        "content": "",
        "extensions": ["sphinx_gallery.gen_gallery"],
    }
    result.update(kwargs)

    return result


@pytest.fixture
def rst_file(request):
    try:
        env = request.node.get_closest_marker("add_rst")
    except AttributeError:  # old pytest
        env = request.node.get_marker("add_rst")
    file = env.kwargs["file"] if env else ""
    return file


class SphinxAppWrapper:
    """Wrapper for sphinx.application.Application.

    This allows control over when the sphinx application is initialized, since
    part of the sphinx-gallery build is done in
    sphinx.application.Application.__init__ and the remainder is done in
    sphinx.application.Application.build.
    """

    def __init__(self, srcdir, confdir, outdir, doctreedir, buildername, **kwargs):
        self.srcdir = srcdir
        self.confdir = confdir
        self.outdir = outdir
        self.doctreedir = doctreedir
        self.buildername = buildername
        self.kwargs = kwargs

    def create_sphinx_app(self):
        """Create Sphinx app."""
        # Avoid warnings about re-registration, see:
        # https://github.com/sphinx-doc/sphinx/issues/5038
        with self.create_sphinx_app_context() as app:
            pass
        return app

    @contextmanager
    def create_sphinx_app_context(self):
        """Create Sphinx app inside context."""
        with docutils_namespace():
            app = Sphinx(
                self.srcdir,
                self.confdir,
                self.outdir,
                self.doctreedir,
                self.buildername,
                **self.kwargs,
            )
            yield app

    def build_sphinx_app(self, *args, **kwargs):
        """Build Sphinx app."""
        with self.create_sphinx_app_context() as app:
            # building should be done in the same docutils_namespace context
            app.build(*args, **kwargs)
        return app


@pytest.fixture
def sphinx_app_wrapper(tmpdir, conf_file, rst_file, req_mpl, req_pil):
    _fixturedir = Path(__file__).parent / "testconfs"
    srcdir = Path(tmpdir) / "config_test"
    shutil.copytree(_fixturedir, srcdir)
    # Copy files to 'examples/' as well because default `examples_dirs` is
    # '../examples' - for tests where we don't update config
    shutil.copytree((_fixturedir / "src"), (Path(tmpdir) / "examples"))
    if rst_file == "own index.rst":
        with open((srcdir / "src" / "index.rst"), "w") as file:
            file.write(INDEX_RST)
    elif rst_file:
        with open((srcdir / "minigallery_test.rst"), "w") as file:
            file.write(rst_file)
        # Add nested gallery
        if "sub_folder/sub_sub_folder" in rst_file:
            dir_path = srcdir / "src" / "sub_folder" / "sub_sub_folder"
            dir_path.mkdir(parents=True)
            with open((dir_path / "plot_nested.py"), "w") as file:
                file.write(NESTED_PY)
            with open((dir_path / "GALLERY_HEADER.rst"), "w") as file:
                file.write(GALLERY_HEADER)

    base_config = f"""
import os
import sphinx_gallery
extensions = {conf_file["extensions"]!r}
exclude_patterns = ['_build', 'src']
source_suffix = '.rst'
master_doc = 'index'
# General information about the project.
project = 'Sphinx-Gallery <Tests>'\n\n
"""
    with open((srcdir / "conf.py"), "w") as conffile:
        conffile.write(base_config + conf_file["content"])

    return SphinxAppWrapper(
        srcdir,
        srcdir,
        (srcdir / "_build"),
        (srcdir / "_build" / "toctree"),
        "html",
        warning=StringIO(),
        status=StringIO(),
    )
