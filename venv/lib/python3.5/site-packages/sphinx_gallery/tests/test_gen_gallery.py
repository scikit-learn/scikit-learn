# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
r"""
Test Sphinx-Gallery
"""

from __future__ import (division, absolute_import, print_function,
                        unicode_literals)
import codecs
import os
import re
import shutil
import pytest
from sphinx.application import Sphinx
from sphinx.errors import ExtensionError
from sphinx_gallery.gen_rst import MixedEncodingStringIO
from sphinx_gallery import sphinx_compatibility
from sphinx_gallery.utils import _TempDir


@pytest.fixture
def conf_file(request):
    env = request.node.get_marker('conf_file')
    kwargs = env.kwargs if env else {}
    result = {
        'content': "",
    }
    result.update(kwargs)

    return result


@pytest.fixture
def tempdir():
    """
    temporary directory that wrapped with `path` class.
    this fixture is for compat with old test implementation.
    """
    return _TempDir()


@pytest.fixture
def config_app(tempdir, conf_file):
    _fixturedir = os.path.join(os.path.dirname(__file__), 'testconfs')
    srcdir = os.path.join(tempdir, "config_test")
    shutil.copytree(_fixturedir, srcdir)
    shutil.copytree(os.path.join(_fixturedir, "src"),
                    os.path.join(tempdir, "examples"))

    with open(os.path.join(srcdir, "conf.py"), "w") as conffile:
        conffile.write(conf_file['content'])

    app = Sphinx(srcdir, srcdir, os.path.join(srcdir, "_build"),
                 os.path.join(srcdir, "_build", "toctree"),
                 "html", warning=MixedEncodingStringIO())

    sphinx_compatibility._app = app
    return app


@pytest.mark.conf_file(content="""
import os
import sphinx_gallery
extensions = ['sphinx_gallery.gen_gallery']
# General information about the project.
project = u'Sphinx-Gallery <Tests>'""")
def test_default_config(config_app):
    """Test the default Sphinx-Gallery configuration is loaded

    if only the extension is added to Sphinx"""

    cfg = config_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    # no duplicate values allowed The config is present already
    with pytest.raises(ExtensionError) as excinfo:
        config_app.add_config_value('sphinx_gallery_conf', 'x', True)
    assert 'already present' in str(excinfo.value)


@pytest.mark.conf_file(content="""
import os
import sphinx_gallery
extensions = ['sphinx_gallery.gen_gallery']
# General information about the project.
project = u'Sphinx-Gallery <Tests>'

sphinx_gallery_conf = {
    'mod_example_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_config_old_backreferences_conf(config_app):
    """Testing Deprecation warning message against old backreference config

    In this case the user is required to update the mod_example_dir config
    variable Sphinx-Gallery should notify the user and also silently update
    the old config to the new one. """

    cfg = config_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    assert cfg.sphinx_gallery_conf['backreferences_dir'] == os.path.join(
        'modules', 'gen')
    build_warn = config_app._warning.getvalue()

    assert "WARNING:" in build_warn
    assert "Old configuration" in build_warn
    assert "mod_example_dir" in build_warn
    assert "Gallery now requires" in build_warn
    assert "For a quick fix" in build_warn
    assert "'backreferences_dir': False" not in build_warn


@pytest.mark.conf_file(content="""
import os
import sphinx_gallery
extensions = ['sphinx_gallery.gen_gallery']
# General information about the project.
project = u'Sphinx-Gallery <Tests>'

sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_config_unset_backreferences(config_app):
    """Testing Deprecation warning message against unset backreference config

    In this case the user is notified to update the set the
    backreferences_dir config variable if such feature is to be enabled or
    otherwise to deactivate the feature. Sphinx-Gallery should notify the
    user and also silently setup the old default config value into the new
    config style. """

    cfg = config_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    assert cfg.sphinx_gallery_conf['backreferences_dir'] == os.path.join(
        'modules', 'generated')
    build_warn = config_app._warning.getvalue()

    assert "Gallery now requires" in build_warn
    assert "'backreferences_dir': False" in build_warn
    assert "WARNING:" in build_warn
    assert "mod_example_dir" not in build_warn


@pytest.mark.conf_file(content="""
import os
import sphinx_gallery
extensions = ['sphinx_gallery.gen_gallery']
# General information about the project.
project = u'Sphinx-Gallery <Tests>'

sphinx_gallery_conf = {
    'backreferences_dir': os.path.join('gen_modules', 'backreferences'),
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_config_backreferences(config_app):
    """Test no warning is issued under the new configuration"""

    cfg = config_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    assert cfg.sphinx_gallery_conf['backreferences_dir'] == os.path.join(
        'gen_modules', 'backreferences')
    build_warn = config_app._warning.getvalue()
    assert build_warn == ""


def _check_order(config_app, key):
    index_fname = os.path.join(config_app.outdir, '..', 'ex', 'index.rst')
    order = list()
    regex = '.*:%s=(.):.*' % key
    with codecs.open(index_fname, 'r', 'utf-8') as fid:
        for line in fid:
            if 'sphx-glr-thumbcontainer' in line:
                order.append(int(re.match(regex, line).group(1)))
    assert len(order) == 3
    assert order == [1, 2, 3]


@pytest.mark.conf_file(content="""
import sphinx_gallery
extensions = ['sphinx_gallery.gen_gallery']
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_example_sorting_default(config_app):
    """Test sorting of examples by default key (number of code lines)."""
    _check_order(config_app, 'lines')


@pytest.mark.conf_file(content="""
import sphinx_gallery
from sphinx_gallery.sorting import FileSizeSortKey
extensions = ['sphinx_gallery.gen_gallery']
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'within_subsection_order': FileSizeSortKey,
}""")
def test_example_sorting_filesize(config_app):
    """Test sorting of examples by filesize."""
    _check_order(config_app, 'filesize')


@pytest.mark.conf_file(content="""
import sphinx_gallery
from sphinx_gallery.sorting import FileNameSortKey
extensions = ['sphinx_gallery.gen_gallery']
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'within_subsection_order': FileNameSortKey,
}""")
def test_example_sorting_filename(config_app):
    """Test sorting of examples by filename."""
    _check_order(config_app, 'filename')


@pytest.mark.conf_file(content="""
import sphinx_gallery
from sphinx_gallery.sorting import ExampleTitleSortKey
extensions = ['sphinx_gallery.gen_gallery']
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'within_subsection_order': ExampleTitleSortKey,
}""")
def test_example_sorting_title(config_app):
    """Test sorting of examples by title."""
    _check_order(config_app, 'title')
