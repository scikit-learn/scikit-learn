# License: 3-clause BSD
"""Test the SG pipeline using Sphinx and tinybuild."""

import os.path as op
import shutil
from io import StringIO

import pytest
from sphinx.application import Sphinx
from sphinx.util.docutils import docutils_namespace


@pytest.fixture(scope="module")
def sphinx_app(tmpdir_factory, req_mpl, req_pil):
    temp_dir = (tmpdir_factory.getbasetemp() / "root_nonexec").strpath
    src_dir = op.join(op.dirname(__file__), "tinybuild")

    def ignore(src, names):
        return ("_build", "gen_modules", "auto_examples")

    shutil.copytree(src_dir, temp_dir, ignore=ignore)
    # For testing iteration, you can get similar behavior just doing `make`
    # inside the tinybuild/doc directory
    src_dir = temp_dir
    conf_dir = op.join(temp_dir, "doc")
    out_dir = op.join(conf_dir, "_build", "html")
    toctrees_dir = op.join(temp_dir, "doc", "_build", "toctrees")
    # Avoid warnings about re-registration, see:
    # https://github.com/sphinx-doc/sphinx/issues/5038
    confoverrides = {
        "sphinx_gallery_conf.plot_gallery": 0,
    }
    with docutils_namespace():
        app = Sphinx(
            conf_dir,
            conf_dir,
            out_dir,
            toctrees_dir,
            buildername="html",
            confoverrides=confoverrides,
            status=StringIO(),
            warning=StringIO(),
        )
        # need to build within the context manager
        # for automodule and backrefs to work
        app.build(False, [])
    return app


def test_dummy_image(sphinx_app):
    """Test that sphinx_gallery_dummy_images are created."""
    img1 = op.join(
        sphinx_app.srcdir, "auto_examples", "images", "sphx_glr_plot_repr_001.png"
    )
    img2 = op.join(
        sphinx_app.srcdir, "auto_examples", "images", "sphx_glr_plot_repr_002.png"
    )
    assert op.isfile(img1)
    assert op.isfile(img2)
