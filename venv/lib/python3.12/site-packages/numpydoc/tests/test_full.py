import os.path as op
import re
import shutil

import pytest
from docutils import __version__ as docutils_version
from packaging import version
from sphinx.application import Sphinx
from sphinx.util.docutils import docutils_namespace


# Test framework adapted from sphinx-gallery (BSD 3-clause)
@pytest.fixture(scope="module")
def sphinx_app(tmpdir_factory):
    temp_dir = (tmpdir_factory.getbasetemp() / "root").strpath
    src_dir = op.join(op.dirname(__file__), "tinybuild")

    def ignore(src, names):
        return ("_build", "generated")

    shutil.copytree(src_dir, temp_dir, ignore=ignore)
    # For testing iteration, you can get similar behavior just doing `make`
    # inside the tinybuild directory
    src_dir = temp_dir
    conf_dir = temp_dir
    out_dir = op.join(temp_dir, "_build", "html")
    toctrees_dir = op.join(temp_dir, "_build", "toctrees")
    kwargs = {"warningiserror": True, "keep_going": True}
    # Avoid warnings about re-registration, see:
    # https://github.com/sphinx-doc/sphinx/issues/5038
    with docutils_namespace():
        app = Sphinx(
            src_dir, conf_dir, out_dir, toctrees_dir, buildername="html", **kwargs
        )
        # need to build within the context manager
        # for automodule and backrefs to work
        app.build(False, [])
    return app


def test_MyClass(sphinx_app):
    """Test that class documentation is reasonable."""
    src_dir, out_dir = sphinx_app.srcdir, sphinx_app.outdir
    class_rst = op.join(src_dir, "generated", "numpydoc_test_module.MyClass.rst")
    with open(class_rst) as fid:
        rst = fid.read()
    assert r"numpydoc\_test\_module" in rst  # properly escaped
    class_html = op.join(out_dir, "generated", "numpydoc_test_module.MyClass.html")
    with open(class_html) as fid:
        html = fid.read()
    # ensure that no autodoc weirdness ($) occurs
    assert "$self" not in html
    assert "/," not in html
    assert "__init__" in html  # inherited
    # escaped * chars should no longer be preceded by \'s,
    # if we see a \* in the output we know it's incorrect:
    assert r"\*" not in html
    # "self" should not be in the parameter list for the class:
    assert "self," not in html
    # check xref was embedded properly (dict should link using xref):
    assert "stdtypes.html#dict" in html


def test_my_function(sphinx_app):
    """Test that function documentation is reasonable."""
    out_dir = sphinx_app.outdir
    function_html = op.join(
        out_dir, "generated", "numpydoc_test_module.my_function.html"
    )
    with open(function_html) as fid:
        html = fid.read()
    assert r"\*args" not in html
    assert "*args" in html
    # check xref (iterable should link using xref):
    assert "glossary.html#term-iterable" in html


@pytest.mark.parametrize(
    ("html_file", "expected_length"),
    (
        (["index.html"], 1),
        (["generated", "numpydoc_test_module.my_function.html"], 1),
        (["generated", "numpydoc_test_module.MyClass.html"], 1),
    ),
)
def test_reference(sphinx_app, html_file, expected_length):
    """Test for bad references"""
    out_dir = sphinx_app.outdir

    with open(op.join(out_dir, *html_file)) as fid:
        html = fid.read()

    # TODO: This check can be removed when the minimum supported docutils version
    # for numpydoc is docutils>=0.18
    pattern = (
        'role="doc-backlink"'
        if version.parse(docutils_version) >= version.parse("0.18")
        else 'class="fn-backref"'
    )
    reference_list = re.findall(rf'<a {pattern} href="\#id\d+">(.*)<\/a>', html)

    assert len(reference_list) == expected_length
    for ref in reference_list:
        assert "-" not in ref  # Bad reference if it contains "-" e.g. R1896e33633d5-1
