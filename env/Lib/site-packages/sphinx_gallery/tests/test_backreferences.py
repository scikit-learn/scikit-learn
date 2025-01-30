# Author: Óscar Nájera
# License: 3-clause BSD
"""Testing the rst files generator."""

import sys
from unittest.mock import MagicMock, patch

import pytest
from sphinx.errors import ExtensionError

import sphinx_gallery.backreferences as sg
from sphinx_gallery.gen_rst import _sanitize_rst
from sphinx_gallery.py_source_parser import Block, split_code_and_text_blocks

REFERENCE = r"""
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="{0}">

.. only:: html

  .. image:: /fake_dir/images/thumb/sphx_glr_test_file_thumb.png
    :alt:

  :ref:`sphx_glr_fake_dir_test_file.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">test title</div>
    </div>

{1}"""


@pytest.mark.parametrize(
    "content, tooltip, is_backref",
    [
        # HTML sanitizing
        ('<"test">', "&lt;&quot;test&quot;&gt;", False),
        # backref support
        ("test formatting", "test formatting", True),
        # reST sanitizing
        (
            "1 :class:`~a.b`. 2 :class:`a.b` 3 :ref:`whatever <better name>`",
            "1 b. 2 a.b 3 better name",
            False,
        ),
        ("use :meth:`mne.io.Raw.plot_psd` to", "use mne.io.Raw.plot_psd to", False),
        (
            "`this` and ``that``; and `these things` and ``those things``",
            "this and that; and these things and those things",
            False,
        ),
        (
            "See `.MyClass` and `~.MyClass.close`",
            "See MyClass and close",
            False,
        ),
    ],
)
def test_thumbnail_div(content, tooltip, is_backref):
    """Test if the thumbnail div generates the correct string."""
    with pytest.raises(ExtensionError, match="internal Sphinx-Gallery thumb"):
        html_div = sg._thumbnail_div(
            "fake_dir", "", "test_file.py", '<"test">', '<"title">'
        )
    content = _sanitize_rst(content)
    title = "test title"
    html_div = sg._thumbnail_div(
        "fake_dir",
        "",
        "test_file.py",
        content,
        title,
        is_backref=is_backref,
        check=False,
    )
    if is_backref:
        extra = """
.. only:: not html

 * :ref:`sphx_glr_fake_dir_test_file.py`
"""
    else:
        extra = ""
    reference = REFERENCE.format(tooltip, extra)
    assert html_div == reference


def test_identify_names(unicode_sample, gallery_conf):
    """Test name identification."""
    expected = {
        "os.path.join": [
            {
                "name": "join",
                "module": "os.path",
                "module_short": "os.path",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        "br.identify_names": [
            {
                "name": "identify_names",
                "module": "sphinx_gallery.back_references",
                "module_short": "sphinx_gallery.back_references",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        "identify_names": [
            {
                "name": "identify_names",
                "module": "sphinx_gallery.back_references",
                "module_short": "sphinx_gallery.back_references",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        # Check `get_short_module_name` points to correct object.
        # Here, `matplotlib.pyplot.figure` (func) can be shortened to
        # `matplotlib.figure` (module) (but should not be)
        "plt.figure": [
            {
                "name": "figure",
                "module": "matplotlib.pyplot",
                "module_short": "matplotlib.pyplot",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        # Check we get a in a().b in `visit_Attribute`
        "DummyClass": [
            {
                "name": "DummyClass",
                "module": "sphinx_gallery._dummy",
                "module_short": "sphinx_gallery._dummy",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        "NestedDummyClass": [
            {
                "name": "NestedDummyClass",
                "module": "sphinx_gallery._dummy.nested",
                "module_short": "sphinx_gallery._dummy",
                "is_class": False,
                "is_explicit": False,
            }
        ],
    }
    _, script_blocks = split_code_and_text_blocks(unicode_sample)
    ref_regex = sg._make_ref_regex()
    res = sg.identify_names(script_blocks, ref_regex)
    assert expected == res


@pytest.mark.parametrize(("mock", "short_module"), [("same", "A"), ("diff", "A.B")])
def test_get_short_module_name(mock, short_module):
    """Check `_get_short_module_name` correctly finds shortest module."""
    if mock == "same":
        mock_mod_1 = mock_mod_2 = MagicMock()
    else:
        mock_mod_1 = MagicMock()
        mock_mod_2 = MagicMock()

    # Mock will return whatever object/attr we ask of it
    # When mock objects are the same, the `obj_name` from "A" is the same as
    # `obj_name` from "A.B", so "A" is an accepted short module
    with patch.dict(sys.modules, {"A": mock_mod_1, "A.B": mock_mod_2}):
        short_mod = sg._get_short_module_name("A.B", "C")
        short_mod_with_attr = sg._get_short_module_name("A.B", "C.D")
        assert short_mod == short_mod_with_attr == short_module
        # Deleting desired attr will cause failure of `hasattr` check (we
        # check this to ensure it is the right class)
        del mock_mod_2.C.D
        short_mod_no_attr = sg._get_short_module_name("A.B", "C.D")
        assert short_mod_no_attr is None


def test_identify_names_implicit(tmpdir, gallery_conf):
    """Test implicit name identification."""
    code_str = b"""
'''
Title
-----

This is an example.
'''
# -*- coding: utf-8 -*-
# \xc3\x9f
from a.b import c
import d as e
import h.i
print(c)
e.HelloWorld().f.g
h.i.j()
"""
    expected = {
        "c": [
            {
                "name": "c",
                "module": "a.b",
                "module_short": "a.b",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        "e.HelloWorld": [
            {
                "name": "HelloWorld",
                "module": "d",
                "module_short": "d",
                "is_class": False,
                "is_explicit": False,
            }
        ],
        "h.i.j": [
            {
                "name": "j",
                "module": "h.i",
                "module_short": "h.i",
                "is_class": False,
                "is_explicit": False,
            }
        ],
    }

    fname = tmpdir.join("identify_names.py")
    fname.write(code_str, "wb")

    _, script_blocks = split_code_and_text_blocks(fname.strpath)
    ref_regex = sg._make_ref_regex()
    res = sg.identify_names(script_blocks, ref_regex)

    assert expected == res


cobj = dict(module="m", module_short="m", name="n", is_class=False, is_explicit=True)


@pytest.mark.parametrize(
    "text, default_role, ref, cobj",
    [
        (":func:`m.n`", None, "m.n", cobj),
        (":func:`~m.n`", "obj", "m.n", cobj),
        (":func:`Title <m.n>`", None, "m.n", cobj),
        (":func:`!m.n` or `!t <~m.n>`", None, None, None),
        ("`m.n`", "obj", "m.n", cobj),
        ("`m.n`", None, None, None),  # see comment below
        (":ref:`m.n`", None, None, None),
        ("`m.n`", "ref", None, None),
        ("``literal``", "obj", None, None),
    ],
    ids=[
        "regular",
        "show only last component",
        "with title",
        "no link for !",
        "default_role obj",
        "no default_role",  # see comment below
        "non-python role",
        "non-python default_role",
        "literal",
    ],
)
# the sphinx default value for default_role is None = no change, the docutils
# default is title-reference (set by the default-role directive), see
# www.sphinx-doc.org/en/master/usage/configuration.html#confval-default_role
# and docutils.sourceforge.io/docs/ref/rst/roles.html
def test_identify_names_explicit(text, default_role, ref, cobj, gallery_conf):
    """Test explicit name identification."""
    default_role = "" if default_role is None else default_role
    script_blocks = [Block("text", text, 1)]
    expected = {ref: [cobj]} if ref else {}
    ref_regex = sg._make_ref_regex(default_role)
    actual = sg.identify_names(script_blocks, ref_regex)
    assert expected == actual
