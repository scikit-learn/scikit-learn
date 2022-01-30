# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
r"""
Testing the Jupyter notebook parser
"""

from __future__ import division, absolute_import, print_function
from collections import defaultdict
from itertools import count
import json
import tempfile
import os
import pytest
import re
import base64
import shutil
import textwrap

from sphinx.errors import ExtensionError

import sphinx_gallery.gen_rst as sg
from sphinx_gallery.notebook import (rst2md, jupyter_notebook, save_notebook,
                                     python_to_jupyter_cli)


def test_latex_conversion(gallery_conf):
    """Latex parsing from rst into Jupyter Markdown"""
    double_inline_rst = r":math:`T<0` and :math:`U>0`"
    double_inline_jmd = r"$T<0$ and $U>0$"
    assert double_inline_jmd == rst2md(double_inline_rst, gallery_conf, "", {})

    align_eq = r"""
.. math::
   \mathcal{H} &= 0 \\
   \mathcal{G} &= D"""

    align_eq_jmd = r"""
\begin{align}\mathcal{H} &= 0 \\
   \mathcal{G} &= D\end{align}"""
    assert align_eq_jmd == rst2md(align_eq, gallery_conf, "", {})


def test_convert(gallery_conf):
    """Test ReST conversion"""
    rst = """hello

.. contents::
    :local:

This is :math:`some` math :math:`stuff`.

.. note::
    Interpolation is a linear operation that can be performed also on
    Raw and Epochs objects.

.. warning::
    Go away

For more details on interpolation see the page :ref:`channel_interpolation`.
.. _foo: bar

.. image:: foobar
  :alt: me
  :whatever: you
  :width: 200px
  :class: img_class
"""

    markdown = """hello

This is $some$ math $stuff$.

<div class="alert alert-info"><h4>Note</h4><p>Interpolation is a linear operation that can be performed also on
    Raw and Epochs objects.</p></div>

<div class="alert alert-danger"><h4>Warning</h4><p>Go away</p></div>

For more details on interpolation see the page `channel_interpolation`.

<img src="file://foobar" alt="me" whatever="you" width="200px" class="img_class">
"""  # noqa
    assert rst2md(rst, gallery_conf, "", {}) == markdown


def test_headings():
    rst = textwrap.dedent("""\
    =========
    Heading 1
    =========

    Heading 2
    =========

    =============
     Heading 1-2
    =============

    Heading 3
    ---------

    =============
    Not a Heading
    -------------
    Mismatch top and bottom

    Not another heading
    -=-=-=-=-=-=-=-=-=-
    Multiple characters

    -------
     Bad heading but okay
    -------------
    Over and under mismatch, not rendered and warning raised by Sphinx

    Another bad heading, but passable
    ^^^^^^^^^^^^^^^^^^^^^
    Too short, warning raised but is rendered by Sphinx

    A
    *

    BC
    **

    Some text
    And then a heading
    ------------------
    Not valid with no blank line above

    =======================
           Centered
    =======================

    ------------------------

    ------------------------
    Blank heading above.

    {0}
    ====================
      White space above
    ====================

    """).format('            ')
    # add whitespace afterward to avoid editors from automatically
    # removing the whitespace on save

    heading_level_counter = count(start=1)
    heading_levels = defaultdict(lambda: next(heading_level_counter))
    text = rst2md(rst, {}, "", heading_levels)

    assert text.startswith("# Heading 1\n")
    assert "\n## Heading 2\n" in text
    assert "\n# Heading 1-2\n" in text
    assert "\n### Heading 3\n" in text
    assert "# Not a Heading" not in text
    assert "# Not another Heading" not in text
    assert "\n#### Bad heading but okay\n" in text
    assert "\n##### Another bad heading, but passable\n" in text
    assert "\n###### A\n" in text
    assert "\n###### BC\n" in text
    assert "# And then a heading\n" not in text
    assert "\n# Centered\n" in text
    assert "#\nBlank heading above." not in text
    assert "# White space above\n" in text


@pytest.mark.parametrize(
    'rst_path,md_path,prefix_enabled',
    (('../_static/image.png', 'file://../_static/image.png', False),
     ('/_static/image.png', 'file://_static/image.png', False),
     ('../_static/image.png', 'https://example.com/_static/image.png', True),
     ('/_static/image.png', 'https://example.com/_static/image.png', True),
     ('https://example.com/image.png', 'https://example.com/image.png', False),
     ('https://example.com/image.png', 'https://example.com/image.png', True)),
    ids=('rel_no_prefix', 'abs_no_prefix', 'rel_prefix', 'abs_prefix',
         'url_no_prefix', 'url_prefix'))
def test_notebook_images_prefix(gallery_conf,
                                rst_path, md_path, prefix_enabled):
    if prefix_enabled:
        gallery_conf = gallery_conf.copy()
        gallery_conf['notebook_images'] = "https://example.com/"
    target_dir = os.path.join(
        gallery_conf['src_dir'], gallery_conf['gallery_dirs'])

    rst = textwrap.dedent("""\
    .. image:: {}
       :alt: My Image
       :width: 100px
       :height: 200px
       :class: image
    """).format(rst_path)
    markdown = rst2md(rst, gallery_conf, target_dir, {})

    assert 'src="{}"'.format(md_path) in markdown
    assert 'alt="My Image"' in markdown
    assert 'width="100px"' in markdown
    assert 'height="200px"' in markdown
    assert 'class="image"' in markdown


def test_notebook_images_data_uri(gallery_conf):
    gallery_conf = gallery_conf.copy()
    gallery_conf['notebook_images'] = True
    target_dir = os.path.join(
        gallery_conf['src_dir'], gallery_conf['gallery_dirs'])

    test_image = os.path.join(
        os.path.dirname(__file__), 'tinybuild',
        '_static_nonstandard', 'demo.png')
    # For windows we need to copy this to tmpdir because if tmpdir and this
    # file are on different drives there is no relpath between them
    dest_dir = os.path.join(gallery_conf['src_dir'], '_static_nonstandard')
    os.mkdir(dest_dir)
    dest_image = os.path.join(dest_dir, 'demo.png')
    shutil.copyfile(test_image, dest_image)
    # Make into "absolute" path from source directory
    test_image_rel = os.path.relpath(dest_image, gallery_conf['src_dir'])
    test_image_abs = '/' + test_image_rel.replace(os.sep, '/')
    rst = textwrap.dedent("""\
    .. image:: {}
       :width: 100px
    """).format(test_image_abs)
    markdown = rst2md(rst, gallery_conf, target_dir, {})

    assert 'data' in markdown
    assert 'src="data:image/png;base64,' in markdown
    with open(test_image, 'rb') as test_file:
        data = base64.b64encode(test_file.read())
    assert data.decode('ascii') in markdown

    rst = textwrap.dedent("""\
    .. image:: /this/image/is/missing.png
       :width: 500px
    """)
    with pytest.raises(ExtensionError):
        rst2md(rst, gallery_conf, target_dir, {})


def test_jupyter_notebook(gallery_conf):
    """Test that written ipython notebook file corresponds to python object."""
    file_conf, blocks = sg.split_code_and_text_blocks(
        'tutorials/plot_parse.py')
    target_dir = 'tutorials'
    example_nb = jupyter_notebook(blocks, gallery_conf, target_dir)

    with tempfile.NamedTemporaryFile('w', delete=False) as f:
        save_notebook(example_nb, f.name)
    try:
        with open(f.name, "r") as fname:
            assert json.load(fname) == example_nb
    finally:
        os.remove(f.name)
    assert example_nb.get('cells')[0]['source'][0] == "%matplotlib inline"

    # Test custom first cell text
    test_text = '# testing\n%matplotlib notebook'
    gallery_conf['first_notebook_cell'] = test_text
    example_nb = jupyter_notebook(blocks, gallery_conf, target_dir)
    assert example_nb.get('cells')[0]['source'][0] == test_text

    # Test empty first cell text
    test_text = None
    gallery_conf['first_notebook_cell'] = test_text
    example_nb = jupyter_notebook(blocks, gallery_conf, target_dir)
    cell_src = example_nb.get('cells')[0]['source'][0]
    assert re.match('^[\n]?# Alternating text and code', cell_src)

    # Test custom last cell text
    test_text = '# testing last cell'
    gallery_conf['last_notebook_cell'] = test_text
    example_nb = jupyter_notebook(blocks, gallery_conf, target_dir)
    assert example_nb.get('cells')[-1]['source'][0] == test_text

    # Test empty first cell text
    test_text = None
    gallery_conf['last_notebook_cell'] = test_text
    example_nb = jupyter_notebook(blocks, gallery_conf, target_dir)
    cell_src = example_nb.get('cells')[-1]['source'][0]
    assert re.match("^Last text block.\n\nThat[\\\\]?'s all folks !", cell_src)


###############################################################################
# Notebook shell utility


def test_with_empty_args():
    """ User passes no args, should fail with SystemExit """
    with pytest.raises(SystemExit):
        python_to_jupyter_cli([])


def test_missing_file():
    """ User passes non existing file, should fail with FileNotFoundError """
    with pytest.raises(FileNotFoundError) as excinfo:
        python_to_jupyter_cli(['nofile.py'])
    excinfo.match(r'No such file or directory.+nofile\.py')


def test_file_is_generated():
    """User passes good python file. Check notebook file is created"""

    python_to_jupyter_cli(['examples/plot_0_sin.py'])
    assert os.path.isfile('examples/plot_0_sin.ipynb')
    os.remove('examples/plot_0_sin.ipynb')
