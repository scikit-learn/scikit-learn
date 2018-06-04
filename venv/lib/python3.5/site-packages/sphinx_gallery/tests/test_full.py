# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
Test the SG pipeline used with Sphinx
"""
from __future__ import division, absolute_import, print_function

import codecs
import os
import os.path as op
from shutil import copytree, copy2

from sphinx.application import Sphinx
from sphinx_gallery.utils import _TempDir


def test_embed_links():
    """Test that links are embedded properly in doc."""
    tempdir = _TempDir()
    srcdir = op.join(op.dirname(__file__), 'tinybuild')
    for name in os.listdir(srcdir):
        if name in ('_build', 'gen_modules', 'auto_examples'):
            continue
        srcname = op.join(srcdir, name)
        destname = op.join(tempdir, name)
        if op.isdir(srcname):
            copytree(srcname, destname)
        else:
            copy2(srcname, destname)
    out_dir = op.join(tempdir, '_build', 'html')
    toctrees = op.join(tempdir, '_build', 'toctrees')
    # For testing iteration, you can get similar behavior just doing `make`
    # inside the tinybulid directory
    app = Sphinx(tempdir, tempdir, out_dir, toctrees, 'html')
    app.build(False, [])
    out_files = os.listdir(out_dir)
    assert 'index.html' in out_files
    assert 'auto_examples' in out_files
    examples_dir = op.join(out_dir, 'auto_examples')
    assert op.isdir(examples_dir)
    example_files = os.listdir(examples_dir)
    assert 'plot_numpy_scipy.html' in example_files
    example_file = op.join(examples_dir, 'plot_numpy_scipy.html')
    with codecs.open(example_file, 'r', 'utf-8') as fid:
        lines = fid.read()
    # ensure we've linked properly
    assert '#module-scipy.signal' in lines
    assert 'scipy.signal.firwin.html' in lines
    assert '#module-numpy' in lines
    assert 'numpy.arange.html' in lines
    # The matplotlib doc download is a bit unsafe, so skip for now:
    # assert '#module-matplotlib.pyplot' in lines
    # assert 'pyplot_api.html' in lines
