"""Tests for IPython.utils.importstring."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import nose.tools as nt

from ..importstring import import_item

def test_import_plain():
    "Test simple imports"
    import os
    os2 = import_item('os')
    nt.assert_true(os is os2)


def test_import_nested():
    "Test nested imports from the stdlib"
    from os import path
    path2 = import_item('os.path')
    nt.assert_true(path is path2)


def test_import_raises():
    "Test that failing imports raise the right exception"
    nt.assert_raises(ImportError, import_item, 'IPython.foobar')
    
