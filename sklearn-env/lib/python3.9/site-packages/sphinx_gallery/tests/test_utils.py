# -*- coding: utf-8 -*-
r"""
Test utility functions
==================


"""
# Author: Nicholas Cain
# License: 3-clause BSD

from __future__ import division, absolute_import, print_function
import sphinx_gallery.utils as utils
import pytest


@pytest.mark.parametrize('file_name', ('some/file/name', '/corner.pycase'))
def test_replace_py_ipynb(file_name):
    # Test behavior of function with expected input:
    assert utils.replace_py_ipynb(file_name + '.py') == file_name + '.ipynb'

    # Test behavior of function with unexpected input:
    with pytest.raises(ValueError, match='Unrecognized file extension'):
        utils.replace_py_ipynb(file_name + '.txt')
