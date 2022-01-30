# -*- coding: utf-8 -*-
r"""
Tests for sorting keys on gallery (sub)sections
===============================================

"""
# Author: Óscar Nájera
# License: 3-clause BSD

from __future__ import division, absolute_import, print_function
import os.path as op
import pytest

from sphinx.errors import ConfigError
from sphinx_gallery.sorting import (ExplicitOrder, NumberOfCodeLinesSortKey,
                                    FileNameSortKey, FileSizeSortKey,
                                    ExampleTitleSortKey)


def test_ExplicitOrder_sorting_key():
    """Test ExplicitOrder"""

    all_folders = ['e', 'f', 'd', 'c', '01b', 'a']
    explicit_folders = ['f', 'd']
    key = ExplicitOrder(explicit_folders)
    sorted_folders = sorted(["d", "f"], key=key)
    assert sorted_folders == explicit_folders

    # Test fails on wrong input
    with pytest.raises(ConfigError) as excinfo:
        ExplicitOrder('nope')
    excinfo.match("ExplicitOrder sorting key takes a list")

    # Test missing folder
    with pytest.raises(ConfigError) as excinfo:
        sorted_folders = sorted(all_folders, key=key)
    excinfo.match('If you use an explicit folder ordering')

    # str(obj) stability for sphinx non-rebuilds
    assert str(key).startswith('<ExplicitOrder : ')
    assert str(key) == str(ExplicitOrder(explicit_folders))
    assert str(key) != str(ExplicitOrder(explicit_folders[::-1]))
    src_dir = op.dirname(__file__)
    for klass, type_ in ((NumberOfCodeLinesSortKey, int),
                         (FileNameSortKey, str),
                         (FileSizeSortKey, int),
                         (ExampleTitleSortKey, str)):
        sorter = klass(src_dir)
        assert str(sorter) == '<%s>' % (klass.__name__,)
        out = sorter(op.basename(__file__))
        assert isinstance(out, type_), type(out)
