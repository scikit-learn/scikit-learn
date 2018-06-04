# -*- coding: utf-8 -*-
r"""
Tests for sorting keys on gallery (sub)sections
===============================================

"""
# Author: Óscar Nájera
# License: 3-clause BSD

from __future__ import division, absolute_import, print_function
import pytest


def test_ExplicitOrder_sorting_key():
    """Test ExplicitOrder"""
    from sphinx_gallery.sorting import ExplicitOrder

    all_folders = ['e', 'f', 'd', 'c', '01b', 'a']
    explicit_folders = ['f', 'd']
    key = ExplicitOrder(explicit_folders)
    sorted_folders = sorted(["d", "f"], key=key)
    assert sorted_folders == explicit_folders

    # Test fails on wrong input
    with pytest.raises(ValueError) as excinfo:
        ExplicitOrder('nope')
    excinfo.match("ExplicitOrder sorting key takes a list")

    # Test missing folder
    with pytest.raises(ValueError) as excinfo:
        sorted_folders = sorted(all_folders, key=key)
    excinfo.match('If you use an explicit folder ordering')
