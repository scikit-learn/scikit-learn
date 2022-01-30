# -*- coding: utf-8 -*-
r"""
Sorters for Sphinx-Gallery (sub)sections
========================================

Sorting key functions for gallery subsection folders and section files.
"""
# Created: Sun May 21 20:38:59 2017
# Author: Óscar Nájera
# License: 3-clause BSD

from __future__ import division, absolute_import, print_function
import os
import types

from sphinx.errors import ConfigError

from .gen_rst import extract_intro_and_title
from .py_source_parser import split_code_and_text_blocks


class ExplicitOrder(object):
    """Sorting key for all gallery subsections.

    This requires all folders to be listed otherwise an exception is raised.

    Parameters
    ----------
    ordered_list : list, tuple, or :term:`python:generator`
        Hold the paths of each galleries' subsections.

    Raises
    ------
    ValueError
        Wrong input type or Subgallery path missing.
    """

    def __init__(self, ordered_list):
        if not isinstance(ordered_list, (list, tuple, types.GeneratorType)):
            raise ConfigError("ExplicitOrder sorting key takes a list, "
                              "tuple or Generator, which hold"
                              "the paths of each gallery subfolder")

        self.ordered_list = list(os.path.normpath(path)
                                 for path in ordered_list)

    def __call__(self, item):
        if item in self.ordered_list:
            return self.ordered_list.index(item)
        else:
            raise ConfigError('If you use an explicit folder ordering, you '
                              'must specify all folders. Explicit order not '
                              'found for {}'.format(item))

    def __repr__(self):
        return '<%s : %s>' % (self.__class__.__name__, self.ordered_list)


class _SortKey(object):
    """Base class for section order key classes."""

    def __init__(self, src_dir):
        self.src_dir = src_dir

    def __repr__(self):
        return '<%s>' % (self.__class__.__name__,)


class NumberOfCodeLinesSortKey(_SortKey):
    """Sort examples in src_dir by the number of code lines.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        file_conf, script_blocks = split_code_and_text_blocks(src_file)
        amount_of_code = sum([len(bcontent)
                              for blabel, bcontent, lineno in script_blocks
                              if blabel == 'code'])
        return amount_of_code


class FileSizeSortKey(_SortKey):
    """Sort examples in src_dir by file size.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        return int(os.stat(src_file).st_size)


class FileNameSortKey(_SortKey):
    """Sort examples in src_dir by file name.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        return filename


class ExampleTitleSortKey(_SortKey):
    """Sort examples in src_dir by example title.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        _, script_blocks = split_code_and_text_blocks(src_file)
        _, title = extract_intro_and_title(src_file, script_blocks[0][1])
        return title
