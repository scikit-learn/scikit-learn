r"""Sorters for Sphinx-Gallery (sub)sections.

Sorting key functions for gallery subsection folders and section files.
"""

# Created: Sun May 21 20:38:59 2017
# Author: Óscar Nájera
# License: 3-clause BSD

import os
import types

from sphinx.errors import ConfigError

from .gen_rst import extract_intro_and_title
from .py_source_parser import split_code_and_text_blocks


class ExplicitOrder:
    """Sorting key for all gallery subsections.

    All subsections folders must be listed, otherwise an exception is raised.
    However, you can add '*' as a placeholder to the list. All not-listed
    subsection folders will be inserted at the given position and no
    exception is raised.

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
            raise ConfigError(
                "ExplicitOrder sorting key takes a list, "
                "tuple or Generator, which hold"
                "the paths of each gallery subfolder"
            )

        self.ordered_list = [
            "*" if path == "*" else os.path.normpath(path) for path in ordered_list
        ]
        try:
            i = ordered_list.index("*")
            self.has_wildcard = True
            self.ordered_list_start = self.ordered_list[:i]
            self.ordered_list_end = self.ordered_list[i + 1 :]
        except ValueError:  # from index("*")
            self.has_wildcard = False
            self.ordered_list_start = []
            self.ordered_list_end = self.ordered_list

    def __call__(self, item):
        """
        Return an integer suited for ordering the items.

        If the ordered_list contains a wildcard "*", items before "*" will return
        negative numbers, items after "*" will have positive numbers, and
        not-listed items will return 0.

        If there is no wildcard, all items with return positive numbers, and
        not-listed items will raise a ConfigError.
        """
        if item in self.ordered_list_start:
            return self.ordered_list_start.index(item) - len(self.ordered_list_start)
        elif item in self.ordered_list_end:
            return self.ordered_list_end.index(item) + 1
        else:
            if self.has_wildcard:
                return 0
            else:
                raise ConfigError(
                    "The subsection folder {!r} was not found in the "
                    "'subsection_order' config. If you use an explicit "
                    "'subsection_order', you must specify all subsection folders "
                    "or add '*' as a wildcard to collect all not-listed subsection "
                    "folders.".format(item)
                )

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.ordered_list}>"


class _SortKey:
    """Base class for section order key classes."""

    def __init__(self, src_dir):
        self.src_dir = src_dir

    def __repr__(self):
        return f"<{self.__class__.__name__}>"


class NumberOfCodeLinesSortKey(_SortKey):
    """Sort examples in src_dir by the number of code lines.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        """Return number of code lines in `filename`."""
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        file_conf, script_blocks = split_code_and_text_blocks(src_file)
        amount_of_code = sum(
            len(block.content) for block in script_blocks if block.type == "code"
        )
        return amount_of_code


class FileSizeSortKey(_SortKey):
    """Sort examples in src_dir by file size.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        """Return file size."""
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
        """Return `filename`."""
        return filename


class ExampleTitleSortKey(_SortKey):
    """Sort examples in src_dir by example title.

    Parameters
    ----------
    src_dir : str
        The source directory.
    """

    def __call__(self, filename):
        """Return title of example."""
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        _, script_blocks = split_code_and_text_blocks(src_file)
        _, title = extract_intro_and_title(src_file, script_blocks[0].content)
        return title


class FunctionSortKey:
    """Sort examples using a function passed through to :py:func:`sorted`.

    Parameters
    ----------
    func : :external+python:term:`callable`
           sorting key function,
           can only take one argument, i.e. lambda func = arg: arg[0] * arg[1]
    r : str, None
        printable representation of object
    """

    def __init__(self, func, r=None):
        self.f = func
        self.r = r

    def __repr__(self):
        return self.r if self.r else "FunctionSortKey"

    def __call__(self, arg):
        """Return func(arg)."""
        return self.f(arg)
