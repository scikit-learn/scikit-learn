# encoding: utf-8
"""
Utilities for working with strings and text.

Inheritance diagram:

.. inheritance-diagram:: IPython.utils.text
   :parts: 3
"""

import os
import re
import sys
import textwrap
from string import Formatter

# datetime.strftime date format for ipython
if sys.platform == 'win32':
    date_format = "%B %d, %Y"
else:
    date_format = "%B %-d, %Y"


def indent(instr,nspaces=4, ntabs=0, flatten=False):
    """Indent a string a given number of spaces or tabstops.

    indent(str,nspaces=4,ntabs=0) -> indent str by ntabs+nspaces.

    Parameters
    ----------

    instr : basestring
        The string to be indented.
    nspaces : int (default: 4)
        The number of spaces to be indented.
    ntabs : int (default: 0)
        The number of tabs to be indented.
    flatten : bool (default: False)
        Whether to scrub existing indentation.  If True, all lines will be
        aligned to the same indentation.  If False, existing indentation will
        be strictly increased.

    Returns
    -------

    str|unicode : string indented by ntabs and nspaces.

    """
    if instr is None:
        return
    ind = '\t'*ntabs+' '*nspaces
    if flatten:
        pat = re.compile(r'^\s*', re.MULTILINE)
    else:
        pat = re.compile(r'^', re.MULTILINE)
    outstr = re.sub(pat, ind, instr)
    if outstr.endswith(os.linesep+ind):
        return outstr[:-len(ind)]
    else:
        return outstr


def dedent(text):
    """Equivalent of textwrap.dedent that ignores unindented first line.

    This means it will still dedent strings like:
    '''foo
    is a bar
    '''

    For use in wrap_paragraphs.
    """

    if text.startswith('\n'):
        # text starts with blank line, don't ignore the first line
        return textwrap.dedent(text)

    # split first line
    splits = text.split('\n',1)
    if len(splits) == 1:
        # only one line
        return textwrap.dedent(text)

    first, rest = splits
    # dedent everything but the first line
    rest = textwrap.dedent(rest)
    return '\n'.join([first, rest])


def wrap_paragraphs(text, ncols=80):
    """Wrap multiple paragraphs to fit a specified width.

    This is equivalent to textwrap.wrap, but with support for multiple
    paragraphs, as separated by empty lines.

    Returns
    -------

    list of complete paragraphs, wrapped to fill `ncols` columns.
    """
    paragraph_re = re.compile(r'\n(\s*\n)+', re.MULTILINE)
    text = dedent(text).strip()
    paragraphs = paragraph_re.split(text)[::2] # every other entry is space
    out_ps = []
    indent_re = re.compile(r'\n\s+', re.MULTILINE)
    for p in paragraphs:
        # presume indentation that survives dedent is meaningful formatting,
        # so don't fill unless text is flush.
        if indent_re.search(p) is None:
            # wrap paragraph
            p = textwrap.fill(p, ncols)
        out_ps.append(p)
    return out_ps


def strip_ansi(source):
    """
    Remove ansi escape codes from text.
    
    Parameters
    ----------
    source : str
        Source to remove the ansi from
    """
    return re.sub(r'\033\[(\d|;)+?m', '', source)


#-----------------------------------------------------------------------------
# Utils to columnize a list of string
#-----------------------------------------------------------------------------

def _chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


def _find_optimal(rlist , separator_size=2 , displaywidth=80):
    """Calculate optimal info to columnize a list of string"""
    for nrow in range(1, len(rlist)+1) :
        chk = list(map(max,_chunks(rlist, nrow)))
        sumlength = sum(chk)
        ncols = len(chk)
        if sumlength+separator_size*(ncols-1) <= displaywidth :
            break;
    return {'columns_numbers' : ncols,
            'optimal_separator_width':(displaywidth - sumlength)/(ncols-1) if (ncols -1) else 0,
            'rows_numbers' : nrow,
            'columns_width' : chk
           }


def _get_or_default(mylist, i, default=None):
    """return list item number, or default if don't exist"""
    if i >= len(mylist):
        return default
    else :
        return mylist[i]


def compute_item_matrix(items, empty=None, *args, **kwargs) :
    """Returns a nested list, and info to columnize items

    Parameters
    ----------

    items
        list of strings to columize
    empty : (default None)
        default value to fill list if needed
    separator_size : int (default=2)
        How much caracters will be used as a separation between each columns.
    displaywidth : int (default=80)
        The width of the area onto wich the columns should enter

    Returns
    -------

    strings_matrix

        nested list of string, the outer most list contains as many list as
        rows, the innermost lists have each as many element as colums. If the
        total number of elements in `items` does not equal the product of
        rows*columns, the last element of some lists are filled with `None`.

    dict_info
        some info to make columnize easier:

        columns_numbers
          number of columns
        rows_numbers
          number of rows
        columns_width
          list of with of each columns
        optimal_separator_width
          best separator width between columns

    Examples
    --------
    ::

        In [1]: l = ['aaa','b','cc','d','eeeee','f','g','h','i','j','k','l']
           ...: compute_item_matrix(l,displaywidth=12)
        Out[1]:
            ([['aaa', 'f', 'k'],
            ['b', 'g', 'l'],
            ['cc', 'h', None],
            ['d', 'i', None],
            ['eeeee', 'j', None]],
            {'columns_numbers': 3,
            'columns_width': [5, 1, 1],
            'optimal_separator_width': 2,
            'rows_numbers': 5})
    """
    info = _find_optimal(list(map(len, items)), *args, **kwargs)
    nrow, ncol = info['rows_numbers'], info['columns_numbers']
    return ([[ _get_or_default(items, c*nrow+i, default=empty) for c in range(ncol) ] for i in range(nrow) ], info)


def columnize(items, separator='  ', displaywidth=80):
    """ Transform a list of strings into a single string with columns.

    Parameters
    ----------
    items : sequence of strings
        The strings to process.

    separator : str, optional [default is two spaces]
        The string that separates columns.

    displaywidth : int, optional [default is 80]
        Width of the display in number of characters.

    Returns
    -------
    The formatted string.
    """
    if not items :
        return '\n'
    matrix, info = compute_item_matrix(items, separator_size=len(separator), displaywidth=displaywidth)
    fmatrix = [filter(None, x) for x in matrix]
    sjoin = lambda x : separator.join([ y.ljust(w, ' ') for y, w in zip(x, info['columns_width'])])
    return '\n'.join(map(sjoin, fmatrix))+'\n'

