"""
Provides utilities to test output reproducibility.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import io
import os
import re
import sys
from subprocess import check_output

import pytest

import matplotlib
from matplotlib import pyplot as plt


def _determinism_save(objects='mhi', format="pdf", usetex=False):
    # save current value of SOURCE_DATE_EPOCH and set it
    # to a constant value, so that time difference is not
    # taken into account
    sde = os.environ.pop('SOURCE_DATE_EPOCH', None)
    os.environ['SOURCE_DATE_EPOCH'] = "946684800"

    matplotlib.rcParams['text.usetex'] = usetex

    fig = plt.figure()

    if 'm' in objects:
        # use different markers...
        ax1 = fig.add_subplot(1, 6, 1)
        x = range(10)
        ax1.plot(x, [1] * 10, marker=u'D')
        ax1.plot(x, [2] * 10, marker=u'x')
        ax1.plot(x, [3] * 10, marker=u'^')
        ax1.plot(x, [4] * 10, marker=u'H')
        ax1.plot(x, [5] * 10, marker=u'v')

    if 'h' in objects:
        # also use different hatch patterns
        ax2 = fig.add_subplot(1, 6, 2)
        bars = (ax2.bar(range(1, 5), range(1, 5)) +
                ax2.bar(range(1, 5), [6] * 4, bottom=range(1, 5)))
        ax2.set_xticks([1.5, 2.5, 3.5, 4.5])

        patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.')
        for bar, pattern in zip(bars, patterns):
            bar.set_hatch(pattern)

    if 'i' in objects:
        # also use different images
        A = [[1, 2, 3], [2, 3, 1], [3, 1, 2]]
        fig.add_subplot(1, 6, 3).imshow(A, interpolation='nearest')
        A = [[1, 3, 2], [1, 2, 3], [3, 1, 2]]
        fig.add_subplot(1, 6, 4).imshow(A, interpolation='bilinear')
        A = [[2, 3, 1], [1, 2, 3], [2, 1, 3]]
        fig.add_subplot(1, 6, 5).imshow(A, interpolation='bicubic')

    x = range(5)
    fig.add_subplot(1, 6, 6).plot(x, x)

    if six.PY2 and format == 'ps':
        stdout = io.StringIO()
    else:
        stdout = getattr(sys.stdout, 'buffer', sys.stdout)
    fig.savefig(stdout, format=format)
    if six.PY2 and format == 'ps':
        sys.stdout.write(stdout.getvalue())

    # Restores SOURCE_DATE_EPOCH
    if sde is None:
        os.environ.pop('SOURCE_DATE_EPOCH', None)
    else:
        os.environ['SOURCE_DATE_EPOCH'] = sde


def _determinism_check(objects='mhi', format="pdf", usetex=False):
    """
    Output three times the same graphs and checks that the outputs are exactly
    the same.

    Parameters
    ----------
    objects : str
        contains characters corresponding to objects to be included in the test
        document: 'm' for markers, 'h' for hatch patterns, 'i' for images. The
        default value is "mhi", so that the test includes all these objects.
    format : str
        format string. The default value is "pdf".
    """
    plots = []
    for i in range(3):
        result = check_output([sys.executable, '-R', '-c',
                               'import matplotlib; '
                               'matplotlib._called_from_pytest = True; '
                               'matplotlib.use(%r); '
                               'from matplotlib.testing.determinism '
                               'import _determinism_save;'
                               '_determinism_save(%r,%r,%r)'
                               % (format, objects, format, usetex)])
        plots.append(result)
    for p in plots[1:]:
        if usetex:
            if p != plots[0]:
                pytest.skip("failed, maybe due to ghostscript timestamps")
        else:
            assert p == plots[0]


def _determinism_source_date_epoch(format, string, keyword=b"CreationDate"):
    """
    Test SOURCE_DATE_EPOCH support. Output a document with the environment
    variable SOURCE_DATE_EPOCH set to 2000-01-01 00:00 UTC and check that the
    document contains the timestamp that corresponds to this date (given as an
    argument).

    Parameters
    ----------
    format : str
        format string, such as "pdf".
    string : str
        timestamp string for 2000-01-01 00:00 UTC.
    keyword : bytes
        a string to look at when searching for the timestamp in the document
        (used in case the test fails).
    """
    buff = check_output([sys.executable, '-R', '-c',
                         'import matplotlib; '
                         'matplotlib._called_from_pytest = True; '
                         'matplotlib.use(%r); '
                         'from matplotlib.testing.determinism '
                         'import _determinism_save;'
                         '_determinism_save(%r,%r)'
                         % (format, "", format)])
    find_keyword = re.compile(b".*" + keyword + b".*")
    key = find_keyword.search(buff)
    if key:
        print(key.group())
    else:
        print("Timestamp keyword (%s) not found!" % keyword)
    assert string in buff
