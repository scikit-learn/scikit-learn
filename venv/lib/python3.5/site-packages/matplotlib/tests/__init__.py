from __future__ import absolute_import, division, print_function

import six

import difflib
import os

from matplotlib import cbook
from matplotlib.testing import setup


# Check that the test directories exist
if not os.path.exists(os.path.join(
        os.path.dirname(__file__), 'baseline_images')):
    raise IOError(
        'The baseline image directory does not exist. '
        'This is most likely because the test data is not installed. '
        'You may need to install matplotlib from source to get the '
        'test data.')


@cbook.deprecated("2.1")
def assert_str_equal(reference_str, test_str,
                     format_str=('String {str1} and {str2} do not '
                                 'match:\n{differences}')):
    """
    Assert the two strings are equal. If not, fail and print their
    diffs using difflib.

    """
    if reference_str != test_str:
        diff = difflib.unified_diff(reference_str.splitlines(1),
                                    test_str.splitlines(1),
                                    'Reference', 'Test result',
                                    '', '', 0)
        raise ValueError(format_str.format(str1=reference_str,
                                           str2=test_str,
                                           differences=''.join(diff)))
