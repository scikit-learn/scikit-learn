from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import functools
import warnings

import matplotlib as mpl
from matplotlib import cbook


def is_called_from_pytest():
    """Returns whether the call was done from pytest"""
    return getattr(mpl, '_called_from_pytest', False)


def _copy_metadata(src_func, tgt_func):
    """Replicates metadata of the function. Returns target function."""
    functools.update_wrapper(tgt_func, src_func)
    tgt_func.__wrapped__ = src_func  # Python2 compatibility.
    return tgt_func


def set_font_settings_for_testing():
    mpl.rcParams['font.family'] = 'DejaVu Sans'
    mpl.rcParams['text.hinting'] = False
    mpl.rcParams['text.hinting_factor'] = 8


def set_reproducibility_for_testing():
    mpl.rcParams['svg.hashsalt'] = 'matplotlib'


def setup():
    # The baseline images are created in this locale, so we should use
    # it during all of the tests.
    import locale
    from matplotlib.backends import backend_agg, backend_pdf, backend_svg

    try:
        locale.setlocale(locale.LC_ALL, str('en_US.UTF-8'))
    except locale.Error:
        try:
            locale.setlocale(locale.LC_ALL, str('English_United States.1252'))
        except locale.Error:
            warnings.warn(
                "Could not set locale to English/United States. "
                "Some date-related tests may fail")

    mpl.use('Agg', warn=False)  # use Agg backend for these tests

    # These settings *must* be hardcoded for running the comparison
    # tests and are not necessarily the default values as specified in
    # rcsetup.py
    mpl.rcdefaults()  # Start with all defaults

    set_font_settings_for_testing()
    set_reproducibility_for_testing()
