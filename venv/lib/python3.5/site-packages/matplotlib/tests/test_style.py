from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile
import warnings
from collections import OrderedDict
from contextlib import contextmanager

import pytest

import matplotlib as mpl
from matplotlib import pyplot as plt, style
from matplotlib.style.core import USER_LIBRARY_PATHS, STYLE_EXTENSION

import six

PARAM = 'image.cmap'
VALUE = 'pink'
DUMMY_SETTINGS = {PARAM: VALUE}


@contextmanager
def temp_style(style_name, settings=None):
    """Context manager to create a style sheet in a temporary directory."""
    if not settings:
        settings = DUMMY_SETTINGS
    temp_file = '%s.%s' % (style_name, STYLE_EXTENSION)

    # Write style settings to file in the temp directory.
    tempdir = tempfile.mkdtemp()
    with open(os.path.join(tempdir, temp_file), 'w') as f:
        for k, v in six.iteritems(settings):
            f.write('%s: %s' % (k, v))

    # Add temp directory to style path and reload so we can access this style.
    USER_LIBRARY_PATHS.append(tempdir)
    style.reload_library()

    try:
        yield
    finally:
        shutil.rmtree(tempdir)
        style.reload_library()


def test_invalid_rc_warning_includes_filename():
    SETTINGS = {'foo': 'bar'}
    basename = 'basename'
    with warnings.catch_warnings(record=True) as warns:
        with temp_style(basename, SETTINGS):
            # style.reload_library() in temp_style() triggers the warning
            pass

    for w in warns:
        assert basename in str(w.message)


def test_available():
    with temp_style('_test_', DUMMY_SETTINGS):
        assert '_test_' in style.available


def test_use():
    mpl.rcParams[PARAM] = 'gray'
    with temp_style('test', DUMMY_SETTINGS):
        with style.context('test'):
            assert mpl.rcParams[PARAM] == VALUE


@pytest.mark.network
def test_use_url():
    with temp_style('test', DUMMY_SETTINGS):
        with style.context('https://gist.github.com/adrn/6590261/raw'):
            assert mpl.rcParams['axes.facecolor'] == "#adeade"


def test_context():
    mpl.rcParams[PARAM] = 'gray'
    with temp_style('test', DUMMY_SETTINGS):
        with style.context('test'):
            assert mpl.rcParams[PARAM] == VALUE
    # Check that this value is reset after the exiting the context.
    assert mpl.rcParams[PARAM] == 'gray'


def test_context_with_dict():
    original_value = 'gray'
    other_value = 'blue'
    mpl.rcParams[PARAM] = original_value
    with style.context({PARAM: other_value}):
        assert mpl.rcParams[PARAM] == other_value
    assert mpl.rcParams[PARAM] == original_value


def test_context_with_dict_after_namedstyle():
    # Test dict after style name where dict modifies the same parameter.
    original_value = 'gray'
    other_value = 'blue'
    mpl.rcParams[PARAM] = original_value
    with temp_style('test', DUMMY_SETTINGS):
        with style.context(['test', {PARAM: other_value}]):
            assert mpl.rcParams[PARAM] == other_value
    assert mpl.rcParams[PARAM] == original_value


def test_context_with_dict_before_namedstyle():
    # Test dict before style name where dict modifies the same parameter.
    original_value = 'gray'
    other_value = 'blue'
    mpl.rcParams[PARAM] = original_value
    with temp_style('test', DUMMY_SETTINGS):
        with style.context([{PARAM: other_value}, 'test']):
            assert mpl.rcParams[PARAM] == VALUE
    assert mpl.rcParams[PARAM] == original_value


def test_context_with_union_of_dict_and_namedstyle():
    # Test dict after style name where dict modifies the a different parameter.
    original_value = 'gray'
    other_param = 'text.usetex'
    other_value = True
    d = {other_param: other_value}
    mpl.rcParams[PARAM] = original_value
    mpl.rcParams[other_param] = (not other_value)
    with temp_style('test', DUMMY_SETTINGS):
        with style.context(['test', d]):
            assert mpl.rcParams[PARAM] == VALUE
            assert mpl.rcParams[other_param] == other_value
    assert mpl.rcParams[PARAM] == original_value
    assert mpl.rcParams[other_param] == (not other_value)


def test_context_with_badparam():
    original_value = 'gray'
    other_value = 'blue'
    d = OrderedDict([(PARAM, original_value), ('badparam', None)])
    with style.context({PARAM: other_value}):
        assert mpl.rcParams[PARAM] == other_value
        x = style.context([d])
        with pytest.raises(KeyError):
            with x:
                pass
        assert mpl.rcParams[PARAM] == other_value


@pytest.mark.parametrize('equiv_styles',
                         [('mpl20', 'default'),
                          ('mpl15', 'classic')],
                         ids=['mpl20', 'mpl15'])
def test_alias(equiv_styles):
    rc_dicts = []
    for sty in equiv_styles:
        with style.context(sty):
            rc_dicts.append(dict(mpl.rcParams))

    rc_base = rc_dicts[0]
    for nm, rc in zip(equiv_styles[1:], rc_dicts[1:]):
        assert rc_base == rc


def test_xkcd_no_cm():
    assert mpl.rcParams["path.sketch"] is None
    plt.xkcd()
    assert mpl.rcParams["path.sketch"] == (1, 100, 2)


def test_xkcd_cm():
    assert mpl.rcParams["path.sketch"] is None
    with plt.xkcd():
        assert mpl.rcParams["path.sketch"] == (1, 100, 2)
    assert mpl.rcParams["path.sketch"] is None
