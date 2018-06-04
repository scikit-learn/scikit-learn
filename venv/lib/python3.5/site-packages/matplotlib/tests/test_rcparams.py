from __future__ import absolute_import, division, print_function

import six

import os
import warnings
from collections import OrderedDict

from cycler import cycler, Cycler
import pytest

try:
    from unittest import mock
except ImportError:
    import mock
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import chain
import numpy as np
from matplotlib.rcsetup import (validate_bool_maybe_none,
                                validate_stringlist,
                                validate_colorlist,
                                validate_color,
                                validate_bool,
                                validate_nseq_int,
                                validate_nseq_float,
                                validate_cycler,
                                validate_hatch,
                                validate_hist_bins,
                                _validate_linestyle)


mpl.rc('text', usetex=False)
mpl.rc('lines', linewidth=22)

fname = os.path.join(os.path.dirname(__file__), 'test_rcparams.rc')


def test_rcparams():
    usetex = mpl.rcParams['text.usetex']
    linewidth = mpl.rcParams['lines.linewidth']

    # test context given dictionary
    with mpl.rc_context(rc={'text.usetex': not usetex}):
        assert mpl.rcParams['text.usetex'] == (not usetex)
    assert mpl.rcParams['text.usetex'] == usetex

    # test context given filename (mpl.rc sets linewdith to 33)
    with mpl.rc_context(fname=fname):
        assert mpl.rcParams['lines.linewidth'] == 33
    assert mpl.rcParams['lines.linewidth'] == linewidth

    # test context given filename and dictionary
    with mpl.rc_context(fname=fname, rc={'lines.linewidth': 44}):
        assert mpl.rcParams['lines.linewidth'] == 44
    assert mpl.rcParams['lines.linewidth'] == linewidth

    # test rc_file
    try:
        mpl.rc_file(fname)
        assert mpl.rcParams['lines.linewidth'] == 33
    finally:
        mpl.rcParams['lines.linewidth'] = linewidth


def test_RcParams_class():
    rc = mpl.RcParams({'font.cursive': ['Apple Chancery',
                                        'Textile',
                                        'Zapf Chancery',
                                        'cursive'],
                       'font.family': 'sans-serif',
                       'font.weight': 'normal',
                       'font.size': 12})

    expected_repr = """
RcParams({'font.cursive': ['Apple Chancery',
                           'Textile',
                           'Zapf Chancery',
                           'cursive'],
          'font.family': ['sans-serif'],
          'font.size': 12.0,
          'font.weight': 'normal'})""".lstrip()

    assert expected_repr == repr(rc)

    expected_str = """
font.cursive: ['Apple Chancery', 'Textile', 'Zapf Chancery', 'cursive']
font.family: ['sans-serif']
font.size: 12.0
font.weight: normal""".lstrip()

    assert expected_str == str(rc)

    # test the find_all functionality
    assert ['font.cursive', 'font.size'] == sorted(rc.find_all('i[vz]'))
    assert ['font.family'] == list(rc.find_all('family'))


def test_rcparams_update():
    rc = mpl.RcParams({'figure.figsize': (3.5, 42)})
    bad_dict = {'figure.figsize': (3.5, 42, 1)}
    # make sure validation happens on input
    with pytest.raises(ValueError):

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',
                                message='.*(validate)',
                                category=UserWarning)
            rc.update(bad_dict)


def test_rcparams_init():
    with pytest.raises(ValueError):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',
                                message='.*(validate)',
                                category=UserWarning)
            mpl.RcParams({'figure.figsize': (3.5, 42, 1)})


def test_Bug_2543():
    # Test that it possible to add all values to itself / deepcopy
    # This was not possible because validate_bool_maybe_none did not
    # accept None as an argument.
    # https://github.com/matplotlib/matplotlib/issues/2543
    # We filter warnings at this stage since a number of them are raised
    # for deprecated rcparams as they should. We don't want these in the
    # printed in the test suite.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                message='.*(deprecated|obsolete)',
                                category=UserWarning)
        with mpl.rc_context():
            _copy = mpl.rcParams.copy()
            for key in _copy:
                mpl.rcParams[key] = _copy[key]
            mpl.rcParams['text.dvipnghack'] = None
        with mpl.rc_context():
            from copy import deepcopy
            _deep_copy = deepcopy(mpl.rcParams)
        # real test is that this does not raise
        assert validate_bool_maybe_none(None) is None
        assert validate_bool_maybe_none("none") is None

    with pytest.raises(ValueError):
        validate_bool_maybe_none("blah")
    with pytest.raises(ValueError):
        validate_bool(None)
    with pytest.raises(ValueError):
        with mpl.rc_context():
            mpl.rcParams['svg.fonttype'] = True


legend_color_tests = [
    ('face', {'color': 'r'}, mcolors.to_rgba('r')),
    ('face', {'color': 'inherit', 'axes.facecolor': 'r'},
     mcolors.to_rgba('r')),
    ('face', {'color': 'g', 'axes.facecolor': 'r'}, mcolors.to_rgba('g')),
    ('edge', {'color': 'r'}, mcolors.to_rgba('r')),
    ('edge', {'color': 'inherit', 'axes.edgecolor': 'r'},
     mcolors.to_rgba('r')),
    ('edge', {'color': 'g', 'axes.facecolor': 'r'}, mcolors.to_rgba('g'))
]
legend_color_test_ids = [
    'same facecolor',
    'inherited facecolor',
    'different facecolor',
    'same edgecolor',
    'inherited edgecolor',
    'different facecolor',
]


@pytest.mark.parametrize('color_type, param_dict, target', legend_color_tests,
                         ids=legend_color_test_ids)
def test_legend_colors(color_type, param_dict, target):
    param_dict['legend.%scolor' % (color_type, )] = param_dict.pop('color')
    get_func = 'get_%scolor' % (color_type, )

    with mpl.rc_context(param_dict):
        _, ax = plt.subplots()
        ax.plot(range(3), label='test')
        leg = ax.legend()
        assert getattr(leg.legendPatch, get_func)() == target


def test_Issue_1713():
    utf32_be = os.path.join(os.path.dirname(__file__),
                           'test_utf32_be_rcparams.rc')
    import locale
    with mock.patch('locale.getpreferredencoding', return_value='UTF-32-BE'):
        rc = mpl.rc_params_from_file(utf32_be, True, False)
    assert rc.get('timezone') == 'UTC'


def generate_validator_testcases(valid):
    validation_tests = (
        {'validator': validate_bool,
         'success': chain(((_, True) for _ in
                           ('t', 'y', 'yes', 'on', 'true', '1', 1, True)),
                           ((_, False) for _ in
                            ('f', 'n', 'no', 'off', 'false', '0', 0, False))),
        'fail': ((_, ValueError)
                 for _ in ('aardvark', 2, -1, [], ))},
        {'validator': validate_stringlist,
         'success': (('', []),
                     ('a,b', ['a', 'b']),
                     ('aardvark', ['aardvark']),
                     ('aardvark, ', ['aardvark']),
                     ('aardvark, ,', ['aardvark']),
                     (['a', 'b'], ['a', 'b']),
                     (('a', 'b'), ['a', 'b']),
                     (iter(['a', 'b']), ['a', 'b']),
                     (np.array(['a', 'b']), ['a', 'b']),
                     ((1, 2), ['1', '2']),
                     (np.array([1, 2]), ['1', '2']),
                    ),
         'fail': ((dict(), ValueError),
                  (1, ValueError),
                 )
        },
        {'validator': validate_nseq_int(2),
         'success': ((_, [1, 2])
                     for _ in ('1, 2', [1.5, 2.5], [1, 2],
                               (1, 2), np.array((1, 2)))),
         'fail': ((_, ValueError)
                  for _ in ('aardvark', ('a', 1),
                            (1, 2, 3)
                            ))
        },
        {'validator': validate_nseq_float(2),
         'success': ((_, [1.5, 2.5])
                     for _ in ('1.5, 2.5', [1.5, 2.5], [1.5, 2.5],
                               (1.5, 2.5), np.array((1.5, 2.5)))),
         'fail': ((_, ValueError)
                  for _ in ('aardvark', ('a', 1),
                            (1, 2, 3)
                            ))
        },
        {'validator': validate_cycler,
         'success': (('cycler("color", "rgb")',
                      cycler("color", 'rgb')),
                     (cycler('linestyle', ['-', '--']),
                      cycler('linestyle', ['-', '--'])),
                     ("""(cycler("color", ["r", "g", "b"]) +
                          cycler("mew", [2, 3, 5]))""",
                      (cycler("color", 'rgb') +
                          cycler("markeredgewidth", [2, 3, 5]))),
                     ("cycler(c='rgb', lw=[1, 2, 3])",
                      cycler('color', 'rgb') + cycler('linewidth', [1, 2, 3])),
                     ("cycler('c', 'rgb') * cycler('linestyle', ['-', '--'])",
                      (cycler('color', 'rgb') *
                          cycler('linestyle', ['-', '--']))),
                     (cycler('ls', ['-', '--']),
                      cycler('linestyle', ['-', '--'])),
                     (cycler(mew=[2, 5]),
                      cycler('markeredgewidth', [2, 5])),
                    ),
         # This is *so* incredibly important: validate_cycler() eval's
         # an arbitrary string! I think I have it locked down enough,
         # and that is what this is testing.
         # TODO: Note that these tests are actually insufficient, as it may
         # be that they raised errors, but still did an action prior to
         # raising the exception. We should devise some additional tests
         # for that...
         'fail': ((4, ValueError),  # Gotta be a string or Cycler object
                  ('cycler("bleh, [])', ValueError),  # syntax error
                  ('Cycler("linewidth", [1, 2, 3])',
                      ValueError),  # only 'cycler()' function is allowed
                  ('1 + 2', ValueError),  # doesn't produce a Cycler object
                  ('os.system("echo Gotcha")', ValueError),  # os not available
                  ('import os', ValueError),  # should not be able to import
                  ('def badjuju(a): return a; badjuju(cycler("color", "rgb"))',
                      ValueError),  # Should not be able to define anything
                                    # even if it does return a cycler
                  ('cycler("waka", [1, 2, 3])', ValueError),  # not a property
                  ('cycler(c=[1, 2, 3])', ValueError),  # invalid values
                  ("cycler(lw=['a', 'b', 'c'])", ValueError),  # invalid values
                  (cycler('waka', [1, 3, 5]), ValueError),  # not a property
                  (cycler('color', ['C1', 'r', 'g']), ValueError)  # no CN
                 )
        },
        {'validator': validate_hatch,
         'success': (('--|', '--|'), ('\\oO', '\\oO'),
                     ('/+*/.x', '/+*/.x'), ('', '')),
         'fail': (('--_', ValueError),
                  (8, ValueError),
                  ('X', ValueError)),
        },
        {'validator': validate_colorlist,
         'success': (('r,g,b', ['r', 'g', 'b']),
                     (['r', 'g', 'b'], ['r', 'g', 'b']),
                     ('r, ,', ['r']),
                     (['', 'g', 'blue'], ['g', 'blue']),
                     ([np.array([1, 0, 0]), np.array([0, 1, 0])],
                         np.array([[1, 0, 0], [0, 1, 0]])),
                     (np.array([[1, 0, 0], [0, 1, 0]]),
                         np.array([[1, 0, 0], [0, 1, 0]])),
                    ),
         'fail': (('fish', ValueError),
                 ),
        },
        {'validator': validate_color,
         'success': (('None', 'none'),
                     ('none', 'none'),
                     ('AABBCC', '#AABBCC'),  # RGB hex code
                     ('AABBCC00', '#AABBCC00'),  # RGBA hex code
                     ('tab:blue', 'tab:blue'),  # named color
                     ('C0', 'C0'),  # color from cycle
                     ('(0, 1, 0)', [0.0, 1.0, 0.0]),  # RGB tuple
                     ((0, 1, 0), (0, 1, 0)),  # non-string version
                     ('(0, 1, 0, 1)', [0.0, 1.0, 0.0, 1.0]),  # RGBA tuple
                     ((0, 1, 0, 1), (0, 1, 0, 1)),  # non-string version
                     ('(0, 1, "0.5")', [0.0, 1.0, 0.5]),  # unusual but valid

                    ),
         'fail': (('tab:veryblue', ValueError),  # invalid name
                  ('C123', ValueError),  # invalid RGB(A) code and cycle index
                  ('(0, 1)', ValueError),  # tuple with length < 3
                  ('(0, 1, 0, 1, 0)', ValueError),  # tuple with length > 4
                  ('(0, 1, none)', ValueError),  # cannot cast none to float
                 ),
        },
        {'validator': validate_hist_bins,
         'success': (('auto', 'auto'),
                     ('10', 10),
                     ('1, 2, 3', [1, 2, 3]),
                     ([1, 2, 3], [1, 2, 3]),
                     (np.arange(15), np.arange(15))
                     ),
         'fail': (('aardvark', ValueError),
                  )
         }
    )

    # The behavior of _validate_linestyle depends on the version of Python.
    # ASCII-compliant bytes arguments should pass on Python 2 because of the
    # automatic conversion between bytes and strings. Python 3 does not
    # perform such a conversion, so the same cases should raise an exception.
    #
    # Common cases:
    ls_test = {'validator': _validate_linestyle,
               'success': (('-', '-'), ('solid', 'solid'),
                           ('--', '--'), ('dashed', 'dashed'),
                           ('-.', '-.'), ('dashdot', 'dashdot'),
                           (':', ':'), ('dotted', 'dotted'),
                           ('', ''), (' ', ' '),
                           ('None', 'none'), ('none', 'none'),
                           ('DoTtEd', 'dotted'),  # case-insensitive
                           (['1.23', '4.56'], (None, [1.23, 4.56])),
                           ([1.23, 456], (None, [1.23, 456.0])),
                           ([1, 2, 3, 4], (None, [1.0, 2.0, 3.0, 4.0])),
                          ),
               'fail': (('aardvark', ValueError),  # not a valid string
                        ('dotted'.encode('utf-16'), ValueError),  # even on PY2
                        ((None, [1, 2]), ValueError),  # (offset, dashes) != OK
                        ((0, [1, 2]), ValueError),  # idem
                        ((-1, [1, 2]), ValueError),  # idem
                        ([1, 2, 3], ValueError),  # sequence with odd length
                        (1.23, ValueError),  # not a sequence
                       )
                }
    # Add some cases of bytes arguments that Python 2 can convert silently:
    ls_bytes_args = (b'dotted', 'dotted'.encode('ascii'))
    if six.PY3:
        ls_test['fail'] += tuple((arg, ValueError) for arg in ls_bytes_args)
    else:
        ls_test['success'] += tuple((arg, 'dotted') for arg in ls_bytes_args)
    # Update the validation test sequence.
    validation_tests += (ls_test,)

    for validator_dict in validation_tests:
        validator = validator_dict['validator']
        if valid:
            for arg, target in validator_dict['success']:
                yield validator, arg, target
        else:
            for arg, error_type in validator_dict['fail']:
                yield validator, arg, error_type


@pytest.mark.parametrize('validator, arg, target',
                         generate_validator_testcases(True))
def test_validator_valid(validator, arg, target):
    res = validator(arg)
    if isinstance(target, np.ndarray):
        assert np.all(res == target)
    elif not isinstance(target, Cycler):
        assert res == target
    else:
        # Cyclers can't simply be asserted equal. They don't implement __eq__
        assert list(res) == list(target)


@pytest.mark.parametrize('validator, arg, exception_type',
                         generate_validator_testcases(False))
def test_validator_invalid(validator, arg, exception_type):
    with pytest.raises(exception_type):
        validator(arg)


def test_keymaps():
    key_list = [k for k in mpl.rcParams if 'keymap' in k]
    for k in key_list:
        assert isinstance(mpl.rcParams[k], list)


def test_rcparams_reset_after_fail():

    # There was previously a bug that meant that if rc_context failed and
    # raised an exception due to issues in the supplied rc parameters, the
    # global rc parameters were left in a modified state.

    with mpl.rc_context(rc={'text.usetex': False}):

        assert mpl.rcParams['text.usetex'] is False

        with pytest.raises(KeyError):
            with mpl.rc_context(rc=OrderedDict([('text.usetex', True),
                                                ('test.blah', True)])):
                pass

        assert mpl.rcParams['text.usetex'] is False


def test_if_rctemplate_is_up_to_date():
    # This tests if the matplotlibrc.template file
    # contains all valid rcParams.
    dep1 = mpl._all_deprecated
    dep2 = mpl._deprecated_set
    deprecated = list(dep1.union(dep2))
    #print(deprecated)
    path_to_rc = mpl.matplotlib_fname()
    with open(path_to_rc, "r") as f:
        rclines = f.readlines()
    missing = {}
    for k,v in mpl.defaultParams.items():
        if k[0] == "_":
            continue
        if k in deprecated:
            continue
        if "verbose" in k:
            continue
        found = False
        for line in rclines:
            if k in line:
                found = True
        if not found:
            missing.update({k:v})
    if missing:
        raise ValueError("The following params are missing " +
                         "in the matplotlibrc.template file: {}"
                         .format(missing.items()))


def test_if_rctemplate_would_be_valid(tmpdir):
    # This tests if the matplotlibrc.template file would result in a valid
    # rc file if all lines are uncommented.
    path_to_rc = mpl.matplotlib_fname()
    with open(path_to_rc, "r") as f:
        rclines = f.readlines()
    newlines = []
    for line in rclines:
        if line[0] == "#":
            newline = line[1:]
        else:
            newline = line
        if "$TEMPLATE_BACKEND" in newline:
            newline = "backend : Agg"
        if "datapath" in newline:
            newline = ""
        newlines.append(newline)
    d = tmpdir.mkdir('test1')
    fname = str(d.join('testrcvalid.temp'))
    with open(fname, "w") as f:
        f.writelines(newlines)
    with pytest.warns(None) as record:
        dic = mpl.rc_params_from_file(fname,
                                      fail_on_error=True,
                                      use_default_template=False)
        assert len(record) == 0
    #d1 = set(dic.keys())
    #d2 = set(matplotlib.defaultParams.keys())
    #print(d2-d1)
