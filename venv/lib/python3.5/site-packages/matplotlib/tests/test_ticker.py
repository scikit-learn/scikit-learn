from __future__ import absolute_import, division, print_function

from numpy.testing import assert_almost_equal
import numpy as np
import pytest

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import warnings


class TestMaxNLocator(object):
    basic_data = [
        (20, 100, np.array([20., 40., 60., 80., 100.])),
        (0.001, 0.0001, np.array([0., 0.0002, 0.0004, 0.0006, 0.0008, 0.001])),
        (-1e15, 1e15, np.array([-1.0e+15, -5.0e+14, 0e+00, 5e+14, 1.0e+15])),
    ]

    integer_data = [
        (-0.1, 1.1, None, np.array([-1, 0, 1, 2])),
        (-0.1, 0.95, None, np.array([-0.25, 0, 0.25, 0.5, 0.75, 1.0])),
        (1, 55, [1, 1.5, 5, 6, 10], np.array([0, 15, 30, 45, 60])),
    ]

    @pytest.mark.parametrize('vmin, vmax, expected', basic_data)
    def test_basic(self, vmin, vmax, expected):
        loc = mticker.MaxNLocator(nbins=5)
        assert_almost_equal(loc.tick_values(vmin, vmax), expected)

    @pytest.mark.parametrize('vmin, vmax, steps, expected', integer_data)
    def test_integer(self, vmin, vmax, steps, expected):
        loc = mticker.MaxNLocator(nbins=5, integer=True, steps=steps)
        assert_almost_equal(loc.tick_values(vmin, vmax), expected)


class TestLinearLocator(object):
    def test_basic(self):
        loc = mticker.LinearLocator(numticks=3)
        test_value = np.array([-0.8, -0.3, 0.2])
        assert_almost_equal(loc.tick_values(-0.8, 0.2), test_value)

    def test_set_params(self):
        """
        Create linear locator with presets={}, numticks=2 and change it to
        something else. See if change was successful. Should not exception.
        """
        loc = mticker.LinearLocator(numticks=2)
        loc.set_params(numticks=8, presets={(0, 1): []})
        assert loc.numticks == 8
        assert loc.presets == {(0, 1): []}


class TestMultipleLocator(object):
    def test_basic(self):
        loc = mticker.MultipleLocator(base=3.147)
        test_value = np.array([-9.441, -6.294, -3.147, 0., 3.147, 6.294,
                               9.441, 12.588])
        assert_almost_equal(loc.tick_values(-7, 10), test_value)

    def test_set_params(self):
        """
        Create multiple locator with 0.7 base, and change it to something else.
        See if change was successful.
        """
        mult = mticker.MultipleLocator(base=0.7)
        mult.set_params(base=1.7)
        assert mult._base == 1.7


class TestAutoMinorLocator(object):
    def test_basic(self):
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1.39)
        ax.minorticks_on()
        test_value = np.array([0.05, 0.1, 0.15, 0.25, 0.3, 0.35, 0.45,
                               0.5, 0.55, 0.65, 0.7, 0.75, 0.85, 0.9,
                               0.95, 1, 1.05, 1.1, 1.15, 1.25, 1.3, 1.35])
        assert_almost_equal(ax.xaxis.get_ticklocs(minor=True), test_value)

    # NB: the following values are assuming that *xlim* is [0, 5]
    params = [
        (0, 0),  # no major tick => no minor tick either
        (1, 0),  # a single major tick => no minor tick
        (2, 4),  # 1 "nice" major step => 1*5 minor **divisions**
        (3, 6)   # 2 "not nice" major steps => 2*4 minor **divisions**
    ]

    @pytest.mark.parametrize('nb_majorticks, expected_nb_minorticks', params)
    def test_low_number_of_majorticks(
            self, nb_majorticks, expected_nb_minorticks):
        # This test is related to issue #8804
        fig, ax = plt.subplots()
        xlims = (0, 5)  # easier to test the different code paths
        ax.set_xlim(*xlims)
        ax.set_xticks(np.linspace(xlims[0], xlims[1], nb_majorticks))
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
        assert len(ax.xaxis.get_minorticklocs()) == expected_nb_minorticks


class TestLogLocator(object):
    def test_basic(self):
        loc = mticker.LogLocator(numticks=5)
        with pytest.raises(ValueError):
            loc.tick_values(0, 1000)

        test_value = np.array([1.00000000e-05, 1.00000000e-03, 1.00000000e-01,
                               1.00000000e+01, 1.00000000e+03, 1.00000000e+05,
                               1.00000000e+07, 1.000000000e+09])
        assert_almost_equal(loc.tick_values(0.001, 1.1e5), test_value)

        loc = mticker.LogLocator(base=2)
        test_value = np.array([0.5, 1., 2., 4., 8., 16., 32., 64., 128., 256.])
        assert_almost_equal(loc.tick_values(1, 100), test_value)

    def test_set_params(self):
        """
        Create log locator with default value, base=10.0, subs=[1.0],
        numdecs=4, numticks=15 and change it to something else.
        See if change was successful. Should not raise exception.
        """
        loc = mticker.LogLocator()
        loc.set_params(numticks=7, numdecs=8, subs=[2.0], base=4)
        assert loc.numticks == 7
        assert loc.numdecs == 8
        assert loc._base == 4
        assert list(loc._subs) == [2.0]


class TestNullLocator(object):
    def test_set_params(self):
        """
        Create null locator, and attempt to call set_params() on it.
        Should not exception, and should raise a warning.
        """
        loc = mticker.NullLocator()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            loc.set_params()
            assert len(w) == 1


class TestLogitLocator(object):
    def test_set_params(self):
        """
        Create logit locator with default minor=False, and change it to
        something else. See if change was successful. Should not exception.
        """
        loc = mticker.LogitLocator()  # Defaults to false.
        loc.set_params(minor=True)
        assert loc.minor


class TestFixedLocator(object):
    def test_set_params(self):
        """
        Create fixed locator with 5 nbins, and change it to something else.
        See if change was successful.
        Should not exception.
        """
        fixed = mticker.FixedLocator(range(0, 24), nbins=5)
        fixed.set_params(nbins=7)
        assert fixed.nbins == 7


class TestIndexLocator(object):
    def test_set_params(self):
        """
        Create index locator with 3 base, 4 offset. and change it to something
        else. See if change was successful.
        Should not exception.
        """
        index = mticker.IndexLocator(base=3, offset=4)
        index.set_params(base=7, offset=7)
        assert index._base == 7
        assert index.offset == 7


class TestSymmetricalLogLocator(object):
    def test_set_params(self):
        """
        Create symmetrical log locator with default subs =[1.0] numticks = 15,
        and change it to something else.
        See if change was successful.
        Should not exception.
        """
        sym = mticker.SymmetricalLogLocator(base=10, linthresh=1)
        sym.set_params(subs=[2.0], numticks=8)
        assert sym._subs == [2.0]
        assert sym.numticks == 8


class TestScalarFormatter(object):
    offset_data = [
        (123, 189, 0),
        (-189, -123, 0),
        (12341, 12349, 12340),
        (-12349, -12341, -12340),
        (99999.5, 100010.5, 100000),
        (-100010.5, -99999.5, -100000),
        (99990.5, 100000.5, 100000),
        (-100000.5, -99990.5, -100000),
        (1233999, 1234001, 1234000),
        (-1234001, -1233999, -1234000),
        (1, 1, 1),
        (123, 123, 120),
        # Test cases courtesy of @WeatherGod
        (.4538, .4578, .45),
        (3789.12, 3783.1, 3780),
        (45124.3, 45831.75, 45000),
        (0.000721, 0.0007243, 0.00072),
        (12592.82, 12591.43, 12590),
        (9., 12., 0),
        (900., 1200., 0),
        (1900., 1200., 0),
        (0.99, 1.01, 1),
        (9.99, 10.01, 10),
        (99.99, 100.01, 100),
        (5.99, 6.01, 6),
        (15.99, 16.01, 16),
        (-0.452, 0.492, 0),
        (-0.492, 0.492, 0),
        (12331.4, 12350.5, 12300),
        (-12335.3, 12335.3, 0),
    ]

    use_offset_data = [True, False]

    @pytest.mark.parametrize('left, right, offset', offset_data)
    def test_offset_value(self, left, right, offset):
        fig, ax = plt.subplots()
        formatter = ax.get_xaxis().get_major_formatter()

        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Attempting to set identical',
                                    UserWarning)
            ax.set_xlim(left, right)
        assert len(w) == (1 if left == right else 0)
        # Update ticks.
        next(ax.get_xaxis().iter_ticks())
        assert formatter.offset == offset

        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Attempting to set identical',
                                    UserWarning)
            ax.set_xlim(right, left)
        assert len(w) == (1 if left == right else 0)
        # Update ticks.
        next(ax.get_xaxis().iter_ticks())
        assert formatter.offset == offset

    @pytest.mark.parametrize('use_offset', use_offset_data)
    def test_use_offset(self, use_offset):
        with matplotlib.rc_context({'axes.formatter.useoffset': use_offset}):
            tmp_form = mticker.ScalarFormatter()
            assert use_offset == tmp_form.get_useOffset()


class FakeAxis(object):
    """Allow Formatter to be called without having a "full" plot set up."""
    def __init__(self, vmin=1, vmax=10):
        self.vmin = vmin
        self.vmax = vmax

    def get_view_interval(self):
        return self.vmin, self.vmax


class TestLogFormatterExponent(object):
    param_data = [
        (True, 4, np.arange(-3, 4.0), np.arange(-3, 4.0),
         ['-3', '-2', '-1', '0', '1', '2', '3']),
        # With labelOnlyBase=False, non-integer powers should be nicely
        # formatted.
        (False, 10, np.array([0.1, 0.00001, np.pi, 0.2, -0.2, -0.00001]),
         range(6), ['0.1', '1e-05', '3.14', '0.2', '-0.2', '-1e-05']),
        (False, 50, np.array([3, 5, 12, 42], dtype='float'), range(6),
         ['3', '5', '12', '42']),
    ]

    base_data = [2.0, 5.0, 10.0, np.pi, np.e]

    @pytest.mark.parametrize(
            'labelOnlyBase, exponent, locs, positions, expected', param_data)
    @pytest.mark.parametrize('base', base_data)
    def test_basic(self, labelOnlyBase, base, exponent, locs, positions,
                   expected):
        formatter = mticker.LogFormatterExponent(base=base,
                                                 labelOnlyBase=labelOnlyBase)
        formatter.axis = FakeAxis(1, base**exponent)
        vals = base**locs
        labels = [formatter(x, pos) for (x, pos) in zip(vals, positions)]
        assert labels == expected

    def test_blank(self):
        # Should be a blank string for non-integer powers if labelOnlyBase=True
        formatter = mticker.LogFormatterExponent(base=10, labelOnlyBase=True)
        formatter.axis = FakeAxis()
        assert formatter(10**0.1) == ''


class TestLogFormatterMathtext():
    fmt = mticker.LogFormatterMathtext()
    test_data = [
        (0, 1, '$\\mathdefault{10^{0}}$'),
        (0, 1e-2, '$\\mathdefault{10^{-2}}$'),
        (0, 1e2, '$\\mathdefault{10^{2}}$'),
        (3, 1, '$\\mathdefault{1}$'),
        (3, 1e-2, '$\\mathdefault{0.01}$'),
        (3, 1e2, '$\\mathdefault{100}$'),
        (3, 1e-3, '$\\mathdefault{10^{-3}}$'),
        (3, 1e3, '$\\mathdefault{10^{3}}$'),
    ]

    @pytest.mark.parametrize('min_exponent, value, expected', test_data)
    def test_min_exponent(self, min_exponent, value, expected):
        with matplotlib.rc_context({'axes.formatter.min_exponent':
                                    min_exponent}):
            assert self.fmt(value) == expected


class TestLogFormatterSciNotation(object):
    test_data = [
        (2, 0.03125, '$\\mathdefault{2^{-5}}$'),
        (2, 1, '$\\mathdefault{2^{0}}$'),
        (2, 32, '$\\mathdefault{2^{5}}$'),
        (2, 0.0375, '$\\mathdefault{1.2\\times2^{-5}}$'),
        (2, 1.2, '$\\mathdefault{1.2\\times2^{0}}$'),
        (2, 38.4, '$\\mathdefault{1.2\\times2^{5}}$'),
        (10, -1, '$\\mathdefault{-10^{0}}$'),
        (10, 1e-05, '$\\mathdefault{10^{-5}}$'),
        (10, 1, '$\\mathdefault{10^{0}}$'),
        (10, 100000, '$\\mathdefault{10^{5}}$'),
        (10, 2e-05, '$\\mathdefault{2\\times10^{-5}}$'),
        (10, 2, '$\\mathdefault{2\\times10^{0}}$'),
        (10, 200000, '$\\mathdefault{2\\times10^{5}}$'),
        (10, 5e-05, '$\\mathdefault{5\\times10^{-5}}$'),
        (10, 5, '$\\mathdefault{5\\times10^{0}}$'),
        (10, 500000, '$\\mathdefault{5\\times10^{5}}$'),
    ]

    @pytest.mark.style('default')
    @pytest.mark.parametrize('base, value, expected', test_data)
    def test_basic(self, base, value, expected):
        formatter = mticker.LogFormatterSciNotation(base=base)
        formatter.sublabel = {1, 2, 5, 1.2}
        with matplotlib.rc_context({'text.usetex': False}):
            assert formatter(value) == expected


class TestLogFormatter(object):
    pprint_data = [
        (3.141592654e-05, 0.001, '3.142e-5'),
        (0.0003141592654, 0.001, '3.142e-4'),
        (0.003141592654, 0.001, '3.142e-3'),
        (0.03141592654, 0.001, '3.142e-2'),
        (0.3141592654, 0.001, '3.142e-1'),
        (3.141592654, 0.001, '3.142'),
        (31.41592654, 0.001, '3.142e1'),
        (314.1592654, 0.001, '3.142e2'),
        (3141.592654, 0.001, '3.142e3'),
        (31415.92654, 0.001, '3.142e4'),
        (314159.2654, 0.001, '3.142e5'),
        (1e-05, 0.001, '1e-5'),
        (0.0001, 0.001, '1e-4'),
        (0.001, 0.001, '1e-3'),
        (0.01, 0.001, '1e-2'),
        (0.1, 0.001, '1e-1'),
        (1, 0.001, '1'),
        (10, 0.001, '10'),
        (100, 0.001, '100'),
        (1000, 0.001, '1000'),
        (10000, 0.001, '1e4'),
        (100000, 0.001, '1e5'),
        (3.141592654e-05, 0.015, '0'),
        (0.0003141592654, 0.015, '0'),
        (0.003141592654, 0.015, '0.003'),
        (0.03141592654, 0.015, '0.031'),
        (0.3141592654, 0.015, '0.314'),
        (3.141592654, 0.015, '3.142'),
        (31.41592654, 0.015, '31.416'),
        (314.1592654, 0.015, '314.159'),
        (3141.592654, 0.015, '3141.593'),
        (31415.92654, 0.015, '31415.927'),
        (314159.2654, 0.015, '314159.265'),
        (1e-05, 0.015, '0'),
        (0.0001, 0.015, '0'),
        (0.001, 0.015, '0.001'),
        (0.01, 0.015, '0.01'),
        (0.1, 0.015, '0.1'),
        (1, 0.015, '1'),
        (10, 0.015, '10'),
        (100, 0.015, '100'),
        (1000, 0.015, '1000'),
        (10000, 0.015, '10000'),
        (100000, 0.015, '100000'),
        (3.141592654e-05, 0.5, '0'),
        (0.0003141592654, 0.5, '0'),
        (0.003141592654, 0.5, '0.003'),
        (0.03141592654, 0.5, '0.031'),
        (0.3141592654, 0.5, '0.314'),
        (3.141592654, 0.5, '3.142'),
        (31.41592654, 0.5, '31.416'),
        (314.1592654, 0.5, '314.159'),
        (3141.592654, 0.5, '3141.593'),
        (31415.92654, 0.5, '31415.927'),
        (314159.2654, 0.5, '314159.265'),
        (1e-05, 0.5, '0'),
        (0.0001, 0.5, '0'),
        (0.001, 0.5, '0.001'),
        (0.01, 0.5, '0.01'),
        (0.1, 0.5, '0.1'),
        (1, 0.5, '1'),
        (10, 0.5, '10'),
        (100, 0.5, '100'),
        (1000, 0.5, '1000'),
        (10000, 0.5, '10000'),
        (100000, 0.5, '100000'),
        (3.141592654e-05, 5, '0'),
        (0.0003141592654, 5, '0'),
        (0.003141592654, 5, '0'),
        (0.03141592654, 5, '0.03'),
        (0.3141592654, 5, '0.31'),
        (3.141592654, 5, '3.14'),
        (31.41592654, 5, '31.42'),
        (314.1592654, 5, '314.16'),
        (3141.592654, 5, '3141.59'),
        (31415.92654, 5, '31415.93'),
        (314159.2654, 5, '314159.27'),
        (1e-05, 5, '0'),
        (0.0001, 5, '0'),
        (0.001, 5, '0'),
        (0.01, 5, '0.01'),
        (0.1, 5, '0.1'),
        (1, 5, '1'),
        (10, 5, '10'),
        (100, 5, '100'),
        (1000, 5, '1000'),
        (10000, 5, '10000'),
        (100000, 5, '100000'),
        (3.141592654e-05, 100, '0'),
        (0.0003141592654, 100, '0'),
        (0.003141592654, 100, '0'),
        (0.03141592654, 100, '0'),
        (0.3141592654, 100, '0.3'),
        (3.141592654, 100, '3.1'),
        (31.41592654, 100, '31.4'),
        (314.1592654, 100, '314.2'),
        (3141.592654, 100, '3141.6'),
        (31415.92654, 100, '31415.9'),
        (314159.2654, 100, '314159.3'),
        (1e-05, 100, '0'),
        (0.0001, 100, '0'),
        (0.001, 100, '0'),
        (0.01, 100, '0'),
        (0.1, 100, '0.1'),
        (1, 100, '1'),
        (10, 100, '10'),
        (100, 100, '100'),
        (1000, 100, '1000'),
        (10000, 100, '10000'),
        (100000, 100, '100000'),
        (3.141592654e-05, 1000000.0, '3.1e-5'),
        (0.0003141592654, 1000000.0, '3.1e-4'),
        (0.003141592654, 1000000.0, '3.1e-3'),
        (0.03141592654, 1000000.0, '3.1e-2'),
        (0.3141592654, 1000000.0, '3.1e-1'),
        (3.141592654, 1000000.0, '3.1'),
        (31.41592654, 1000000.0, '3.1e1'),
        (314.1592654, 1000000.0, '3.1e2'),
        (3141.592654, 1000000.0, '3.1e3'),
        (31415.92654, 1000000.0, '3.1e4'),
        (314159.2654, 1000000.0, '3.1e5'),
        (1e-05, 1000000.0, '1e-5'),
        (0.0001, 1000000.0, '1e-4'),
        (0.001, 1000000.0, '1e-3'),
        (0.01, 1000000.0, '1e-2'),
        (0.1, 1000000.0, '1e-1'),
        (1, 1000000.0, '1'),
        (10, 1000000.0, '10'),
        (100, 1000000.0, '100'),
        (1000, 1000000.0, '1000'),
        (10000, 1000000.0, '1e4'),
        (100000, 1000000.0, '1e5'),
    ]

    @pytest.mark.parametrize('value, domain, expected', pprint_data)
    def test_pprint(self, value, domain, expected):
        fmt = mticker.LogFormatter()
        label = fmt.pprint_val(value, domain)
        assert label == expected

    def _sub_labels(self, axis, subs=()):
        "Test whether locator marks subs to be labeled"
        fmt = axis.get_minor_formatter()
        minor_tlocs = axis.get_minorticklocs()
        fmt.set_locs(minor_tlocs)
        coefs = minor_tlocs / 10**(np.floor(np.log10(minor_tlocs)))
        label_expected = [np.round(c) in subs for c in coefs]
        label_test = [fmt(x) != '' for x in minor_tlocs]
        assert label_test == label_expected

    @pytest.mark.style('default')
    def test_sublabel(self):
        # test label locator
        fig, ax = plt.subplots()
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mticker.LogLocator(base=10, subs=[]))
        ax.xaxis.set_minor_locator(mticker.LogLocator(base=10,
                                                      subs=np.arange(2, 10)))
        ax.xaxis.set_major_formatter(mticker.LogFormatter(labelOnlyBase=True))
        ax.xaxis.set_minor_formatter(mticker.LogFormatter(labelOnlyBase=False))
        # axis range above 3 decades, only bases are labeled
        ax.set_xlim(1, 1e4)
        fmt = ax.xaxis.get_major_formatter()
        fmt.set_locs(ax.xaxis.get_majorticklocs())
        show_major_labels = [fmt(x) != ''
                             for x in ax.xaxis.get_majorticklocs()]
        assert np.all(show_major_labels)
        self._sub_labels(ax.xaxis, subs=[])

        # For the next two, if the numdec threshold in LogFormatter.set_locs
        # were 3, then the label sub would be 3 for 2-3 decades and (2,5)
        # for 1-2 decades.  With a threshold of 1, subs are not labeled.
        # axis range at 2 to 3 decades
        ax.set_xlim(1, 800)
        self._sub_labels(ax.xaxis, subs=[])

        # axis range at 1 to 2 decades
        ax.set_xlim(1, 80)
        self._sub_labels(ax.xaxis, subs=[])

        # axis range at 0.4 to 1 decades, label subs 2, 3, 4, 6
        ax.set_xlim(1, 8)
        self._sub_labels(ax.xaxis, subs=[2, 3, 4, 6])

        # axis range at 0 to 0.4 decades, label all
        ax.set_xlim(0.5, 0.9)
        self._sub_labels(ax.xaxis, subs=np.arange(2, 10, dtype=int))

    @pytest.mark.parametrize('val', [1, 10, 100, 1000])
    def test_LogFormatter_call(self, val):
        # test _num_to_string method used in __call__
        temp_lf = mticker.LogFormatter()
        temp_lf.axis = FakeAxis()
        assert temp_lf(val) == str(val)


class TestFormatStrFormatter(object):
    def test_basic(self):
        # test % style formatter
        tmp_form = mticker.FormatStrFormatter('%05d')
        assert '00002' == tmp_form(2)


class TestStrMethodFormatter(object):
    test_data = [
        ('{x:05d}', (2,), '00002'),
        ('{x:03d}-{pos:02d}', (2, 1), '002-01'),
    ]

    @pytest.mark.parametrize('format, input, expected', test_data)
    def test_basic(self, format, input, expected):
        fmt = mticker.StrMethodFormatter(format)
        assert fmt(*input) == expected


class TestEngFormatter(object):
    # (input, expected) where ''expected'' corresponds to the outputs
    # respectively returned when (places=None, places=0, places=2)
    raw_format_data = [
        (-1234.56789, ('-1.23457 k', '-1 k', '-1.23 k')),
        (-1.23456789, ('-1.23457', '-1', '-1.23')),
        (-0.123456789, ('-123.457 m', '-123 m', '-123.46 m')),
        (-0.00123456789, ('-1.23457 m', '-1 m', '-1.23 m')),
        (-0.0, ('0', '0', '0.00')),
        (-0, ('0', '0', '0.00')),
        (0, ('0', '0', '0.00')),
        (1.23456789e-6, (u'1.23457 \u03bc', u'1 \u03bc', u'1.23 \u03bc')),
        (0.123456789, ('123.457 m', '123 m', '123.46 m')),
        (0.1, ('100 m', '100 m', '100.00 m')),
        (1, ('1', '1', '1.00')),
        (1.23456789, ('1.23457', '1', '1.23')),
        (999.9, ('999.9', '1 k', '999.90')),  # places=0: corner-case rounding
        (999.9999, ('1 k', '1 k', '1.00 k')),  # corner-case roudning for all
        (1000, ('1 k', '1 k', '1.00 k')),
        (1001, ('1.001 k', '1 k', '1.00 k')),
        (100001, ('100.001 k', '100 k', '100.00 k')),
        (987654.321, ('987.654 k', '988 k', '987.65 k')),
        (1.23e27, ('1230 Y', '1230 Y', '1230.00 Y'))  # OoR value (> 1000 Y)
    ]

    @pytest.mark.parametrize('input, expected', raw_format_data)
    def test_params(self, input, expected):
        """
        Test the formatting of EngFormatter for various values of the 'places'
        argument, in several cases:
            0. without a unit symbol but with a (default) space separator;
            1. with both a unit symbol and a (default) space separator;
            2. with both a unit symbol and some non default separators;
            3. without a unit symbol but with some non default separators.
        Note that cases 2. and 3. are looped over several separator strings.
        """

        UNIT = 's'  # seconds
        DIGITS = '0123456789'  # %timeit showed 10-20% faster search than set

        # Case 0: unit='' (default) and sep=' ' (default).
        # 'expected' already corresponds to this reference case.
        exp_outputs = expected
        formatters = (
            mticker.EngFormatter(),  # places=None (default)
            mticker.EngFormatter(places=0),
            mticker.EngFormatter(places=2)
        )
        for _formatter, _exp_output in zip(formatters, exp_outputs):
            assert _formatter(input) == _exp_output

        # Case 1: unit=UNIT and sep=' ' (default).
        # Append a unit symbol to the reference case.
        # Beware of the values in [1, 1000), where there is no prefix!
        exp_outputs = (_s + " " + UNIT if _s[-1] in DIGITS  # case w/o prefix
                       else _s + UNIT for _s in expected)
        formatters = (
            mticker.EngFormatter(unit=UNIT),  # places=None (default)
            mticker.EngFormatter(unit=UNIT, places=0),
            mticker.EngFormatter(unit=UNIT, places=2)
        )
        for _formatter, _exp_output in zip(formatters, exp_outputs):
            assert _formatter(input) == _exp_output

        # Test several non default separators: no separator, a narrow
        # no-break space (unicode character) and an extravagant string.
        for _sep in ("", "\N{NARROW NO-BREAK SPACE}", "@_@"):
            # Case 2: unit=UNIT and sep=_sep.
            # Replace the default space separator from the reference case
            # with the tested one `_sep` and append a unit symbol to it.
            exp_outputs = (_s + _sep + UNIT if _s[-1] in DIGITS  # no prefix
                           else _s.replace(" ", _sep) + UNIT
                           for _s in expected)
            formatters = (
                mticker.EngFormatter(unit=UNIT, sep=_sep),  # places=None
                mticker.EngFormatter(unit=UNIT, places=0, sep=_sep),
                mticker.EngFormatter(unit=UNIT, places=2, sep=_sep)
            )
            for _formatter, _exp_output in zip(formatters, exp_outputs):
                assert _formatter(input) == _exp_output

            # Case 3: unit='' (default) and sep=_sep.
            # Replace the default space separator from the reference case
            # with the tested one `_sep`. Reference case is already unitless.
            exp_outputs = (_s.replace(" ", _sep) for _s in expected)
            formatters = (
                mticker.EngFormatter(sep=_sep),  # places=None (default)
                mticker.EngFormatter(places=0, sep=_sep),
                mticker.EngFormatter(places=2, sep=_sep)
            )
            for _formatter, _exp_output in zip(formatters, exp_outputs):
                assert _formatter(input) == _exp_output


class TestPercentFormatter(object):
    percent_data = [
        # Check explicitly set decimals over different intervals and values
        (100, 0, '%', 120, 100, '120%'),
        (100, 0, '%', 100, 90, '100%'),
        (100, 0, '%', 90, 50, '90%'),
        (100, 0, '%', -1.7, 40, '-2%'),
        (100, 1, '%', 90.0, 100, '90.0%'),
        (100, 1, '%', 80.1, 90, '80.1%'),
        (100, 1, '%', 70.23, 50, '70.2%'),
        # 60.554 instead of 60.55: see https://bugs.python.org/issue5118
        (100, 1, '%', -60.554, 40, '-60.6%'),
        # Check auto decimals over different intervals and values
        (100, None, '%', 95, 1, '95.00%'),
        (1.0, None, '%', 3, 6, '300%'),
        (17.0, None, '%', 1, 8.5, '6%'),
        (17.0, None, '%', 1, 8.4, '5.9%'),
        (5, None, '%', -100, 0.000001, '-2000.00000%'),
        # Check percent symbol
        (1.0, 2, None, 1.2, 100, '120.00'),
        (75, 3, '', 50, 100, '66.667'),
        (42, None, '^^Foobar$$', 21, 12, '50.0^^Foobar$$'),
    ]

    percent_ids = [
        # Check explicitly set decimals over different intervals and values
        'decimals=0, x>100%',
        'decimals=0, x=100%',
        'decimals=0, x<100%',
        'decimals=0, x<0%',
        'decimals=1, x>100%',
        'decimals=1, x=100%',
        'decimals=1, x<100%',
        'decimals=1, x<0%',
        # Check auto decimals over different intervals and values
        'autodecimal, x<100%, display_range=1',
        'autodecimal, x>100%, display_range=6 (custom xmax test)',
        'autodecimal, x<100%, display_range=8.5 (autodecimal test 1)',
        'autodecimal, x<100%, display_range=8.4 (autodecimal test 2)',
        'autodecimal, x<-100%, display_range=1e-6 (tiny display range)',
        # Check percent symbol
        'None as percent symbol',
        'Empty percent symbol',
        'Custom percent symbol',
    ]

    latex_data = [
        (False, False, r'50\{t}%'),
        (False, True, r'50\\\{t\}\%'),
        (True, False, r'50\{t}%'),
        (True, True, r'50\{t}%'),
    ]

    @pytest.mark.parametrize(
            'xmax, decimals, symbol, x, display_range, expected',
            percent_data, ids=percent_ids)
    def test_basic(self, xmax, decimals, symbol,
                   x, display_range, expected):
        formatter = mticker.PercentFormatter(xmax, decimals, symbol)
        with matplotlib.rc_context(rc={'text.usetex': False}):
            assert formatter.format_pct(x, display_range) == expected

    @pytest.mark.parametrize('is_latex, usetex, expected', latex_data)
    def test_latex(self, is_latex, usetex, expected):
        fmt = mticker.PercentFormatter(symbol='\\{t}%', is_latex=is_latex)
        with matplotlib.rc_context(rc={'text.usetex': usetex}):
            assert fmt.format_pct(50, 100) == expected
