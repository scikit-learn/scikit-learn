from __future__ import division, print_function, absolute_import

import itertools

import numpy as np
from numpy.testing import assert_
from scipy.special._testutils import FuncData
import pytest

from scipy.special import smirnov, smirnovi, kolmogorov, kolmogi

_rtol = 1e-10

class TestSmirnov(object):
    def test_nan(self):
        assert_(np.isnan(smirnov(1, np.nan)))

    def test_basic(self):
        dataset = [(1, 0.1, 0.9),
                   (1, 0.875, 0.125),
                   (2, 0.875, 0.125 * 0.125),
                   (3, 0.875, 0.125 * 0.125 * 0.125)]

        dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_x_equals_0(self):
        dataset = [(n, 0, 1) for n in itertools.chain(range(2, 20), range(1010, 1020))]
        dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_x_equals_1(self):
        dataset = [(n, 1, 0) for n in itertools.chain(range(2, 20), range(1010, 1020))]
        dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_x_equals_0point5(self):
        dataset = [(1, 0.5, 0.5),
                   (2, 0.5, 0.25),
                   (3, 0.5, 0.166666666667),
                   (4, 0.5, 0.09375),
                   (5, 0.5, 0.056),
                   (6, 0.5, 0.0327932098765),
                   (7, 0.5, 0.0191958707681),
                   (8, 0.5, 0.0112953186035),
                   (9, 0.5, 0.00661933257355),
                   (10, 0.5, 0.003888705)]

        dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_n_equals_1(self):
        x = np.linspace(0, 1, 101, endpoint=True)
        dataset = np.column_stack([[1]*len(x), x, 1-x])
        # dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_n_equals_2(self):
        x = np.linspace(0.5, 1, 101, endpoint=True)
        p = np.power(1-x, 2)
        n = np.array([2] * len(x))
        dataset = np.column_stack([n, x, p])
        # dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_n_equals_3(self):
        x = np.linspace(0.7, 1, 31, endpoint=True)
        p = np.power(1-x, 3)
        n = np.array([3] * len(x))
        dataset = np.column_stack([n, x, p])
        # dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_n_large(self):
        # test for large values of n
        # Probabilities should go down as n goes up
        x = 0.4
        pvals = np.array([smirnov(n, x) for n in range(400, 1100, 20)])
        dfs = np.diff(pvals)
        assert_(np.all(dfs <= 0), msg='Not all diffs negative %s' % dfs)

        dataset = [(1000, 1 - 1.0/2000, np.power(2000.0, -1000))]
        dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=_rtol).check()

        # Check asymptotic behaviour
        dataset = [(n, 1.0 / np.sqrt(n), np.exp(-2)) for n in range(1000, 5000, 1000)]
        dataset = np.asarray(dataset)
        FuncData(smirnov, dataset, (0, 1), 2, rtol=.05).check()


class TestSmirnovi(object):
    def test_nan(self):
        assert_(np.isnan(smirnovi(1, np.nan)))

    @pytest.mark.xfail(reason="test fails; smirnovi() is not always accurate")
    def test_basic(self):
        dataset = [(1, 0.4, 0.6),
                   (1, 0.6, 0.4),
                   (1, 0.99, 0.01),
                   (1, 0.01, 0.99),
                   (2, 0.125 * 0.125, 0.875),
                   (3, 0.125 * 0.125 * 0.125, 0.875),
                   (10, 1.0 / 16 ** 10, 1 - 1.0 / 16)]

        dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; smirnovi(_,0) is not accurate")
    def test_x_equals_0(self):
        dataset = [(n, 0, 1) for n in itertools.chain(range(2, 20), range(1010, 1020))]
        dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_x_equals_1(self):
        dataset = [(n, 1, 0) for n in itertools.chain(range(2, 20), range(1010, 1020))]
        dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; smirnovi(1,) is not accurate")
    def test_n_equals_1(self):
        pp = np.linspace(0, 1, 101, endpoint=True)
        dataset = [(1, p, 1-p) for p in pp]
        dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; smirnovi(2,_) is not accurate")
    def test_n_equals_2(self):
        x = np.linspace(0.5, 1, 101, endpoint=True)
        p = np.power(1-x, 2)
        n = np.array([2] * len(x))
        dataset = np.column_stack([n, p, x])
        # dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; smirnovi(3,_) is not accurate")
    def test_n_equals_3(self):
        x = np.linspace(0.7, 1, 31, endpoint=True)
        p = np.power(1-x, 3)
        n = np.array([3] * len(x))
        dataset = np.column_stack([n, p, x])
        # dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; smirnovi(_,_) is not accurate")
    def test_round_trip(self):
        def _sm_smi(n, p):
            return smirnov(n, smirnovi(n, p))

        dataset = [(1, 0.4, 0.4),
                   (1, 0.6, 0.6),
                   (2, 0.875, 0.875),
                   (3, 0.875, 0.875),
                   (3, 0.125, 0.125),
                   (10, 0.999, 0.999),
                   (10, 0.0001, 0.0001)]

        dataset = np.asarray(dataset)
        FuncData(_sm_smi, dataset, (0, 1), 2, rtol=_rtol).check()

    def test_x_equals_0point5(self):
        dataset = [(1, 0.5, 0.5),
                   (2, 0.5, 0.366025403784),
                   (2, 0.25, 0.5),
                   (3, 0.5, 0.297156508177),
                   (4, 0.5, 0.255520481121),
                   (5, 0.5, 0.234559536069),
                   (6, 0.5, 0.21715965898),
                   (7, 0.5, 0.202722580034),
                   (8, 0.5, 0.190621765256),
                   (9, 0.5, 0.180363501362),
                   (10, 0.5, 0.17157867006)]

        dataset = np.asarray(dataset)
        FuncData(smirnovi, dataset, (0, 1), 2, rtol=_rtol).check()


class TestKolmogorov(object):
    def test_nan(self):
        assert_(np.isnan(kolmogorov(np.nan)))

    def test_basic(self):
        dataset = [(0, 1.0),
                   (0.5, 0.96394524366487511),
                   (1, 0.26999967167735456),
                   (2, 0.00067092525577969533)]

        dataset = np.asarray(dataset)
        FuncData(kolmogorov, dataset, (0,), 1, rtol=_rtol).check()

    def test_smallx(self):
        epsilon = 0.1 ** np.arange(1, 14)
        x = np.array([0.571173265106, 0.441027698518, 0.374219690278, 0.331392659217,
                      0.300820537459, 0.277539353999, 0.259023494805, 0.243829561254,
                      0.231063086389, 0.220135543236, 0.210641372041, 0.202290283658,
                      0.19487060742])

        dataset = np.column_stack([x, 1-epsilon])
        FuncData(kolmogorov, dataset, (0,), 1, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; kolmogi() is not accurate for small p")
    def test_round_trip(self):
        def _ki_k(_x):
            return kolmogi(kolmogorov(_x))

        x = np.linspace(0.0, 2.0, 21, endpoint=True)
        dataset = np.column_stack([x, x])
        FuncData(_ki_k, dataset, (0,), 1, rtol=_rtol).check()


class TestKolmogi(object):
    def test_nan(self):
        assert_(np.isnan(kolmogi(np.nan)))

    @pytest.mark.xfail(reason="test fails; kolmogi() is not accurate for small p")
    def test_basic(self):
        dataset = [(1.0, 0),
                   (0.96394524366487511, 0.5),
                   (0.26999967167735456, 1),
                   (0.00067092525577969533, 2)]

        dataset = np.asarray(dataset)
        FuncData(kolmogi, dataset, (0,), 1, rtol=_rtol).check()

    @pytest.mark.xfail(reason="test fails; kolmogi() is not accurate for small p")
    def test_smallp(self):
        epsilon = 0.1 ** np.arange(1, 14)
        x = np.array([0.571173265106, 0.441027698518, 0.374219690278, 0.331392659217,
                      0.300820537459, 0.277539353999, 0.259023494805, 0.243829561254,
                      0.231063086389, 0.220135543236, 0.210641372041, 0.202290283658,
                      0.19487060742])

        dataset = np.column_stack([1-epsilon, x])
        FuncData(kolmogi, dataset, (0,), 1, rtol=_rtol).check()

    def test_round_trip(self):
        def _k_ki(_p):
            return kolmogorov(kolmogi(_p))

        p = np.linspace(0.1, 1.0, 10, endpoint=True)
        dataset = np.column_stack([p, p])
        FuncData(_k_ki, dataset, (0,), 1, rtol=_rtol).check()
