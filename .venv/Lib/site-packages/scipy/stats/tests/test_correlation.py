import pytest
import numpy as np

from scipy._lib._array_api import make_xp_test_case, xp_default_dtype, is_jax
from scipy._lib._array_api_no_0d import xp_assert_close
from scipy import stats
from scipy.stats._axis_nan_policy import SmallSampleWarning


@make_xp_test_case(stats.chatterjeexi)
class TestChatterjeeXi:
    @pytest.mark.parametrize('case', [
        dict(y_cont=True, statistic=-0.303030303030303, pvalue=0.9351329808526656),
        dict(y_cont=False, statistic=0.07407407407407396, pvalue=0.3709859367123997)])
    @pytest.mark.parametrize('dtype', ['float32', 'float64', None])
    def test_against_R_XICOR(self, case, dtype, xp):
        # Test against R package XICOR, e.g.
        # library(XICOR)
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558,
        #       0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495)
        # y = c(0.8402081904493103, 0.5946972833914318, 0.23481606164114155,
        #       0.49754786197715384, 0.9146460831206026, 0.5848057749217579,
        #       0.7620801065573549, 0.31410063302647495, 0.7935620302236199,
        #       0.5423085761365468)
        # xicor(x, y, ties=FALSE, pvalue=TRUE)
        if dtype == 'float32' and np.__version__ < "2":
            pytest.skip("Scalar dtypes only respected after NEP 50.")
        dtype = xp_default_dtype(xp) if dtype is None else getattr(xp, dtype)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=10)
        y = (rng.random(size=10) if case['y_cont']
             else rng.integers(0, 5, size=10))

        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.chatterjeexi(x, y, y_continuous=case['y_cont'])

        xp_assert_close(res.statistic, xp.asarray(case['statistic'], dtype=dtype))
        xp_assert_close(res.pvalue, xp.asarray(case['pvalue'], dtype=dtype))

    @pytest.mark.parametrize('y_continuous', (False, True))
    def test_permutation_asymptotic(self, y_continuous):
        # XICOR doesn't seem to perform the permutation test as advertised, so
        # compare the result of a permutation test against an asymptotic test.
        rng = np.random.default_rng(2524579827426)
        n = np.floor(rng.uniform(100, 150)).astype(int)
        shape = (2, n)
        x = rng.random(size=shape)
        y = (rng.random(size=shape) if y_continuous
             else rng.integers(0, 10, size=shape))
        method = stats.PermutationMethod(rng=rng)
        res = stats.chatterjeexi(x, y, method=method,
                                 y_continuous=y_continuous, axis=-1)
        ref = stats.chatterjeexi(x, y, y_continuous=y_continuous, axis=-1)
        np.testing.assert_allclose(res.statistic, ref.statistic, rtol=1e-15)
        np.testing.assert_allclose(res.pvalue, ref.pvalue, rtol=2e-2)

    def test_input_validation(self, xp):
        rng = np.random.default_rng(25932435798274926)
        x, y = rng.random(size=(2, 10))
        x, y = xp.asarray(x), xp.asarray(y)

        message = 'Array shapes are incompatible for broadcasting.|Incompatible shapes'
        with pytest.raises((ValueError, TypeError), match=message):
            stats.chatterjeexi(x, y[:-1])

        if not is_jax(xp):
            # jax misses out on some input validation from _axis_nan_policy decorator
            message = '...axis 10 is out of bounds for array...|out of range'
            with pytest.raises((ValueError, IndexError), match=message):
                stats.chatterjeexi(x, y, axis=10)

        message = '`y_continuous` must be boolean.'
        with pytest.raises(ValueError, match=message):
            stats.chatterjeexi(x, y, y_continuous='a herring')

        message = "`method` must be 'asymptotic' or"
        with pytest.raises(ValueError, match=message):
            stats.chatterjeexi(x, y, method='ekki ekii')

    @pytest.mark.skip_xp_backends('jax.numpy', reason='no SmallSampleWarning (lazy)')
    def test_special_cases(self, xp):
        message = 'One or more sample arguments is too small...'
        with pytest.warns(SmallSampleWarning, match=message):
            res = stats.chatterjeexi(xp.asarray([1]), xp.asarray([2]))

        assert xp.isnan(res.statistic)
        assert xp.isnan(res.pvalue)


@make_xp_test_case(stats.spearmanrho)
class TestSpearmanRho:
    @pytest.mark.parametrize('alternative, statistic, pvalue', [
        ('two-sided', -0.2727272727272727, 0.4458383415428),
        ('greater', -0.2727272727272727, 0.7770808292286),
        ('less', -0.2727272727272727, 0.2229191707714),
    ])
    @pytest.mark.parametrize('dtype', ['float32', 'float64', None])
    def test_against_R_cor_test(self, alternative, statistic, pvalue, dtype, xp):
        # Test against R cor.test, e.g.
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558,
        #       0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495)
        # y = c(0.8402081904493103, 0.5946972833914318, 0.23481606164114155,
        #       0.49754786197715384, 0.9146460831206026, 0.5848057749217579,
        #       0.7620801065573549, 0.31410063302647495, 0.7935620302236199,
        #       0.5423085761365468)
        # cor.test(x, y, method='spearman', alternative='t', exact=FALSE)
        if dtype == 'float32' and np.__version__ < "2":
            pytest.skip("Scalar dtypes only respected after NEP 50.")
        dtype = xp_default_dtype(xp) if dtype is None else getattr(xp, dtype)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=10)
        y = rng.random(size=10)

        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.spearmanrho(x, y, alternative=alternative)

        xp_assert_close(res.statistic, xp.asarray(statistic, dtype=dtype))
        xp_assert_close(res.pvalue, xp.asarray(pvalue, dtype=dtype))


    @pytest.mark.parametrize('alternative, statistic, pvalue', [
        ('two-sided', -0.4857142857142857, 0.3555555555556),
        ('greater', -0.4857142857142857, 0.8513888888889),
        ('less', -0.4857142857142857, 0.1777777777778),
    ])
    def test_against_R_cor_test_exact(self, alternative, statistic, pvalue):
        xp = np  # will test for multiple backends when gh-23772 merges
        # Test against R cor.test exact=TRUE, e.g.
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558)
        # y = c(0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495, 0.8402081904493103, 0.5946972833914318)
        # cor.test(x, y, method='spearman', alternative='t', exact=TRUE)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=6)
        y = rng.random(size=6)
        x, y = xp.asarray(x), xp.asarray(y)

        method = stats.PermutationMethod()
        res = stats.spearmanrho(x, y, method=method, alternative=alternative)

        xp_assert_close(res.statistic, xp.asarray(statistic))
        xp_assert_close(res.pvalue, xp.asarray(pvalue))


    @pytest.mark.parametrize('alternative', ('two-sided', 'greater', 'less'))
    @pytest.mark.parametrize('n', [9, 99, 999])
    def test_against_scipy_spearmanr(self, alternative, n, xp):
        rng = np.random.default_rng(5982435982346983)
        x = rng.integers(n//2, size=n)
        y = rng.integers(n//2, size=n)
        dtype = xp_default_dtype(xp)

        ref = stats.spearmanr(x, y, alternative=alternative)
        ref_statistic = xp.asarray(ref.statistic, dtype=dtype)
        ref_pvalue = xp.asarray(ref.pvalue, dtype=dtype)

        x = xp.asarray(x, dtype=dtype)
        y = xp.asarray(y, dtype=dtype)
        res = stats.spearmanrho(x, y, alternative=alternative)

        xp_assert_close(res.statistic, ref_statistic)
        xp_assert_close(res.pvalue, ref_pvalue)


    @pytest.mark.parametrize('axis', [-1, 0, 1])
    def test_other_backends_nd(self, axis, xp):
        # NumPy n-d behavior is tested by `test_axis_nan_policy`;
        # check other backends against it.
        rng = np.random.default_rng(8243598346983259)
        shape = (8, 9, 10)
        x = rng.standard_normal(size=shape)
        y = rng.standard_normal(size=shape)
        res = stats.spearmanrho(xp.asarray(x), xp.asarray(y), axis=axis)
        ref = stats.spearmanrho(x, y, axis=axis)
        xp_assert_close(res.statistic, xp.asarray(ref.statistic), atol=1e-16)
        xp_assert_close(res.pvalue, xp.asarray(ref.pvalue), atol=1e-16)


    def test_input_validation(self, xp):
        rng = np.random.default_rng(25932435798274926)
        x, y = rng.random(size=(2, 10))
        x, y = xp.asarray(x), xp.asarray(y)

        msg = 'incompatible for broadcasting|Incompatible shapes|must be broadcastable'
        with pytest.raises((ValueError, TypeError), match=msg):
            stats.spearmanrho(x, y[:-1])

        if not is_jax(xp):
            # jax misses out on some input validation from _axis_nan_policy decorator
            message = '...axis 10 is out of bounds for array...|out of range'
            with pytest.raises((ValueError, IndexError), match=message):
                stats.spearmanrho(x, y, axis=10)

        message = "`alternative` must be 'less', 'greater', or 'two-sided'."
        with pytest.raises(ValueError, match=message):
            stats.spearmanrho(x, y, alternative='alternative')

        message = ("`method` must be...")
        with pytest.raises(ValueError, match=message):
            stats.spearmanrho(x, y, method='method')


    @pytest.mark.skip_xp_backends('jax.numpy', reason='no SmallSampleWarning (lazy)')
    def test_special_cases(self, xp):
        def check_nan(res):
            assert xp.isnan(res.statistic)
            assert xp.isnan(res.pvalue)

        message = 'One or more sample arguments is too small...'
        with pytest.warns(SmallSampleWarning, match=message):
            res = stats.spearmanrho(xp.asarray([1]), xp.asarray([2]))
            check_nan(res)

        x = xp.asarray([1, 1, 1, 1, 1])
        y = xp.asarray([1, 2, 3, 4, 5])
        message = 'An input array is constant; the correlation coefficient...'
        with pytest.warns(stats.ConstantInputWarning, match=message):
            res = stats.spearmanrho(x, y)
            check_nan(res)
        with pytest.warns(stats.ConstantInputWarning, match=message):
            res = stats.spearmanrho(y, x)
            check_nan(res)
