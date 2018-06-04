# -*- coding: utf-8 -*-
""" Test printing of scalar types.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import assert_, assert_equal, run_module_suite
import sys, tempfile


class TestRealScalars(object):
    def test_str(self):
        svals = [0.0, -0.0, 1, -1, np.inf, -np.inf, np.nan]
        styps = [np.float16, np.float32, np.float64, np.longdouble]
        wanted = [
             ['0.0',  '0.0',  '0.0',  '0.0' ],
             ['-0.0', '-0.0', '-0.0', '-0.0'],
             ['1.0',  '1.0',  '1.0',  '1.0' ],
             ['-1.0', '-1.0', '-1.0', '-1.0'],
             ['inf',  'inf',  'inf',  'inf' ],
             ['-inf', '-inf', '-inf', '-inf'],
             ['nan',  'nan',  'nan',  'nan']]

        for wants, val in zip(wanted, svals):
            for want, styp in zip(wants, styps):
                msg = 'for str({}({}))'.format(np.dtype(styp).name, repr(val))
                assert_equal(str(styp(val)), want, err_msg=msg)

    def test_scalar_cutoffs(self):
        # test that both the str and repr of np.float64 behaves
        # like python floats in python3. Note that in python2
        # the str has truncated digits, but we do not do this
        def check(v):
            # we compare str to repr, to avoid python2 truncation behavior
            assert_equal(str(np.float64(v)), repr(v))
            assert_equal(repr(np.float64(v)), repr(v))

        # check we use the same number of significant digits
        check(1.12345678901234567890)
        check(0.0112345678901234567890)

        # check switch from scientific output to positional and back
        check(1e-5)
        check(1e-4)
        check(1e15)
        check(1e16)

    def test_py2_float_print(self):
        # gh-10753
        # In python2, the python float type implements an obsolte method
        # tp_print, which overrides tp_repr and tp_str when using the "print"
        # keyword/method to output to a "real file" (ie, not a StringIO). Make
        # sure we don't inherit it.
        x = np.double(0.1999999999999)
        f = tempfile.TemporaryFile('r+t') # must output to real file, not StringIO
        print(x, file=f)
        f.seek(0)
        output = f.read()
        f.close()
        assert_equal(output, '0.1999999999999\n')
        # (compare to "print 0.1999999999999" printing "0.2" in python2)

    def test_dragon4(self):
        # these tests are adapted from Ryan Juckett's dragon4 implementation,
        # see dragon4.c for details.

        fpos32 = lambda x, **k: np.format_float_positional(np.float32(x), **k)
        fsci32 = lambda x, **k: np.format_float_scientific(np.float32(x), **k)
        fpos64 = lambda x, **k: np.format_float_positional(np.float64(x), **k)
        fsci64 = lambda x, **k: np.format_float_scientific(np.float64(x), **k)

        preckwd = lambda prec: {'unique': False, 'precision': prec}

        assert_equal(fpos32('1.0'), "1.")
        assert_equal(fsci32('1.0'), "1.e+00")
        assert_equal(fpos32('10.234'), "10.234")
        assert_equal(fpos32('-10.234'), "-10.234")
        assert_equal(fsci32('10.234'), "1.0234e+01")
        assert_equal(fsci32('-10.234'), "-1.0234e+01")
        assert_equal(fpos32('1000.0'), "1000.")
        assert_equal(fpos32('1.0', precision=0), "1.")
        assert_equal(fsci32('1.0', precision=0), "1.e+00")
        assert_equal(fpos32('10.234', precision=0), "10.")
        assert_equal(fpos32('-10.234', precision=0), "-10.")
        assert_equal(fsci32('10.234', precision=0), "1.e+01")
        assert_equal(fsci32('-10.234', precision=0), "-1.e+01")
        assert_equal(fpos32('10.234', precision=2), "10.23")
        assert_equal(fsci32('-10.234', precision=2), "-1.02e+01")
        assert_equal(fsci64('9.9999999999999995e-08', **preckwd(16)),
                            '9.9999999999999995e-08')
        assert_equal(fsci64('9.8813129168249309e-324', **preckwd(16)),
                            '9.8813129168249309e-324')
        assert_equal(fsci64('9.9999999999999694e-311', **preckwd(16)),
                            '9.9999999999999694e-311')


        # test rounding
        # 3.1415927410 is closest float32 to np.pi
        assert_equal(fpos32('3.14159265358979323846', **preckwd(10)),
                            "3.1415927410")
        assert_equal(fsci32('3.14159265358979323846', **preckwd(10)),
                            "3.1415927410e+00")
        assert_equal(fpos64('3.14159265358979323846', **preckwd(10)),
                            "3.1415926536")
        assert_equal(fsci64('3.14159265358979323846', **preckwd(10)),
                            "3.1415926536e+00")
        # 299792448 is closest float32 to 299792458
        assert_equal(fpos32('299792458.0', **preckwd(5)), "299792448.00000")
        assert_equal(fsci32('299792458.0', **preckwd(5)), "2.99792e+08")
        assert_equal(fpos64('299792458.0', **preckwd(5)), "299792458.00000")
        assert_equal(fsci64('299792458.0', **preckwd(5)), "2.99792e+08")

        assert_equal(fpos32('3.14159265358979323846', **preckwd(25)),
                            "3.1415927410125732421875000")
        assert_equal(fpos64('3.14159265358979323846', **preckwd(50)),
                         "3.14159265358979311599796346854418516159057617187500")
        assert_equal(fpos64('3.14159265358979323846'), "3.141592653589793")


        # smallest numbers
        assert_equal(fpos32(0.5**(126 + 23), unique=False, precision=149),
                    "0.00000000000000000000000000000000000000000000140129846432"
                    "4817070923729583289916131280261941876515771757068283889791"
                    "08268586060148663818836212158203125")
        assert_equal(fpos64(0.5**(1022 + 52), unique=False, precision=1074),
                    "0.00000000000000000000000000000000000000000000000000000000"
                    "0000000000000000000000000000000000000000000000000000000000"
                    "0000000000000000000000000000000000000000000000000000000000"
                    "0000000000000000000000000000000000000000000000000000000000"
                    "0000000000000000000000000000000000000000000000000000000000"
                    "0000000000000000000000000000000000049406564584124654417656"
                    "8792868221372365059802614324764425585682500675507270208751"
                    "8652998363616359923797965646954457177309266567103559397963"
                    "9877479601078187812630071319031140452784581716784898210368"
                    "8718636056998730723050006387409153564984387312473397273169"
                    "6151400317153853980741262385655911710266585566867681870395"
                    "6031062493194527159149245532930545654440112748012970999954"
                    "1931989409080416563324524757147869014726780159355238611550"
                    "1348035264934720193790268107107491703332226844753335720832"
                    "4319360923828934583680601060115061698097530783422773183292"
                    "4790498252473077637592724787465608477820373446969953364701"
                    "7972677717585125660551199131504891101451037862738167250955"
                    "8373897335989936648099411642057026370902792427675445652290"
                    "87538682506419718265533447265625")

        # largest numbers
        assert_equal(fpos32(np.finfo(np.float32).max, **preckwd(0)),
                    "340282346638528859811704183484516925440.")
        assert_equal(fpos64(np.finfo(np.float64).max, **preckwd(0)),
                    "1797693134862315708145274237317043567980705675258449965989"
                    "1747680315726078002853876058955863276687817154045895351438"
                    "2464234321326889464182768467546703537516986049910576551282"
                    "0762454900903893289440758685084551339423045832369032229481"
                    "6580855933212334827479782620414472316873817718091929988125"
                    "0404026184124858368.")
        # Warning: In unique mode only the integer digits necessary for
        # uniqueness are computed, the rest are 0. Should we change this?
        assert_equal(fpos32(np.finfo(np.float32).max, precision=0),
                    "340282350000000000000000000000000000000.")

        # test trailing zeros
        assert_equal(fpos32('1.0', unique=False, precision=3), "1.000")
        assert_equal(fpos64('1.0', unique=False, precision=3), "1.000")
        assert_equal(fsci32('1.0', unique=False, precision=3), "1.000e+00")
        assert_equal(fsci64('1.0', unique=False, precision=3), "1.000e+00")
        assert_equal(fpos32('1.5', unique=False, precision=3), "1.500")
        assert_equal(fpos64('1.5', unique=False, precision=3), "1.500")
        assert_equal(fsci32('1.5', unique=False, precision=3), "1.500e+00")
        assert_equal(fsci64('1.5', unique=False, precision=3), "1.500e+00")
        # gh-10713
        assert_equal(fpos64('324', unique=False, precision=5, fractional=False), "324.00")

    def test_dragon4_interface(self):
        tps = [np.float16, np.float32, np.float64]
        if hasattr(np, 'float128'):
            tps.append(np.float128)

        fpos = np.format_float_positional
        fsci = np.format_float_scientific

        for tp in tps:
            # test padding
            assert_equal(fpos(tp('1.0'), pad_left=4, pad_right=4), "   1.    ")
            assert_equal(fpos(tp('-1.0'), pad_left=4, pad_right=4), "  -1.    ")
            assert_equal(fpos(tp('-10.2'),
                         pad_left=4, pad_right=4), " -10.2   ")

            # test exp_digits
            assert_equal(fsci(tp('1.23e1'), exp_digits=5), "1.23e+00001")

            # test fixed (non-unique) mode
            assert_equal(fpos(tp('1.0'), unique=False, precision=4), "1.0000")
            assert_equal(fsci(tp('1.0'), unique=False, precision=4),
                         "1.0000e+00")

            # test trimming
            # trim of 'k' or '.' only affects non-unique mode, since unique
            # mode will not output trailing 0s.
            assert_equal(fpos(tp('1.'), unique=False, precision=4, trim='k'),
                         "1.0000")

            assert_equal(fpos(tp('1.'), unique=False, precision=4, trim='.'),
                         "1.")
            assert_equal(fpos(tp('1.2'), unique=False, precision=4, trim='.'),
                         "1.2" if tp != np.float16 else "1.2002")

            assert_equal(fpos(tp('1.'), unique=False, precision=4, trim='0'),
                         "1.0")
            assert_equal(fpos(tp('1.2'), unique=False, precision=4, trim='0'),
                         "1.2" if tp != np.float16 else "1.2002")
            assert_equal(fpos(tp('1.'), trim='0'), "1.0")

            assert_equal(fpos(tp('1.'), unique=False, precision=4, trim='-'),
                         "1")
            assert_equal(fpos(tp('1.2'), unique=False, precision=4, trim='-'),
                         "1.2" if tp != np.float16 else "1.2002")
            assert_equal(fpos(tp('1.'), trim='-'), "1")

    def float32_roundtrip(self):
        # gh-9360
        x = np.float32(1024 - 2**-14)
        y = np.float32(1024 - 2**-13)
        assert_(repr(x) != repr(y))
        assert_equal(np.float32(repr(x)), x)
        assert_equal(np.float32(repr(y)), y)

    def float64_vs_python(self):
        # gh-2643, gh-6136, gh-6908
        assert_equal(repr(np.float64(0.1)), repr(0.1))
        assert_(repr(np.float64(0.20000000000000004)) != repr(0.2))

if __name__ == "__main__":
    run_module_suite()
