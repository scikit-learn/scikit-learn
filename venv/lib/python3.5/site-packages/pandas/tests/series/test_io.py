# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime
import collections
import pytest

import numpy as np
import pandas as pd

from pandas import Series, DataFrame

from pandas.compat import StringIO, u
from pandas.util.testing import (assert_series_equal, assert_almost_equal,
                                 assert_frame_equal, ensure_clean)
import pandas.util.testing as tm

from .common import TestData


class TestSeriesToCSV(TestData):

    def read_csv(self, path, **kwargs):
        params = dict(squeeze=True, index_col=0,
                      header=None, parse_dates=True)
        params.update(**kwargs)

        header = params.get("header")
        out = pd.read_csv(path, **params)

        if header is None:
            out.name = out.index.name = None

        return out

    def test_from_csv_deprecation(self):
        # see gh-17812
        with ensure_clean() as path:
            self.ts.to_csv(path)

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                ts = self.read_csv(path)
                depr_ts = Series.from_csv(path)
                assert_series_equal(depr_ts, ts)

    def test_from_csv(self):

        with ensure_clean() as path:
            self.ts.to_csv(path)
            ts = self.read_csv(path)
            assert_series_equal(self.ts, ts, check_names=False)

            assert ts.name is None
            assert ts.index.name is None

            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                depr_ts = Series.from_csv(path)
                assert_series_equal(depr_ts, ts)

            # see gh-10483
            self.ts.to_csv(path, header=True)
            ts_h = self.read_csv(path, header=0)
            assert ts_h.name == "ts"

            self.series.to_csv(path)
            series = self.read_csv(path)
            assert_series_equal(self.series, series, check_names=False)

            assert series.name is None
            assert series.index.name is None

            self.series.to_csv(path, header=True)
            series_h = self.read_csv(path, header=0)
            assert series_h.name == "series"

            outfile = open(path, "w")
            outfile.write("1998-01-01|1.0\n1999-01-01|2.0")
            outfile.close()

            series = self.read_csv(path, sep="|")
            check_series = Series({datetime(1998, 1, 1): 1.0,
                                   datetime(1999, 1, 1): 2.0})
            assert_series_equal(check_series, series)

            series = self.read_csv(path, sep="|", parse_dates=False)
            check_series = Series({"1998-01-01": 1.0, "1999-01-01": 2.0})
            assert_series_equal(check_series, series)

    def test_to_csv(self):
        import io

        with ensure_clean() as path:
            self.ts.to_csv(path)

            with io.open(path, newline=None) as f:
                lines = f.readlines()
            assert (lines[1] != '\n')

            self.ts.to_csv(path, index=False)
            arr = np.loadtxt(path)
            assert_almost_equal(arr, self.ts.values)

    def test_to_csv_unicode_index(self):
        buf = StringIO()
        s = Series([u("\u05d0"), "d2"], index=[u("\u05d0"), u("\u05d1")])

        s.to_csv(buf, encoding="UTF-8")
        buf.seek(0)

        s2 = self.read_csv(buf, index_col=0, encoding="UTF-8")
        assert_series_equal(s, s2)

    def test_to_csv_float_format(self):

        with ensure_clean() as filename:
            ser = Series([0.123456, 0.234567, 0.567567])
            ser.to_csv(filename, float_format="%.2f")

            rs = self.read_csv(filename)
            xp = Series([0.12, 0.23, 0.57])
            assert_series_equal(rs, xp)

    def test_to_csv_list_entries(self):
        s = Series(['jack and jill', 'jesse and frank'])

        split = s.str.split(r'\s+and\s+')

        buf = StringIO()
        split.to_csv(buf)

    def test_to_csv_path_is_none(self):
        # GH 8215
        # Series.to_csv() was returning None, inconsistent with
        # DataFrame.to_csv() which returned string
        s = Series([1, 2, 3])
        csv_str = s.to_csv(path=None)
        assert isinstance(csv_str, str)

    def test_to_csv_compression(self, compression):

        s = Series([0.123456, 0.234567, 0.567567], index=['A', 'B', 'C'],
                   name='X')

        with ensure_clean() as filename:

            s.to_csv(filename, compression=compression, header=True)

            # test the round trip - to_csv -> read_csv
            rs = pd.read_csv(filename, compression=compression,
                             index_col=0, squeeze=True)
            assert_series_equal(s, rs)

            # explicitly ensure file was compressed
            with tm.decompress_file(filename, compression) as fh:
                text = fh.read().decode('utf8')
                assert s.name in text

            with tm.decompress_file(filename, compression) as fh:
                assert_series_equal(s, pd.read_csv(fh,
                                                   index_col=0,
                                                   squeeze=True))


class TestSeriesIO(TestData):

    def test_to_frame(self):
        self.ts.name = None
        rs = self.ts.to_frame()
        xp = pd.DataFrame(self.ts.values, index=self.ts.index)
        assert_frame_equal(rs, xp)

        self.ts.name = 'testname'
        rs = self.ts.to_frame()
        xp = pd.DataFrame(dict(testname=self.ts.values), index=self.ts.index)
        assert_frame_equal(rs, xp)

        rs = self.ts.to_frame(name='testdifferent')
        xp = pd.DataFrame(
            dict(testdifferent=self.ts.values), index=self.ts.index)
        assert_frame_equal(rs, xp)

    def test_timeseries_periodindex(self):
        # GH2891
        from pandas import period_range
        prng = period_range('1/1/2011', '1/1/2012', freq='M')
        ts = Series(np.random.randn(len(prng)), prng)
        new_ts = tm.round_trip_pickle(ts)
        assert new_ts.index.freq == 'M'

    def test_pickle_preserve_name(self):
        for n in [777, 777., 'name', datetime(2001, 11, 11), (1, 2)]:
            unpickled = self._pickle_roundtrip_name(tm.makeTimeSeries(name=n))
            assert unpickled.name == n

    def _pickle_roundtrip_name(self, obj):

        with ensure_clean() as path:
            obj.to_pickle(path)
            unpickled = pd.read_pickle(path)
            return unpickled

    def test_to_frame_expanddim(self):
        # GH 9762

        class SubclassedSeries(Series):

            @property
            def _constructor_expanddim(self):
                return SubclassedFrame

        class SubclassedFrame(DataFrame):
            pass

        s = SubclassedSeries([1, 2, 3], name='X')
        result = s.to_frame()
        assert isinstance(result, SubclassedFrame)
        expected = SubclassedFrame({'X': [1, 2, 3]})
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('mapping', (
        dict,
        collections.defaultdict(list),
        collections.OrderedDict))
    def test_to_dict(self, mapping):
        # GH16122
        ts = TestData().ts
        tm.assert_series_equal(
            Series(ts.to_dict(mapping), name='ts'), ts)
        from_method = Series(ts.to_dict(collections.Counter))
        from_constructor = Series(collections.Counter(ts.iteritems()))
        tm.assert_series_equal(from_method, from_constructor)
