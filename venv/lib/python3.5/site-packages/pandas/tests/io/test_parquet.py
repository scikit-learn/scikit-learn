""" test parquet compat """

import pytest
import datetime
from distutils.version import LooseVersion
from warnings import catch_warnings

import numpy as np
import pandas as pd
from pandas.compat import PY3, is_platform_windows, is_platform_mac
from pandas.io.parquet import (to_parquet, read_parquet, get_engine,
                               PyArrowImpl, FastParquetImpl)
from pandas.util import testing as tm

try:
    import pyarrow  # noqa
    _HAVE_PYARROW = True
except ImportError:
    _HAVE_PYARROW = False

try:
    import fastparquet  # noqa
    _HAVE_FASTPARQUET = True
except ImportError:
    _HAVE_FASTPARQUET = False


# setup engines & skips
@pytest.fixture(params=[
    pytest.param('fastparquet',
                 marks=pytest.mark.skipif(not _HAVE_FASTPARQUET,
                                          reason='fastparquet is '
                                                 'not installed')),
    pytest.param('pyarrow',
                 marks=pytest.mark.skipif(not _HAVE_PYARROW,
                                          reason='pyarrow is '
                                                 'not installed'))])
def engine(request):
    return request.param


@pytest.fixture
def pa():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    return 'pyarrow'


@pytest.fixture
def pa_lt_070():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    if LooseVersion(pyarrow.__version__) >= LooseVersion('0.7.0'):
        pytest.skip("pyarrow is >= 0.7.0")
    return 'pyarrow'


@pytest.fixture
def pa_ge_070():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    if LooseVersion(pyarrow.__version__) < LooseVersion('0.7.0'):
        pytest.skip("pyarrow is < 0.7.0")
    return 'pyarrow'


@pytest.fixture
def fp():
    if not _HAVE_FASTPARQUET:
        pytest.skip("fastparquet is not installed")
    return 'fastparquet'


@pytest.fixture
def fp_lt_014():
    if not _HAVE_FASTPARQUET:
        pytest.skip("fastparquet is not installed")
    if LooseVersion(fastparquet.__version__) >= LooseVersion('0.1.4'):
        pytest.skip("fastparquet is >= 0.1.4")
    return 'fastparquet'


@pytest.fixture
def df_compat():
    return pd.DataFrame({'A': [1, 2, 3], 'B': 'foo'})


@pytest.fixture
def df_cross_compat():
    df = pd.DataFrame({'a': list('abc'),
                       'b': list(range(1, 4)),
                       # 'c': np.arange(3, 6).astype('u1'),
                       'd': np.arange(4.0, 7.0, dtype='float64'),
                       'e': [True, False, True],
                       'f': pd.date_range('20130101', periods=3),
                       # 'g': pd.date_range('20130101', periods=3,
                       #                    tz='US/Eastern'),
                       # 'h': pd.date_range('20130101', periods=3, freq='ns')
                       })
    return df


@pytest.fixture
def df_full():
    return pd.DataFrame(
        {'string': list('abc'),
         'string_with_nan': ['a', np.nan, 'c'],
         'string_with_none': ['a', None, 'c'],
         'bytes': [b'foo', b'bar', b'baz'],
         'unicode': [u'foo', u'bar', u'baz'],
         'int': list(range(1, 4)),
         'uint': np.arange(3, 6).astype('u1'),
         'float': np.arange(4.0, 7.0, dtype='float64'),
         'float_with_nan': [2., np.nan, 3.],
         'bool': [True, False, True],
         'datetime': pd.date_range('20130101', periods=3),
         'datetime_with_nat': [pd.Timestamp('20130101'),
                               pd.NaT,
                               pd.Timestamp('20130103')]})


def check_round_trip(df, engine=None, path=None,
                     write_kwargs=None, read_kwargs=None,
                     expected=None, check_names=True,
                     repeat=2):
    """Verify parquet serializer and deserializer produce the same results.

    Performs a pandas to disk and disk to pandas round trip,
    then compares the 2 resulting DataFrames to verify equality.

    Parameters
    ----------
    df: Dataframe
    engine: str, optional
        'pyarrow' or 'fastparquet'
    path: str, optional
    write_kwargs: dict of str:str, optional
    read_kwargs: dict of str:str, optional
    expected: DataFrame, optional
        Expected deserialization result, otherwise will be equal to `df`
    check_names: list of str, optional
        Closed set of column names to be compared
    repeat: int, optional
        How many times to repeat the test
    """

    write_kwargs = write_kwargs or {'compression': None}
    read_kwargs = read_kwargs or {}

    if expected is None:
        expected = df

    if engine:
        write_kwargs['engine'] = engine
        read_kwargs['engine'] = engine

    def compare(repeat):
        for _ in range(repeat):
            df.to_parquet(path, **write_kwargs)
            with catch_warnings(record=True):
                actual = read_parquet(path, **read_kwargs)
            tm.assert_frame_equal(expected, actual,
                                  check_names=check_names)

    if path is None:
        with tm.ensure_clean() as path:
            compare(repeat)
    else:
        compare(repeat)


def test_invalid_engine(df_compat):
    with pytest.raises(ValueError):
        check_round_trip(df_compat, 'foo', 'bar')


def test_options_py(df_compat, pa):
    # use the set option

    with pd.option_context('io.parquet.engine', 'pyarrow'):
        check_round_trip(df_compat)


def test_options_fp(df_compat, fp):
    # use the set option

    with pd.option_context('io.parquet.engine', 'fastparquet'):
        check_round_trip(df_compat)


def test_options_auto(df_compat, fp, pa):
    # use the set option

    with pd.option_context('io.parquet.engine', 'auto'):
        check_round_trip(df_compat)


def test_options_get_engine(fp, pa):
    assert isinstance(get_engine('pyarrow'), PyArrowImpl)
    assert isinstance(get_engine('fastparquet'), FastParquetImpl)

    with pd.option_context('io.parquet.engine', 'pyarrow'):
        assert isinstance(get_engine('auto'), PyArrowImpl)
        assert isinstance(get_engine('pyarrow'), PyArrowImpl)
        assert isinstance(get_engine('fastparquet'), FastParquetImpl)

    with pd.option_context('io.parquet.engine', 'fastparquet'):
        assert isinstance(get_engine('auto'), FastParquetImpl)
        assert isinstance(get_engine('pyarrow'), PyArrowImpl)
        assert isinstance(get_engine('fastparquet'), FastParquetImpl)

    with pd.option_context('io.parquet.engine', 'auto'):
        assert isinstance(get_engine('auto'), PyArrowImpl)
        assert isinstance(get_engine('pyarrow'), PyArrowImpl)
        assert isinstance(get_engine('fastparquet'), FastParquetImpl)


@pytest.mark.xfail(is_platform_windows() or is_platform_mac(),
                   reason="reading pa metadata failing on Windows/mac")
def test_cross_engine_pa_fp(df_cross_compat, pa, fp):
    # cross-compat with differing reading/writing engines

    df = df_cross_compat
    with tm.ensure_clean() as path:
        df.to_parquet(path, engine=pa, compression=None)

        result = read_parquet(path, engine=fp)
        tm.assert_frame_equal(result, df)

        result = read_parquet(path, engine=fp, columns=['a', 'd'])
        tm.assert_frame_equal(result, df[['a', 'd']])


def test_cross_engine_fp_pa(df_cross_compat, pa, fp):
    # cross-compat with differing reading/writing engines

    df = df_cross_compat
    with tm.ensure_clean() as path:
        df.to_parquet(path, engine=fp, compression=None)

        with catch_warnings(record=True):
            result = read_parquet(path, engine=pa)
            tm.assert_frame_equal(result, df)

            result = read_parquet(path, engine=pa, columns=['a', 'd'])
            tm.assert_frame_equal(result, df[['a', 'd']])


class Base(object):

    def check_error_on_write(self, df, engine, exc):
        # check that we are raising the exception on writing
        with tm.ensure_clean() as path:
            with pytest.raises(exc):
                to_parquet(df, path, engine, compression=None)


class TestBasic(Base):

    def test_error(self, engine):
        for obj in [pd.Series([1, 2, 3]), 1, 'foo', pd.Timestamp('20130101'),
                    np.array([1, 2, 3])]:
            self.check_error_on_write(obj, engine, ValueError)

    def test_columns_dtypes(self, engine):
        df = pd.DataFrame({'string': list('abc'),
                           'int': list(range(1, 4))})

        # unicode
        df.columns = [u'foo', u'bar']
        check_round_trip(df, engine)

    def test_columns_dtypes_invalid(self, engine):
        df = pd.DataFrame({'string': list('abc'),
                           'int': list(range(1, 4))})

        # numeric
        df.columns = [0, 1]
        self.check_error_on_write(df, engine, ValueError)

        if PY3:
            # bytes on PY3, on PY2 these are str
            df.columns = [b'foo', b'bar']
            self.check_error_on_write(df, engine, ValueError)

        # python object
        df.columns = [datetime.datetime(2011, 1, 1, 0, 0),
                      datetime.datetime(2011, 1, 1, 1, 1)]
        self.check_error_on_write(df, engine, ValueError)

    @pytest.mark.parametrize('compression', [None, 'gzip', 'snappy', 'brotli'])
    def test_compression(self, engine, compression):

        if compression == 'snappy':
            pytest.importorskip('snappy')

        elif compression == 'brotli':
            pytest.importorskip('brotli')

        df = pd.DataFrame({'A': [1, 2, 3]})
        check_round_trip(df, engine, write_kwargs={'compression': compression})

    def test_read_columns(self, engine):
        # GH18154
        df = pd.DataFrame({'string': list('abc'),
                           'int': list(range(1, 4))})

        expected = pd.DataFrame({'string': list('abc')})
        check_round_trip(df, engine, expected=expected,
                         read_kwargs={'columns': ['string']})

    def test_write_index(self, engine):
        check_names = engine != 'fastparquet'

        if engine == 'pyarrow':
            import pyarrow
            if LooseVersion(pyarrow.__version__) < LooseVersion('0.7.0'):
                pytest.skip("pyarrow is < 0.7.0")

        df = pd.DataFrame({'A': [1, 2, 3]})
        check_round_trip(df, engine)

        indexes = [
            [2, 3, 4],
            pd.date_range('20130101', periods=3),
            list('abc'),
            [1, 3, 4],
        ]
        # non-default index
        for index in indexes:
            df.index = index
            check_round_trip(df, engine, check_names=check_names)

        # index with meta-data
        df.index = [0, 1, 2]
        df.index.name = 'foo'
        check_round_trip(df, engine)

    def test_write_multiindex(self, pa_ge_070):
        # Not suppoprted in fastparquet as of 0.1.3 or older pyarrow version
        engine = pa_ge_070

        df = pd.DataFrame({'A': [1, 2, 3]})
        index = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)])
        df.index = index
        check_round_trip(df, engine)

    def test_write_column_multiindex(self, engine):
        # column multi-index
        mi_columns = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)])
        df = pd.DataFrame(np.random.randn(4, 3), columns=mi_columns)
        self.check_error_on_write(df, engine, ValueError)

    def test_multiindex_with_columns(self, pa_ge_070):
        engine = pa_ge_070
        dates = pd.date_range('01-Jan-2018', '01-Dec-2018', freq='MS')
        df = pd.DataFrame(np.random.randn(2 * len(dates), 3),
                          columns=list('ABC'))
        index1 = pd.MultiIndex.from_product(
            [['Level1', 'Level2'], dates],
            names=['level', 'date'])
        index2 = index1.copy(names=None)
        for index in [index1, index2]:
            df.index = index

            check_round_trip(df, engine)
            check_round_trip(df, engine, read_kwargs={'columns': ['A', 'B']},
                             expected=df[['A', 'B']])


class TestParquetPyArrow(Base):

    def test_basic(self, pa, df_full):

        df = df_full

        # additional supported types for pyarrow
        import pyarrow
        if LooseVersion(pyarrow.__version__) >= LooseVersion('0.7.0'):
            df['datetime_tz'] = pd.date_range('20130101', periods=3,
                                              tz='Europe/Brussels')
        df['bool_with_none'] = [True, None, True]

        check_round_trip(df, pa)

    @pytest.mark.xfail(reason="pyarrow fails on this (ARROW-1883)")
    def test_basic_subset_columns(self, pa, df_full):
        # GH18628

        df = df_full
        # additional supported types for pyarrow
        df['datetime_tz'] = pd.date_range('20130101', periods=3,
                                          tz='Europe/Brussels')

        check_round_trip(df, pa, expected=df[['string', 'int']],
                         read_kwargs={'columns': ['string', 'int']})

    def test_duplicate_columns(self, pa):
        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3),
                          columns=list('aaa')).copy()
        self.check_error_on_write(df, pa, ValueError)

    def test_unsupported(self, pa):
        # period
        df = pd.DataFrame({'a': pd.period_range('2013', freq='M', periods=3)})
        self.check_error_on_write(df, pa, ValueError)

        # timedelta
        df = pd.DataFrame({'a': pd.timedelta_range('1 day',
                                                   periods=3)})
        self.check_error_on_write(df, pa, NotImplementedError)

        # mixed python objects
        df = pd.DataFrame({'a': ['a', 1, 2.0]})
        self.check_error_on_write(df, pa, ValueError)

    def test_categorical(self, pa_ge_070):
        pa = pa_ge_070

        # supported in >= 0.7.0
        df = pd.DataFrame({'a': pd.Categorical(list('abc'))})

        # de-serialized as object
        expected = df.assign(a=df.a.astype(object))
        check_round_trip(df, pa, expected=expected)

    def test_categorical_unsupported(self, pa_lt_070):
        pa = pa_lt_070

        # supported in >= 0.7.0
        df = pd.DataFrame({'a': pd.Categorical(list('abc'))})
        self.check_error_on_write(df, pa, NotImplementedError)

    def test_s3_roundtrip(self, df_compat, s3_resource, pa):
        # GH #19134
        check_round_trip(df_compat, pa,
                         path='s3://pandas-test/pyarrow.parquet')


class TestParquetFastParquet(Base):

    def test_basic(self, fp, df_full):
        df = df_full

        # additional supported types for fastparquet
        if LooseVersion(fastparquet.__version__) >= LooseVersion('0.1.4'):
            df['datetime_tz'] = pd.date_range('20130101', periods=3,
                                              tz='US/Eastern')
        df['timedelta'] = pd.timedelta_range('1 day', periods=3)
        check_round_trip(df, fp)

    @pytest.mark.skip(reason="not supported")
    def test_duplicate_columns(self, fp):

        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3),
                          columns=list('aaa')).copy()
        self.check_error_on_write(df, fp, ValueError)

    def test_bool_with_none(self, fp):
        df = pd.DataFrame({'a': [True, None, False]})
        expected = pd.DataFrame({'a': [1.0, np.nan, 0.0]}, dtype='float16')
        check_round_trip(df, fp, expected=expected)

    def test_unsupported(self, fp):

        # period
        df = pd.DataFrame({'a': pd.period_range('2013', freq='M', periods=3)})
        self.check_error_on_write(df, fp, ValueError)

        # mixed
        df = pd.DataFrame({'a': ['a', 1, 2.0]})
        self.check_error_on_write(df, fp, ValueError)

    def test_categorical(self, fp):
        if LooseVersion(fastparquet.__version__) < LooseVersion("0.1.3"):
            pytest.skip("CategoricalDtype not supported for older fp")
        df = pd.DataFrame({'a': pd.Categorical(list('abc'))})
        check_round_trip(df, fp)

    def test_datetime_tz(self, fp_lt_014):

        # fastparquet<0.1.4 doesn't preserve tz
        df = pd.DataFrame({'a': pd.date_range('20130101', periods=3,
                                              tz='US/Eastern')})
        # warns on the coercion
        with catch_warnings(record=True):
            check_round_trip(df, fp_lt_014,
                             expected=df.astype('datetime64[ns]'))

    def test_filter_row_groups(self, fp):
        d = {'a': list(range(0, 3))}
        df = pd.DataFrame(d)
        with tm.ensure_clean() as path:
            df.to_parquet(path, fp, compression=None,
                          row_group_offsets=1)
            result = read_parquet(path, fp, filters=[('a', '==', 0)])
        assert len(result) == 1

    def test_s3_roundtrip(self, df_compat, s3_resource, fp):
        # GH #19134
        check_round_trip(df_compat, fp,
                         path='s3://pandas-test/fastparquet.parquet')
