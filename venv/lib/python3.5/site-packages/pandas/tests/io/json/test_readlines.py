# -*- coding: utf-8 -*-
import pytest
import pandas as pd
from pandas import DataFrame, read_json
from pandas.compat import StringIO
from pandas.io.json.json import JsonReader
import pandas.util.testing as tm
from pandas.util.testing import (assert_frame_equal, assert_series_equal,
                                 ensure_clean)


@pytest.fixture
def lines_json_df():
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    return df.to_json(lines=True, orient="records")


def test_read_jsonl():
    # GH9180
    result = read_json('{"a": 1, "b": 2}\n{"b":2, "a" :1}\n', lines=True)
    expected = DataFrame([[1, 2], [1, 2]], columns=['a', 'b'])
    assert_frame_equal(result, expected)


def test_read_jsonl_unicode_chars():
    # GH15132: non-ascii unicode characters
    # \u201d == RIGHT DOUBLE QUOTATION MARK

    # simulate file handle
    json = '{"a": "foo”", "b": "bar"}\n{"a": "foo", "b": "bar"}\n'
    json = StringIO(json)
    result = read_json(json, lines=True)
    expected = DataFrame([[u"foo\u201d", "bar"], ["foo", "bar"]],
                         columns=['a', 'b'])
    assert_frame_equal(result, expected)

    # simulate string
    json = '{"a": "foo”", "b": "bar"}\n{"a": "foo", "b": "bar"}\n'
    result = read_json(json, lines=True)
    expected = DataFrame([[u"foo\u201d", "bar"], ["foo", "bar"]],
                         columns=['a', 'b'])
    assert_frame_equal(result, expected)


def test_to_jsonl():
    # GH9180
    df = DataFrame([[1, 2], [1, 2]], columns=['a', 'b'])
    result = df.to_json(orient="records", lines=True)
    expected = '{"a":1,"b":2}\n{"a":1,"b":2}'
    assert result == expected

    df = DataFrame([["foo}", "bar"], ['foo"', "bar"]], columns=['a', 'b'])
    result = df.to_json(orient="records", lines=True)
    expected = '{"a":"foo}","b":"bar"}\n{"a":"foo\\"","b":"bar"}'
    assert result == expected
    assert_frame_equal(read_json(result, lines=True), df)

    # GH15096: escaped characters in columns and data
    df = DataFrame([["foo\\", "bar"], ['foo"', "bar"]],
                   columns=["a\\", 'b'])
    result = df.to_json(orient="records", lines=True)
    expected = ('{"a\\\\":"foo\\\\","b":"bar"}\n'
                '{"a\\\\":"foo\\"","b":"bar"}')
    assert result == expected
    assert_frame_equal(read_json(result, lines=True), df)


@pytest.mark.parametrize("chunksize", [1, 1.0])
def test_readjson_chunks(lines_json_df, chunksize):
    # Basic test that read_json(chunks=True) gives the same result as
    # read_json(chunks=False)
    # GH17048: memory usage when lines=True

    unchunked = read_json(StringIO(lines_json_df), lines=True)
    reader = read_json(StringIO(lines_json_df), lines=True,
                       chunksize=chunksize)
    chunked = pd.concat(reader)

    assert_frame_equal(chunked, unchunked)


def test_readjson_chunksize_requires_lines(lines_json_df):
    msg = "chunksize can only be passed if lines=True"
    with tm.assert_raises_regex(ValueError, msg):
        pd.read_json(StringIO(lines_json_df), lines=False, chunksize=2)


def test_readjson_chunks_series():
    # Test reading line-format JSON to Series with chunksize param
    s = pd.Series({'A': 1, 'B': 2})

    strio = StringIO(s.to_json(lines=True, orient="records"))
    unchunked = pd.read_json(strio, lines=True, typ='Series')

    strio = StringIO(s.to_json(lines=True, orient="records"))
    chunked = pd.concat(pd.read_json(
        strio, lines=True, typ='Series', chunksize=1
    ))

    assert_series_equal(chunked, unchunked)


def test_readjson_each_chunk(lines_json_df):
    # Other tests check that the final result of read_json(chunksize=True)
    # is correct. This checks the intermediate chunks.
    chunks = list(
        pd.read_json(StringIO(lines_json_df), lines=True, chunksize=2)
    )
    assert chunks[0].shape == (2, 2)
    assert chunks[1].shape == (1, 2)


def test_readjson_chunks_from_file():
    with ensure_clean('test.json') as path:
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        df.to_json(path, lines=True, orient="records")
        chunked = pd.concat(pd.read_json(path, lines=True, chunksize=1))
        unchunked = pd.read_json(path, lines=True)
        assert_frame_equal(unchunked, chunked)


@pytest.mark.parametrize("chunksize", [None, 1])
def test_readjson_chunks_closes(chunksize):
    with ensure_clean('test.json') as path:
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        df.to_json(path, lines=True, orient="records")
        reader = JsonReader(
            path, orient=None, typ="frame", dtype=True, convert_axes=True,
            convert_dates=True, keep_default_dates=True, numpy=False,
            precise_float=False, date_unit=None, encoding=None,
            lines=True, chunksize=chunksize, compression=None)
        reader.read()
        assert reader.open_stream.closed, "didn't close stream with \
            chunksize = {chunksize}".format(chunksize=chunksize)


@pytest.mark.parametrize("chunksize", [0, -1, 2.2, "foo"])
def test_readjson_invalid_chunksize(lines_json_df, chunksize):
    msg = r"'chunksize' must be an integer >=1"

    with tm.assert_raises_regex(ValueError, msg):
        pd.read_json(StringIO(lines_json_df), lines=True,
                     chunksize=chunksize)


@pytest.mark.parametrize("chunksize", [None, 1, 2])
def test_readjson_chunks_multiple_empty_lines(chunksize):
    j = """

    {"A":1,"B":4}



    {"A":2,"B":5}







    {"A":3,"B":6}
    """
    orig = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    test = pd.read_json(j, lines=True, chunksize=chunksize)
    if chunksize is not None:
        test = pd.concat(test)
    tm.assert_frame_equal(
        orig, test, obj="chunksize: {chunksize}".format(chunksize=chunksize))
