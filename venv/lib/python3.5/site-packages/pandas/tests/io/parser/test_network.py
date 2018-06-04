# -*- coding: utf-8 -*-

"""
Tests parsers ability to read and parse non-local files
and hence require a network connection to be read.
"""
import logging

import pytest
import numpy as np

import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas import DataFrame
from pandas.io.parsers import read_csv, read_table
from pandas.compat import BytesIO, StringIO


@pytest.mark.network
@pytest.mark.parametrize(
    "compress_type, extension", [
        ('gzip', '.gz'), ('bz2', '.bz2'), ('zip', '.zip'),
        pytest.param('xz', '.xz', marks=td.skip_if_no_lzma)
    ]
)
@pytest.mark.parametrize('mode', ['explicit', 'infer'])
@pytest.mark.parametrize('engine', ['python', 'c'])
def test_compressed_urls(salaries_table, compress_type, extension, mode,
                         engine):
    check_compressed_urls(salaries_table, compress_type, extension, mode,
                          engine)


@tm.network
def check_compressed_urls(salaries_table, compression, extension, mode,
                          engine):
    # test reading compressed urls with various engines and
    # extension inference
    base_url = ('https://github.com/pandas-dev/pandas/raw/master/'
                'pandas/tests/io/parser/data/salaries.csv')

    url = base_url + extension

    if mode != 'explicit':
        compression = mode

    url_table = read_table(url, compression=compression, engine=engine)
    tm.assert_frame_equal(url_table, salaries_table)


@pytest.mark.usefixtures("s3_resource")
class TestS3(object):

    def test_parse_public_s3_bucket(self):
        pytest.importorskip('s3fs')
        # more of an integration test due to the not-public contents portion
        # can probably mock this though.
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df = read_csv('s3://pandas-test/tips.csv' +
                          ext, compression=comp)
            assert isinstance(df, DataFrame)
            assert not df.empty
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')), df)

        # Read public file from bucket with not-public contents
        df = read_csv('s3://cant_get_it/tips.csv')
        assert isinstance(df, DataFrame)
        assert not df.empty
        tm.assert_frame_equal(read_csv(tm.get_data_path('tips.csv')), df)

    def test_parse_public_s3n_bucket(self):

        # Read from AWS s3 as "s3n" URL
        df = read_csv('s3n://pandas-test/tips.csv', nrows=10)
        assert isinstance(df, DataFrame)
        assert not df.empty
        tm.assert_frame_equal(read_csv(
            tm.get_data_path('tips.csv')).iloc[:10], df)

    def test_parse_public_s3a_bucket(self):
        # Read from AWS s3 as "s3a" URL
        df = read_csv('s3a://pandas-test/tips.csv', nrows=10)
        assert isinstance(df, DataFrame)
        assert not df.empty
        tm.assert_frame_equal(read_csv(
            tm.get_data_path('tips.csv')).iloc[:10], df)

    def test_parse_public_s3_bucket_nrows(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df = read_csv('s3://pandas-test/tips.csv' +
                          ext, nrows=10, compression=comp)
            assert isinstance(df, DataFrame)
            assert not df.empty
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')).iloc[:10], df)

    def test_parse_public_s3_bucket_chunked(self):
        # Read with a chunksize
        chunksize = 5
        local_tips = read_csv(tm.get_data_path('tips.csv'))
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df_reader = read_csv('s3://pandas-test/tips.csv' + ext,
                                 chunksize=chunksize, compression=comp)
            assert df_reader.chunksize == chunksize
            for i_chunk in [0, 1, 2]:
                # Read a couple of chunks and make sure we see them
                # properly.
                df = df_reader.get_chunk()
                assert isinstance(df, DataFrame)
                assert not df.empty
                true_df = local_tips.iloc[
                    chunksize * i_chunk: chunksize * (i_chunk + 1)]
                tm.assert_frame_equal(true_df, df)

    def test_parse_public_s3_bucket_chunked_python(self):
        # Read with a chunksize using the Python parser
        chunksize = 5
        local_tips = read_csv(tm.get_data_path('tips.csv'))
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df_reader = read_csv('s3://pandas-test/tips.csv' + ext,
                                 chunksize=chunksize, compression=comp,
                                 engine='python')
            assert df_reader.chunksize == chunksize
            for i_chunk in [0, 1, 2]:
                # Read a couple of chunks and make sure we see them properly.
                df = df_reader.get_chunk()
                assert isinstance(df, DataFrame)
                assert not df.empty
                true_df = local_tips.iloc[
                    chunksize * i_chunk: chunksize * (i_chunk + 1)]
                tm.assert_frame_equal(true_df, df)

    def test_parse_public_s3_bucket_python(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df = read_csv('s3://pandas-test/tips.csv' + ext, engine='python',
                          compression=comp)
            assert isinstance(df, DataFrame)
            assert not df.empty
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')), df)

    def test_infer_s3_compression(self):
        for ext in ['', '.gz', '.bz2']:
            df = read_csv('s3://pandas-test/tips.csv' + ext,
                          engine='python', compression='infer')
            assert isinstance(df, DataFrame)
            assert not df.empty
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')), df)

    def test_parse_public_s3_bucket_nrows_python(self):
        for ext, comp in [('', None), ('.gz', 'gzip'), ('.bz2', 'bz2')]:
            df = read_csv('s3://pandas-test/tips.csv' + ext, engine='python',
                          nrows=10, compression=comp)
            assert isinstance(df, DataFrame)
            assert not df.empty
            tm.assert_frame_equal(read_csv(
                tm.get_data_path('tips.csv')).iloc[:10], df)

    def test_s3_fails(self):
        with pytest.raises(IOError):
            read_csv('s3://nyqpug/asdf.csv')

        # Receive a permission error when trying to read a private bucket.
        # It's irrelevant here that this isn't actually a table.
        with pytest.raises(IOError):
            read_csv('s3://cant_get_it/')

    def test_read_csv_handles_boto_s3_object(self,
                                             s3_resource,
                                             tips_file):
        # see gh-16135

        s3_object = s3_resource.meta.client.get_object(
            Bucket='pandas-test',
            Key='tips.csv')

        result = read_csv(BytesIO(s3_object["Body"].read()), encoding='utf8')
        assert isinstance(result, DataFrame)
        assert not result.empty

        expected = read_csv(tips_file)
        tm.assert_frame_equal(result, expected)

    def test_read_csv_chunked_download(self, s3_resource, caplog):
        # 8 MB, S3FS usees 5MB chunks
        df = DataFrame(np.random.randn(100000, 4), columns=list('abcd'))
        buf = BytesIO()
        str_buf = StringIO()

        df.to_csv(str_buf)

        buf = BytesIO(str_buf.getvalue().encode('utf-8'))

        s3_resource.Bucket("pandas-test").put_object(
            Key="large-file.csv",
            Body=buf)

        with caplog.at_level(logging.DEBUG, logger='s3fs.core'):
            read_csv("s3://pandas-test/large-file.csv", nrows=5)
            # log of fetch_range (start, stop)
            assert ((0, 5505024) in set(x.args[-2:] for x in caplog.records))
