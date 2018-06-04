from pandas.io.common import is_s3_url


class TestS3URL(object):

    def test_is_s3_url(self):
        assert is_s3_url("s3://pandas/somethingelse.com")
        assert not is_s3_url("s4://pandas/somethingelse.com")
