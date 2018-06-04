import pandas.util.testing as tm


class BaseExtensionTests(object):
    assert_series_equal = staticmethod(tm.assert_series_equal)
    assert_frame_equal = staticmethod(tm.assert_frame_equal)
    assert_extension_array_equal = staticmethod(
        tm.assert_extension_array_equal
    )
