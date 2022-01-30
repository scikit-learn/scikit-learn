"""
Helpers for sharing tests between DataFrame/Series
"""

from pandas import DataFrame


def get_dtype(obj):
    if isinstance(obj, DataFrame):
        # Note: we are assuming only one column
        return obj.dtypes.iat[0]
    else:
        return obj.dtype


def get_obj(df: DataFrame, klass):
    """
    For sharing tests using frame_or_series, either return the DataFrame
    unchanged or return it's first column as a Series.
    """
    if klass is DataFrame:
        return df
    return df._ixs(0, axis=1)
