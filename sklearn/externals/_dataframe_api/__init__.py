from ...utils.validation import _is_pandas_df
from ...utils.validation import _is_polars_df

__all__ = ["get_dataframe_standard", "has_supported_dataframe_standards"]


def get_dataframe_standard(df):
    if hasattr(df, "__dataframe_standard__"):
        return df.__dataframe_standard__()
    elif _is_pandas_df(df):
        from .pandas_standard import dataframe_standard as pandas_dataframe_standard

        return pandas_dataframe_standard(df)
    elif _is_polars_df(df):
        from .polars_standard import dataframe_standard as polars_dataframe_standard

        return polars_dataframe_standard(df)
    else:
        raise ValueError("Only pandas and polars DataFrames are supported.")


def has_supported_dataframe_standards(df):
    return (
        hasattr(df, "__dataframe_standard__") or _is_pandas_df(df) or _is_polars_df(df)
    )
