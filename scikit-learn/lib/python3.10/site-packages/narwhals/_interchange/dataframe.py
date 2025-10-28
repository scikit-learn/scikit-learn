from __future__ import annotations

import enum
from typing import TYPE_CHECKING, Any, NoReturn

from narwhals._utils import Version, parse_version

if TYPE_CHECKING:
    import pandas as pd
    import pyarrow as pa
    from typing_extensions import Self, TypeIs

    from narwhals._interchange.series import InterchangeSeries
    from narwhals.dtypes import DType
    from narwhals.stable.v1.typing import DataFrameLike


class DtypeKind(enum.IntEnum):
    # https://data-apis.org/dataframe-protocol/latest/API.html
    INT = 0
    UINT = 1
    FLOAT = 2
    BOOL = 20
    STRING = 21  # UTF-8
    DATETIME = 22
    CATEGORICAL = 23


def map_interchange_dtype_to_narwhals_dtype(  # noqa: C901, PLR0911, PLR0912
    interchange_dtype: tuple[DtypeKind, int, Any, Any],
) -> DType:
    dtypes = Version.V1.dtypes
    if interchange_dtype[0] == DtypeKind.INT:
        if interchange_dtype[1] == 64:
            return dtypes.Int64()
        if interchange_dtype[1] == 32:
            return dtypes.Int32()
        if interchange_dtype[1] == 16:
            return dtypes.Int16()
        if interchange_dtype[1] == 8:
            return dtypes.Int8()
        msg = "Invalid bit width for INT"  # pragma: no cover
        raise AssertionError(msg)
    if interchange_dtype[0] == DtypeKind.UINT:
        if interchange_dtype[1] == 64:
            return dtypes.UInt64()
        if interchange_dtype[1] == 32:
            return dtypes.UInt32()
        if interchange_dtype[1] == 16:
            return dtypes.UInt16()
        if interchange_dtype[1] == 8:
            return dtypes.UInt8()
        msg = "Invalid bit width for UINT"  # pragma: no cover
        raise AssertionError(msg)
    if interchange_dtype[0] == DtypeKind.FLOAT:
        if interchange_dtype[1] == 64:
            return dtypes.Float64()
        if interchange_dtype[1] == 32:
            return dtypes.Float32()
        msg = "Invalid bit width for FLOAT"  # pragma: no cover
        raise AssertionError(msg)
    if interchange_dtype[0] == DtypeKind.BOOL:
        return dtypes.Boolean()
    if interchange_dtype[0] == DtypeKind.STRING:
        return dtypes.String()
    if interchange_dtype[0] == DtypeKind.DATETIME:
        return dtypes.Datetime()
    if interchange_dtype[0] == DtypeKind.CATEGORICAL:  # pragma: no cover
        # upstream issue: https://github.com/ibis-project/ibis/issues/9570
        return dtypes.Categorical()
    msg = f"Invalid dtype, got: {interchange_dtype}"  # pragma: no cover
    raise AssertionError(msg)


class InterchangeFrame:
    _version = Version.V1

    def __init__(self, df: DataFrameLike) -> None:
        self._interchange_frame = df.__dataframe__()

    def __narwhals_dataframe__(self) -> Self:
        return self

    def __native_namespace__(self) -> NoReturn:
        msg = (
            "Cannot access native namespace for interchange-level dataframes with unknown backend."
            "If you would like to see this kind of object supported in Narwhals, please "
            "open a feature request at https://github.com/narwhals-dev/narwhals/issues."
        )
        raise NotImplementedError(msg)

    def get_column(self, name: str) -> InterchangeSeries:
        from narwhals._interchange.series import InterchangeSeries

        return InterchangeSeries(self._interchange_frame.get_column_by_name(name))

    def to_pandas(self) -> pd.DataFrame:
        import pandas as pd  # ignore-banned-import()

        if parse_version(pd) < (1, 5, 0):  # pragma: no cover
            msg = (
                "Conversion to pandas is achieved via interchange protocol which requires"
                f" 'pandas>=1.5.0' to be installed, found {pd.__version__}"
            )
            raise NotImplementedError(msg)
        return pd.api.interchange.from_dataframe(self._interchange_frame)

    def to_arrow(self) -> pa.Table:
        from pyarrow.interchange.from_dataframe import (  # ignore-banned-import()
            from_dataframe,
        )

        return from_dataframe(self._interchange_frame)

    @property
    def schema(self) -> dict[str, DType]:
        return {
            column_name: map_interchange_dtype_to_narwhals_dtype(
                self._interchange_frame.get_column_by_name(column_name).dtype
            )
            for column_name in self._interchange_frame.column_names()
        }

    @property
    def columns(self) -> list[str]:
        return list(self._interchange_frame.column_names())

    def __getattr__(self, attr: str) -> NoReturn:
        msg = (
            f"Attribute {attr} is not supported for interchange-level dataframes.\n\n"
            "Hint: you probably called `nw.from_native` on an object which isn't fully "
            "supported by Narwhals, yet implements `__dataframe__`. If you would like to "
            "see this kind of object supported in Narwhals, please open a feature request "
            "at https://github.com/narwhals-dev/narwhals/issues."
        )
        raise NotImplementedError(msg)

    def simple_select(self, *column_names: str) -> Self:
        frame = self._interchange_frame.select_columns_by_name(list(column_names))
        if not hasattr(frame, "_df"):  # pragma: no cover
            msg = (
                "Expected interchange object to implement `_df` property to allow for recovering original object.\n"
                "See https://github.com/data-apis/dataframe-api/issues/360."
            )
            raise NotImplementedError(msg)
        return self.__class__(frame._df)

    def select(self, *exprs: str) -> Self:  # pragma: no cover
        msg = (
            "`select`-ing not by name is not supported for interchange-only level.\n\n"
            "If you would like to see this kind of object better supported in "
            "Narwhals, please open a feature request "
            "at https://github.com/narwhals-dev/narwhals/issues."
        )
        raise NotImplementedError(msg)


def supports_dataframe_interchange(obj: Any) -> TypeIs[DataFrameLike]:
    return hasattr(obj, "__dataframe__")
