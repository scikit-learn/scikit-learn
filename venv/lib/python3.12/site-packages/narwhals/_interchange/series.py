from __future__ import annotations

from typing import TYPE_CHECKING, Any, NoReturn

from narwhals._interchange.dataframe import map_interchange_dtype_to_narwhals_dtype
from narwhals._utils import Version

if TYPE_CHECKING:
    from typing_extensions import Self

    from narwhals.dtypes import DType


class InterchangeSeries:
    _version = Version.V1

    def __init__(self, df: Any) -> None:
        self._native_series = df

    def __narwhals_series__(self) -> Self:
        return self

    def __native_namespace__(self) -> NoReturn:
        msg = (
            "Cannot access native namespace for interchange-level series with unknown backend. "
            "If you would like to see this kind of object supported in Narwhals, please "
            "open a feature request at https://github.com/narwhals-dev/narwhals/issues."
        )
        raise NotImplementedError(msg)

    @property
    def dtype(self) -> DType:
        return map_interchange_dtype_to_narwhals_dtype(self._native_series.dtype)

    @property
    def native(self) -> Any:
        return self._native_series

    def __getattr__(self, attr: str) -> NoReturn:
        msg = (  # pragma: no cover
            f"Attribute {attr} is not supported for interchange-level dataframes.\n\n"
            "Hint: you probably called `nw.from_native` on an object which isn't fully "
            "supported by Narwhals, yet implements `__dataframe__`. If you would like to "
            "see this kind of object supported in Narwhals, please open a feature request "
            "at https://github.com/narwhals-dev/narwhals/issues."
        )
        raise NotImplementedError(msg)
