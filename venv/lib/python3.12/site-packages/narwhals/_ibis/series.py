from __future__ import annotations

from typing import TYPE_CHECKING, Any, NoReturn

from narwhals._ibis.utils import native_to_narwhals_dtype
from narwhals.dependencies import get_ibis

if TYPE_CHECKING:
    from types import ModuleType

    from typing_extensions import Self

    from narwhals._utils import Version
    from narwhals.dtypes import DType


class IbisInterchangeSeries:
    def __init__(self, df: Any, version: Version) -> None:
        self._native_series = df
        self._version = version

    def __narwhals_series__(self) -> Self:
        return self

    def __native_namespace__(self) -> ModuleType:
        return get_ibis()

    @property
    def dtype(self) -> DType:
        return native_to_narwhals_dtype(
            self._native_series.schema().types[0], self._version
        )

    def __getattr__(self, attr: str) -> NoReturn:
        msg = (
            f"Attribute {attr} is not supported for interchange-level dataframes.\n\n"
            "If you would like to see this kind of object better supported in "
            "Narwhals, please open a feature request "
            "at https://github.com/narwhals-dev/narwhals/issues."
        )
        raise NotImplementedError(msg)
