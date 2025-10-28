from __future__ import annotations

from typing import TYPE_CHECKING, Any

from narwhals._compliant.any_namespace import StringNamespace
from narwhals._pandas_like.utils import PandasLikeSeriesNamespace, is_dtype_pyarrow

if TYPE_CHECKING:
    from narwhals._pandas_like.series import PandasLikeSeries


class PandasLikeSeriesStringNamespace(
    PandasLikeSeriesNamespace, StringNamespace["PandasLikeSeries"]
):
    def len_chars(self) -> PandasLikeSeries:
        return self.with_native(self.native.str.len())

    def replace(
        self, pattern: str, value: str, *, literal: bool, n: int
    ) -> PandasLikeSeries:
        try:
            series = self.native.str.replace(
                pat=pattern, repl=value, n=n, regex=not literal
            )
        except TypeError as e:
            if not isinstance(value, str):
                msg = f"{self.compliant._implementation} backed `.str.replace` only supports str replacement values"
                raise TypeError(msg) from e
            raise
        return self.with_native(series)

    def replace_all(self, pattern: str, value: str, *, literal: bool) -> PandasLikeSeries:
        return self.replace(pattern, value, literal=literal, n=-1)

    def strip_chars(self, characters: str | None) -> PandasLikeSeries:
        return self.with_native(self.native.str.strip(characters))

    def starts_with(self, prefix: str) -> PandasLikeSeries:
        return self.with_native(self.native.str.startswith(prefix))

    def ends_with(self, suffix: str) -> PandasLikeSeries:
        return self.with_native(self.native.str.endswith(suffix))

    def contains(self, pattern: str, *, literal: bool) -> PandasLikeSeries:
        return self.with_native(self.native.str.contains(pat=pattern, regex=not literal))

    def slice(self, offset: int, length: int | None) -> PandasLikeSeries:
        stop = offset + length if length else None
        return self.with_native(self.native.str.slice(start=offset, stop=stop))

    def split(self, by: str) -> PandasLikeSeries:
        implementation = self.implementation
        if not implementation.is_cudf() and not is_dtype_pyarrow(self.native.dtype):
            msg = (
                "This operation requires a pyarrow-backed series. "
                "Please refer to https://narwhals-dev.github.io/narwhals/api-reference/narwhals/#narwhals.maybe_convert_dtypes "
                "and ensure you are using dtype_backend='pyarrow'. "
                "Additionally, make sure you have pandas version 1.5+ and pyarrow installed. "
            )
            raise TypeError(msg)
        return self.with_native(self.native.str.split(pat=by))

    def to_datetime(self, format: str | None) -> PandasLikeSeries:
        # If we know inputs are timezone-aware, we can pass `utc=True` for better performance.
        if format and any(x in format for x in ("%z", "Z")):
            return self.with_native(self._to_datetime(format, utc=True))
        result = self.with_native(self._to_datetime(format, utc=False))
        if (tz := getattr(result.dtype, "time_zone", None)) and tz != "UTC":
            return result.dt.convert_time_zone("UTC")
        return result

    def _to_datetime(self, format: str | None, *, utc: bool) -> Any:
        result = self.implementation.to_native_namespace().to_datetime(
            self.native, format=format, utc=utc
        )
        return (
            result.convert_dtypes(dtype_backend="pyarrow")
            if is_dtype_pyarrow(self.native.dtype)
            else result
        )

    def to_date(self, format: str | None) -> PandasLikeSeries:
        return self.to_datetime(format=format).dt.date()

    def to_uppercase(self) -> PandasLikeSeries:
        return self.with_native(self.native.str.upper())

    def to_lowercase(self) -> PandasLikeSeries:
        return self.with_native(self.native.str.lower())

    def to_titlecase(self) -> PandasLikeSeries:
        return self.with_native(self.native.str.title())

    def zfill(self, width: int) -> PandasLikeSeries:
        return self.with_native(self.native.str.zfill(width))
