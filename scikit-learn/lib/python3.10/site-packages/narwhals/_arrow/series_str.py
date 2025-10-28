from __future__ import annotations

import string
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.utils import ArrowSeriesNamespace, lit, parse_datetime_format
from narwhals._compliant.any_namespace import StringNamespace

if TYPE_CHECKING:
    from narwhals._arrow.series import ArrowSeries
    from narwhals._arrow.typing import Incomplete


class ArrowSeriesStringNamespace(ArrowSeriesNamespace, StringNamespace["ArrowSeries"]):
    def len_chars(self) -> ArrowSeries:
        return self.with_native(pc.utf8_length(self.native))

    def replace(self, pattern: str, value: str, *, literal: bool, n: int) -> ArrowSeries:
        fn = pc.replace_substring if literal else pc.replace_substring_regex
        try:
            arr = fn(self.native, pattern, replacement=value, max_replacements=n)
        except TypeError as e:
            if not isinstance(value, str):
                msg = "PyArrow backed `.str.replace` only supports str replacement values"
                raise TypeError(msg) from e
            raise
        return self.with_native(arr)

    def replace_all(self, pattern: str, value: str, *, literal: bool) -> ArrowSeries:
        try:
            return self.replace(pattern, value, literal=literal, n=-1)
        except TypeError as e:
            if not isinstance(value, str):
                msg = "PyArrow backed `.str.replace_all` only supports str replacement values."
                raise TypeError(msg) from e
            raise

    def strip_chars(self, characters: str | None) -> ArrowSeries:
        return self.with_native(
            pc.utf8_trim(self.native, characters or string.whitespace)
        )

    def starts_with(self, prefix: str) -> ArrowSeries:
        return self.with_native(pc.equal(self.slice(0, len(prefix)).native, lit(prefix)))

    def ends_with(self, suffix: str) -> ArrowSeries:
        return self.with_native(
            pc.equal(self.slice(-len(suffix), None).native, lit(suffix))
        )

    def contains(self, pattern: str, *, literal: bool) -> ArrowSeries:
        check_func = pc.match_substring if literal else pc.match_substring_regex
        return self.with_native(check_func(self.native, pattern))

    def slice(self, offset: int, length: int | None) -> ArrowSeries:
        stop = offset + length if length is not None else None
        return self.with_native(
            pc.utf8_slice_codeunits(self.native, start=offset, stop=stop)
        )

    def split(self, by: str) -> ArrowSeries:
        split_series = pc.split_pattern(self.native, by)  # type: ignore[call-overload]
        return self.with_native(split_series)

    def to_datetime(self, format: str | None) -> ArrowSeries:
        format = parse_datetime_format(self.native) if format is None else format
        timestamp_array = pc.strptime(self.native, format=format, unit="us")
        return self.with_native(timestamp_array)

    def to_date(self, format: str | None) -> ArrowSeries:
        return self.to_datetime(format=format).dt.date()

    def to_uppercase(self) -> ArrowSeries:
        return self.with_native(pc.utf8_upper(self.native))

    def to_lowercase(self) -> ArrowSeries:
        return self.with_native(pc.utf8_lower(self.native))

    def to_titlecase(self) -> ArrowSeries:
        return self.with_native(pc.utf8_title(self.native))

    def zfill(self, width: int) -> ArrowSeries:
        binary_join: Incomplete = pc.binary_join_element_wise
        native = self.native
        hyphen, plus = lit("-"), lit("+")
        first_char, remaining_chars = (
            self.slice(0, 1).native,
            self.slice(1, None).native,
        )

        # Conditions
        less_than_width = pc.less(pc.utf8_length(native), lit(width))
        starts_with_hyphen = pc.equal(first_char, hyphen)
        starts_with_plus = pc.equal(first_char, plus)

        conditions = pc.make_struct(
            pc.and_(starts_with_hyphen, less_than_width),
            pc.and_(starts_with_plus, less_than_width),
            less_than_width,
        )

        # Cases
        padded_remaining_chars = pc.utf8_lpad(remaining_chars, width - 1, padding="0")

        result = pc.case_when(
            conditions,
            binary_join(
                pa.repeat(hyphen, len(native)), padded_remaining_chars, ""
            ),  # starts with hyphen and less than width
            binary_join(
                pa.repeat(plus, len(native)), padded_remaining_chars, ""
            ),  # starts with plus and less than width
            pc.utf8_lpad(native, width=width, padding="0"),  # less than width
            native,
        )
        return self.with_native(result)
