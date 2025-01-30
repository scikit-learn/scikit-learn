from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars import DataFrame


def _check_arg_is_1byte(
    arg_name: str, arg: str | None, *, can_be_empty: bool = False
) -> None:
    if isinstance(arg, str):
        arg_byte_length = len(arg.encode("utf-8"))
        if can_be_empty:
            if arg_byte_length > 1:
                msg = (
                    f'{arg_name}="{arg}" should be a single byte character or empty,'
                    f" but is {arg_byte_length} bytes long"
                )
                raise ValueError(msg)
        elif arg_byte_length != 1:
            msg = (
                f'{arg_name}="{arg}" should be a single byte character, but is'
                f" {arg_byte_length} bytes long"
            )
            raise ValueError(msg)


def _update_columns(df: DataFrame, new_columns: Sequence[str]) -> DataFrame:
    if df.width > len(new_columns):
        cols = df.columns
        for i, name in enumerate(new_columns):
            cols[i] = name
        new_columns = cols
    df.columns = list(new_columns)
    return df
