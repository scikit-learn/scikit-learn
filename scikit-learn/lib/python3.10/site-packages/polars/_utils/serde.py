"""Utility for serializing Polars objects."""

from __future__ import annotations

from io import BytesIO, StringIO
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Literal, overload

from polars._utils.various import normalize_filepath

if TYPE_CHECKING:
    from io import IOBase

    from polars._typing import SerializationFormat


@overload
def serialize_polars_object(
    serializer: Callable[[IOBase | str], None], file: None, format: Literal["binary"]
) -> bytes: ...
@overload
def serialize_polars_object(
    serializer: Callable[[IOBase | str], None], file: None, format: Literal["json"]
) -> str: ...
@overload
def serialize_polars_object(
    serializer: Callable[[IOBase | str], None],
    file: IOBase | str | Path,
    format: SerializationFormat,
) -> None: ...


def serialize_polars_object(
    serializer: Callable[[IOBase | str], None],
    file: IOBase | str | Path | None,
    format: SerializationFormat,
) -> bytes | str | None:
    """Serialize a Polars object (DataFrame/LazyFrame/Expr)."""

    def serialize_to_bytes() -> bytes:
        with BytesIO() as buf:
            serializer(buf)
            serialized = buf.getvalue()
        return serialized

    if file is None:
        serialized = serialize_to_bytes()
        return serialized.decode() if format == "json" else serialized
    elif isinstance(file, StringIO):
        serialized_str = serialize_to_bytes().decode()
        file.write(serialized_str)
        return None
    elif isinstance(file, BytesIO):
        serialized = serialize_to_bytes()
        file.write(serialized)
        return None
    elif isinstance(file, (str, Path)):
        file = normalize_filepath(file)
        serializer(file)
        return None
    else:
        serializer(file)
        return None
