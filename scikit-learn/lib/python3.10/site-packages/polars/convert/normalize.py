# This code is partially forked and adapted from pandas.
# Some parts are distributed under: https://github.com/pandas-dev/pandas/blob/main/LICENSE
from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

from polars._utils.unstable import unstable
from polars.dataframe import DataFrame
from polars.datatypes.constants import N_INFER_DEFAULT

if TYPE_CHECKING:
    from polars._typing import JSONEncoder
    from polars.schema import Schema


def _simple_json_normalize(
    data: dict[Any, Any] | Sequence[dict[Any, Any] | Any],
    separator: str,
    max_level: int,
    encoder: JSONEncoder,
) -> dict[Any, Any] | list[dict[Any, Any]] | Any:
    if max_level > 0:
        # expect dict or list (both are valid JSON objects)
        normalized_json_object = {}
        if isinstance(data, dict):
            normalized_json_object = _normalize_json_ordered(
                data=data,
                separator=separator,
                max_level=max_level,
                encoder=encoder,
            )
        elif isinstance(data, list):
            normalized_json_list = [
                _simple_json_normalize(
                    row,
                    separator=separator,
                    max_level=max_level,
                    encoder=encoder,
                )
                for row in data
            ]
            return normalized_json_list
        return normalized_json_object
    else:
        return data


def _normalize_json(
    data: Any,
    key_string: str,
    normalized_dict: dict[str, Any],
    separator: str,
    max_level: int,
    encoder: JSONEncoder,
) -> dict[str, Any]:
    """
    Main recursive function.

    Designed for the most basic use case of `pl.json_normalize(data)`,
    intended as a performance improvement.

    Parameters
    ----------
    data : Any
        Type dependent on types contained within nested Json
    key_string : str
        New key (with separator(s) in) for data
    normalized_dict : dict
        The new normalized/flattened Json dict
    separator : str, default '.'
        Nested records will generate names separated by sep,
        e.g., for sep='.', { 'foo' : { 'bar' : 0 } } -> foo.bar
    max_level
        recursion depth
    encoder
        Custom JSON encoder; if not given, `json.dumps` is used.
    """
    if isinstance(data, dict):
        if max_level > 0:
            key_root = f"{key_string}{separator}" if key_string else ""
            nested_max_level = max_level - 1

            for key, value in data.items():
                new_key = f"{key_root}{key}" if key_root else key
                _normalize_json(
                    data=value,
                    key_string=new_key,
                    normalized_dict=normalized_dict,
                    separator=separator,
                    max_level=nested_max_level,
                    encoder=encoder,
                )
        else:
            normalized_dict[key_string] = encoder(data)
            return normalized_dict
    else:
        normalized_dict[key_string] = data
    return normalized_dict


def _normalize_json_ordered(
    data: dict[str, Any],
    separator: str,
    max_level: int,
    encoder: JSONEncoder,
) -> dict[str, Any]:
    """
    Order the top level keys and then recursively go to depth.

    Parameters
    ----------
    data
        Deserialized JSON objects (dict or list of dicts)
    separator
        Nested records will generate names separated by sep. e.g.,
        for `separator=".", {"foo": {"bar": 0}}` -> foo.bar.
    max_level
        Max number of levels(depth of dict) to normalize.
    encoder
        Custom JSON encoder; if not given, `json.dumps` is used.

    Returns
    -------
    dict or list of dicts, matching `normalized_json_object`
    """
    top_, nested_data = {}, {}
    for k, v in data.items():
        if isinstance(v, dict):
            nested_data[k] = v
        else:
            top_[k] = v

    nested_ = _normalize_json(
        data=nested_data,
        key_string="",
        normalized_dict={},
        separator=separator,
        max_level=max_level,
        encoder=encoder,
    )
    return {**top_, **nested_}


@unstable()
def json_normalize(
    data: dict[Any, Any] | Sequence[dict[Any, Any] | Any],
    *,
    separator: str = ".",
    max_level: int | None = None,
    schema: Schema | None = None,
    strict: bool = True,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    encoder: JSONEncoder | None = None,
) -> DataFrame:
    """
    Normalize semi-structured deserialized JSON data into a flat table.

    Dictionary objects that will not be unnested/normalized are encoded
    as json string data. Unlike it pandas' counterpart, this function will
    not encode dictionaries as objects at any level.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    Parameters
    ----------
    data
        Deserialized JSON objects.
    separator
        Nested records will generate names separated by sep. e.g.,
        for `separator=".", {"foo": {"bar": 0}}` -> foo.bar.
    max_level
        Max number of levels(depth of dict) to normalize.
        If None, normalizes all levels.
    schema
        Overwrite the `Schema` when the normalized data is passed to
        the `DataFrame` constructor.
    strict
        Whether Polars should be strict when constructing the DataFrame.
    infer_schema_length
        Number of rows to take into consideration to determine the schema.
    encoder
        Custom JSON encoder function; if not given, `json.dumps` is used.

    Examples
    --------
    >>> data = [
    ...     {
    ...         "id": 1,
    ...         "name": "Cole Volk",
    ...         "fitness": {"height": 180, "weight": 85},
    ...     },
    ...     {
    ...         "id": 2,
    ...         "name": "Faye Raker",
    ...         "fitness": {"height": 155, "weight": 58},
    ...     },
    ...     {
    ...         "name": "Mark Reg",
    ...         "fitness": {"height": 170, "weight": 78},
    ...     },
    ... ]
    >>> pl.json_normalize(data, max_level=1)
    shape: (3, 4)
    ┌──────┬────────────┬────────────────┬────────────────┐
    │ id   ┆ name       ┆ fitness.height ┆ fitness.weight │
    │ ---  ┆ ---        ┆ ---            ┆ ---            │
    │ i64  ┆ str        ┆ i64            ┆ i64            │
    ╞══════╪════════════╪════════════════╪════════════════╡
    │ 1    ┆ Cole Volk  ┆ 180            ┆ 85             │
    │ 2    ┆ Faye Raker ┆ 155            ┆ 58             │
    │ null ┆ Mark Reg   ┆ 170            ┆ 78             │
    └──────┴────────────┴────────────────┴────────────────┘

    Normalize to a specific depth, using a custom JSON encoder
    (note that `orson.dumps` encodes to bytes, not str).

    >>> import orjson
    >>> pl.json_normalize(data, max_level=0, encoder=orjson.dumps)
    shape: (3, 3)
    ┌──────┬────────────┬───────────────────────────────┐
    │ id   ┆ name       ┆ fitness                       │
    │ ---  ┆ ---        ┆ ---                           │
    │ i64  ┆ str        ┆ binary                        │
    ╞══════╪════════════╪═══════════════════════════════╡
    │ 1    ┆ Cole Volk  ┆ b"{"height":180,"weight":85}" │
    │ 2    ┆ Faye Raker ┆ b"{"height":155,"weight":58}" │
    │ null ┆ Mark Reg   ┆ b"{"height":170,"weight":78}" │
    └──────┴────────────┴───────────────────────────────┘
    """
    if max_level is None:
        max_level = 1 << 32  # eg: u32
    max_level += 1

    if isinstance(data, Sequence) and len(data) == 0:
        return DataFrame(schema=schema)
    elif isinstance(data, Mapping):
        data = [data]
    elif isinstance(data, Iterable) and not isinstance(data, str):  # type: ignore[redundant-expr]
        data = list(data)
    else:
        msg = "expected list or dict of objects"
        raise ValueError(msg)

    if encoder is None:
        encoder = json.dumps

    return DataFrame(
        _simple_json_normalize(
            data,
            separator=separator,
            max_level=max_level,
            encoder=encoder,
        ),
        schema=schema,
        strict=strict,
        infer_schema_length=infer_schema_length,
    )
