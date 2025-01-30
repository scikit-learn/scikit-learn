# This code is partially forked and adapted from pandas.
# Some parts are distributed under: https://github.com/pandas-dev/pandas/blob/main/LICENSE
from __future__ import annotations

import json
from collections import abc
from typing import TYPE_CHECKING, Any

from polars._utils.unstable import unstable
from polars.dataframe import DataFrame
from polars.datatypes.constants import N_INFER_DEFAULT

if TYPE_CHECKING:
    from collections.abc import Sequence

    from polars.schema import Schema


def _simple_json_normalize(
    data: dict[Any, Any] | Sequence[dict[Any, Any] | Any],
    separator: str,
    max_level: int,
) -> dict[Any, Any] | list[dict[Any, Any]] | Any:
    if max_level > 0:
        normalized_json_object = {}
        # expect a dictionary, as most jsons are. However, lists are perfectly valid
        if isinstance(data, dict):
            normalized_json_object = _normalize_json_ordered(
                data=data, separator=separator, max_level=max_level
            )
        elif isinstance(data, list):
            normalized_json_list = [
                _simple_json_normalize(row, separator=separator, max_level=max_level)
                for row in data
            ]
            return normalized_json_list
        return normalized_json_object
    else:
        return data


def _normalize_json_ordered(
    data: dict[str, Any], separator: str, max_level: int
) -> dict[str, Any]:
    """
    Order the top level keys and then recursively go to depth.

    Parameters
    ----------
    data
        dict or list of dicts
    separator
        str, default '.'
        Nested records will generate names separated by sep,
        e.g., for sep='.', { 'foo' : { 'bar' : 0 } } -> foo.bar
    max_level
        max recursing level

    Returns
    -------
    dict or list of dicts, matching `normalized_json_object`
    """
    top_dict_ = {k: v for k, v in data.items() if not isinstance(v, dict)}
    nested_dict_ = normalize_json(
        data={k: v for k, v in data.items() if isinstance(v, dict)},
        key_string="",
        normalized_dict={},
        separator=separator,
        max_level=max_level,
    )
    return {**top_dict_, **nested_dict_}


@unstable()
def json_normalize(
    data: dict[Any, Any] | Sequence[dict[Any, Any] | Any],
    *,
    separator: str = ".",
    max_level: int | None = None,
    schema: Schema | None = None,
    strict: bool = True,
    infer_schema_length: int | None = N_INFER_DEFAULT,
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

    Examples
    --------
    >>> data = [
    ...     {
    ...         "id": 1,
    ...         "name": "Cole Volk",
    ...         "fitness": {"height": 130, "weight": 60},
    ...     },
    ...     {"name": "Mark Reg", "fitness": {"height": 130, "weight": 60}},
    ...     {
    ...         "id": 2,
    ...         "name": "Faye Raker",
    ...         "fitness": {"height": 130, "weight": 60},
    ...     },
    ... ]
    >>> pl.json_normalize(data, max_level=1)
    shape: (3, 4)
    ┌──────┬────────────┬────────────────┬────────────────┐
    │ id   ┆ name       ┆ fitness.height ┆ fitness.weight │
    │ ---  ┆ ---        ┆ ---            ┆ ---            │
    │ i64  ┆ str        ┆ i64            ┆ i64            │
    ╞══════╪════════════╪════════════════╪════════════════╡
    │ 1    ┆ Cole Volk  ┆ 130            ┆ 60             │
    │ null ┆ Mark Reg   ┆ 130            ┆ 60             │
    │ 2    ┆ Faye Raker ┆ 130            ┆ 60             │
    └──────┴────────────┴────────────────┴────────────────┘
    >>> pl.json_normalize(data, max_level=0)
    shape: (3, 3)
    ┌──────┬────────────┬───────────────────────────────┐
    │ id   ┆ name       ┆ fitness                       │
    │ ---  ┆ ---        ┆ ---                           │
    │ i64  ┆ str        ┆ str                           │
    ╞══════╪════════════╪═══════════════════════════════╡
    │ 1    ┆ Cole Volk  ┆ {"height": 130, "weight": 60} │
    │ null ┆ Mark Reg   ┆ {"height": 130, "weight": 60} │
    │ 2    ┆ Faye Raker ┆ {"height": 130, "weight": 60} │
    └──────┴────────────┴───────────────────────────────┘

    """
    if max_level is None:
        max_level = 1 << 32
    max_level += 1
    if isinstance(data, list) and len(data) == 0:
        return DataFrame()
    elif isinstance(data, dict):
        data = [data]
    elif isinstance(data, abc.Iterable) and not isinstance(data, str):  # type: ignore[redundant-expr]
        data = list(data)
    else:
        msg = "expected list of objects"
        raise ValueError(msg)
    return DataFrame(
        _simple_json_normalize(data, separator=separator, max_level=max_level),
        schema=schema,
        strict=strict,
        infer_schema_length=infer_schema_length,
    )


def normalize_json(
    data: Any,
    key_string: str,
    normalized_dict: dict[str, Any],
    separator: str,
    max_level: int,
) -> dict[str, Any]:
    """
    Main recursive function.

    Designed for the most basic use case of pl.json_normalize(data)
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
    """
    if isinstance(data, dict):
        if max_level > 0:
            for key, value in data.items():
                new_key = f"{key_string}{separator}{key}"

                if not key_string:
                    new_key = new_key.removeprefix(separator)

                normalize_json(
                    data=value,
                    key_string=new_key,
                    normalized_dict=normalized_dict,
                    separator=separator,
                    max_level=max_level - 1,
                )
        else:
            normalized_dict[key_string] = json.dumps(data)
            return normalized_dict
    else:
        normalized_dict[key_string] = data
    return normalized_dict
