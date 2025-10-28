from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


def _parquet_field_overwrites_dict_to_dict_list(
    pqo: dict[str, ParquetFieldOverwrites],
) -> list[dict[str, Any]]:
    children = []
    for name, child in pqo.items():
        if child.name is not None:
            msg = "ParquetFieldOverwrites has both a name in the dictionary and in the overwrites"
            raise ValueError(msg)
        child.name = name
        children.append(_parquet_field_overwrites_to_dict(child))
    return children


def _parquet_field_overwrites_to_dict(pqo: ParquetFieldOverwrites) -> dict[str, Any]:
    d: dict[str, Any] = {}

    # Name
    if pqo.name is not None:
        d["name"] = pqo.name

    # Children
    if pqo.children is not None:
        if isinstance(pqo.children, ParquetFieldOverwrites):
            d["children"] = _parquet_field_overwrites_to_dict(pqo.children)
        elif isinstance(pqo.children, dict):
            d["children"] = _parquet_field_overwrites_dict_to_dict_list(pqo.children)
        elif isinstance(pqo.children, list):
            d["children"] = [_parquet_field_overwrites_to_dict(c) for c in pqo.children]
        else:
            msg = "invalid ParquetFieldOverwrites children type"
            raise TypeError(msg)

    if pqo.field_id is not None:
        d["field_id"] = pqo.field_id

    # Metadata
    if pqo.metadata is not None:
        d["metadata"] = list(pqo.metadata.items())

    if pqo.required is not None:
        d["required"] = pqo.required

    return d


class ParquetFieldOverwrites:
    """
    Write-option overwrites for individual Parquet fields.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.


    Examples
    --------
    >>> lf = pl.LazyFrame(
    ...     {
    ...         "a": [None, 2, 3, 4],
    ...         "b": [[1, 2, 3], [42], [13], [37]],
    ...         "c": [
    ...             {"x": "a", "y": 42},
    ...             {"x": "b", "y": 13},
    ...             {"x": "X", "y": 37},
    ...             {"x": "Y", "y": 15},
    ...         ],
    ...     }
    ... )  # doctest: +SKIP
    >>> lf.sink_parquet(
    ...     "./out/parquet",
    ...     field_overwrites={
    ...         "a": ParquetFieldOverwrites(metadata={"flat_from_polars": "yes"}),
    ...         "b": ParquetFieldOverwrites(
    ...             children=ParquetFieldOverwrites(metadata={"listitem": "yes"}),
    ...             metadata={"list": "true"},
    ...         ),
    ...         "c": ParquetFieldOverwrites(
    ...             children=[
    ...                 ParquetFieldOverwrites(name="x", metadata={"md": "yes"}),
    ...                 ParquetFieldOverwrites(name="y", metadata={"md2": "Yes!"}),
    ...             ],
    ...             metadata={"struct": "true"},
    ...         ),
    ...     },
    ... )  # doctest: +SKIP
    """

    name: None | str  #: Name of the column or field
    children: (
        None
        | ParquetFieldOverwrites
        | list[ParquetFieldOverwrites]
        | dict[str, ParquetFieldOverwrites]
    )  #: Children of the column or field.
    #
    # For flat types (e.g. `Int32`), this should be `None`. For lists, this can be a
    # unnamed `ParquetFieldOverwrites`. For structs, this can be a dict or list of named
    # overwrites.

    field_id: int | None = None  #: The field ID used in the Parquet schema
    metadata: (
        dict[str, None | str] | None
    )  #: Arrow metadata added to the field before writing
    required: bool | None = None  #: Is the field not allowed to have missing values

    def __init__(
        self,
        *,
        name: str | None = None,
        children: (
            None
            | ParquetFieldOverwrites
            | Sequence[ParquetFieldOverwrites]
            | Mapping[str, ParquetFieldOverwrites]
        ) = None,
        field_id: int | None = None,
        metadata: Mapping[str, None | str] | None = None,
        required: bool | None = None,
    ) -> None:
        self.name = name

        if isinstance(children, Mapping):
            self.children = dict(children)
        elif isinstance(children, Sequence):
            self.children = list(children)
        else:
            self.children = children

        self.field_id = field_id
        if isinstance(metadata, Mapping):
            self.metadata = dict(metadata)
        else:
            self.metadata = metadata
        self.required = required
