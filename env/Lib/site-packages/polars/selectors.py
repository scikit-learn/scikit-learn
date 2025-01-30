from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from datetime import timezone
from functools import reduce
from operator import or_
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    NoReturn,
    overload,
)

from polars import functions as F
from polars._utils.parse.expr import _parse_inputs_as_iterable
from polars._utils.various import is_column, re_escape
from polars.datatypes import (
    Binary,
    Boolean,
    Categorical,
    Date,
    Datetime,
    Decimal,
    Duration,
    Object,
    String,
    Time,
    is_polars_dtype,
)
from polars.datatypes.group import (
    FLOAT_DTYPES,
    INTEGER_DTYPES,
    NUMERIC_DTYPES,
    SIGNED_INTEGER_DTYPES,
    TEMPORAL_DTYPES,
    UNSIGNED_INTEGER_DTYPES,
)
from polars.expr import Expr

if TYPE_CHECKING:
    import sys
    from collections.abc import Iterable

    from polars import DataFrame, LazyFrame
    from polars._typing import PolarsDataType, PythonDataType, SelectorType, TimeUnit

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

__all__ = [
    "all",
    "alpha",
    "alphanumeric",
    "binary",
    "boolean",
    "by_dtype",
    "by_index",
    "by_name",
    "categorical",
    "contains",
    "date",
    "datetime",
    "decimal",
    "digit",
    "duration",
    "ends_with",
    "exclude",
    "expand_selector",
    "first",
    "float",
    "integer",
    "is_selector",
    "last",
    "matches",
    "numeric",
    "signed_integer",
    "starts_with",
    "string",
    "temporal",
    "time",
    "unsigned_integer",
]


@overload
def is_selector(obj: _selector_proxy_) -> Literal[True]: ...


@overload
def is_selector(obj: Any) -> Literal[False]: ...


def is_selector(obj: Any) -> bool:
    """
    Indicate whether the given object/expression is a selector.

    Examples
    --------
    >>> from polars.selectors import is_selector
    >>> import polars.selectors as cs
    >>> is_selector(pl.col("colx"))
    False
    >>> is_selector(cs.first() | cs.last())
    True
    """
    # note: don't want to expose the "_selector_proxy_" object
    return isinstance(obj, _selector_proxy_) and hasattr(obj, "_attrs")


# TODO: Don't use this as it collects a schema (can be very expensive for LazyFrame).
#  This should move to IR conversion / Rust.
def expand_selector(
    target: DataFrame | LazyFrame | Mapping[str, PolarsDataType],
    selector: SelectorType | Expr,
    *,
    strict: bool = True,
) -> tuple[str, ...]:
    """
    Expand selector to column names, with respect to a specific frame or target schema.

    .. versionadded:: 0.20.30
        The `strict` parameter was added.

    Parameters
    ----------
    target
        A polars DataFrame, LazyFrame or schema.
    selector
        An arbitrary polars selector (or compound selector).
    strict
        Setting False will additionally allow for a broader range of column selection
        expressions (such as bare columns or use of `.exclude()`) to be expanded, not
        just the dedicated selectors.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "colx": ["a", "b", "c"],
    ...         "coly": [123, 456, 789],
    ...         "colz": [2.0, 5.5, 8.0],
    ...     }
    ... )

    Expand selector with respect to an existing `DataFrame`:

    >>> cs.expand_selector(df, cs.numeric())
    ('coly', 'colz')
    >>> cs.expand_selector(df, cs.first() | cs.last())
    ('colx', 'colz')

    This also works with `LazyFrame`:

    >>> cs.expand_selector(df.lazy(), ~(cs.first() | cs.last()))
    ('coly',)

    Expand selector with respect to a standalone schema:

    >>> schema = {
    ...     "id": pl.Int64,
    ...     "desc": pl.String,
    ...     "count": pl.UInt32,
    ...     "value": pl.Float64,
    ... }
    >>> cs.expand_selector(schema, cs.string() | cs.float())
    ('desc', 'value')

    Allow for non-strict selection expressions (such as those
    including use of an `.exclude()` constraint) to be expanded:

    >>> cs.expand_selector(schema, cs.numeric().exclude("id"), strict=False)
    ('count', 'value')
    """
    if isinstance(target, Mapping):
        from polars.dataframe import DataFrame

        target = DataFrame(schema=target)

    if not (
        is_selector(selector)
        if strict
        else selector.meta.is_column_selection(allow_aliasing=False)
    ):
        msg = f"expected a selector; found {selector!r} instead."
        raise TypeError(msg)

    return tuple(target.select(selector).collect_schema())


# TODO: Don't use this as it collects a schema (can be very expensive for LazyFrame).
#  This should move to IR conversion / Rust.
def _expand_selectors(frame: DataFrame | LazyFrame, *items: Any) -> list[Any]:
    """
    Internal function that expands any selectors to column names in the given input.

    Non-selector values are left as-is.

    Examples
    --------
    >>> from polars.selectors import _expand_selectors
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "colw": ["a", "b"],
    ...         "colx": ["x", "y"],
    ...         "coly": [123, 456],
    ...         "colz": [2.0, 5.5],
    ...     }
    ... )
    >>> _expand_selectors(df, ["colx", cs.numeric()])
    ['colx', 'coly', 'colz']
    >>> _expand_selectors(df, cs.string(), cs.float())
    ['colw', 'colx', 'colz']
    """
    items_iter = _parse_inputs_as_iterable(items)

    expanded: list[Any] = []
    for item in items_iter:
        if is_selector(item):
            selector_cols = expand_selector(frame, item)
            expanded.extend(selector_cols)
        else:
            expanded.append(item)
    return expanded


def _expand_selector_dicts(
    df: DataFrame,
    d: Mapping[Any, Any] | None,
    *,
    expand_keys: bool,
    expand_values: bool,
    tuple_keys: bool = False,
) -> dict[str, Any]:
    """Expand dict key/value selectors into their underlying column names."""
    expanded = {}
    for key, value in (d or {}).items():
        if expand_values and is_selector(value):
            expanded[key] = expand_selector(df, selector=value)
            value = expanded[key]
        if expand_keys and is_selector(key):
            cols = expand_selector(df, selector=key)
            if tuple_keys:
                expanded[cols] = value
            else:
                expanded.update(dict.fromkeys(cols, value))
        else:
            expanded[key] = value
    return expanded


def _combine_as_selector(
    items: (
        str
        | Expr
        | PolarsDataType
        | SelectorType
        | Collection[str | Expr | PolarsDataType | SelectorType]
    ),
    *more_items: str | Expr | PolarsDataType | SelectorType,
) -> SelectorType:
    """Create a combined selector from cols, names, dtypes, and/or other selectors."""
    names, regexes, dtypes, selectors = [], [], [], []  # type: ignore[var-annotated]
    for item in (
        *(
            items
            if isinstance(items, Collection) and not isinstance(items, str)
            else [items]
        ),
        *more_items,
    ):
        if is_selector(item):
            selectors.append(item)
        elif is_polars_dtype(item):
            dtypes.append(item)
        elif isinstance(item, str):
            if item.startswith("^") and item.endswith("$"):
                regexes.append(item)
            else:
                names.append(item)
        elif is_column(item):
            names.append(item.meta.output_name())  # type: ignore[union-attr]
        else:
            msg = f"expected one or more `str`, `DataType` or selector; found {item!r} instead."
            raise TypeError(msg)

    selected = []
    if names:
        selected.append(by_name(*names))
    if dtypes:
        selected.append(by_dtype(*dtypes))
    if regexes:
        selected.append(
            matches(
                regexes[0]
                if len(regexes) > 1
                else "|".join(f"({rx})" for rx in regexes)
            )
        )
    if selectors:
        selected.extend(selectors)

    return reduce(or_, selected)


class _selector_proxy_(Expr):
    """Base column selector expression/proxy."""

    _attrs: dict[str, Any]
    _repr_override: str

    def __init__(
        self,
        expr: Expr,
        name: str,
        parameters: dict[str, Any] | None = None,
    ) -> None:
        self._pyexpr = expr._pyexpr
        self._attrs = {
            "params": parameters,
            "name": name,
        }

    def __hash__(self) -> int:
        # note: this is a suitable hash for selectors (but NOT expressions in general),
        # as the repr is guaranteed to be unique across all selector/param permutations
        return hash(repr(self))

    def __invert__(self) -> Self:
        """Invert the selector."""
        if is_selector(self):
            inverted = all() - self
            inverted._repr_override = f"~{self!r}"
        else:
            inverted = ~self.as_expr()
        return inverted  # type: ignore[return-value]

    def __repr__(self) -> str:
        if not hasattr(self, "_attrs"):
            return repr(self.as_expr())
        elif hasattr(self, "_repr_override"):
            return self._repr_override
        else:
            selector_name, params = self._attrs["name"], self._attrs["params"] or {}
            set_ops = {"and": "&", "or": "|", "sub": "-", "xor": "^"}
            if selector_name in set_ops:
                op = set_ops[selector_name]
                return "({})".format(f" {op} ".join(repr(p) for p in params.values()))
            else:
                str_params = ", ".join(
                    (repr(v)[1:-1] if k.startswith("*") else f"{k}={v!r}")
                    for k, v in params.items()
                ).rstrip(",")
                return f"cs.{selector_name}({str_params})"

    @overload
    def __sub__(self, other: SelectorType) -> SelectorType: ...

    @overload
    def __sub__(self, other: Any) -> SelectorType | Expr: ...

    def __sub__(self, other: Any) -> Expr:
        if is_selector(other):
            return _selector_proxy_(
                self.meta._as_selector().meta._selector_sub(other),
                parameters={"self": self, "other": other},
                name="sub",
            )
        else:
            return self.as_expr().__sub__(other)

    def __rsub__(self, other: Any) -> NoReturn:
        msg = "unsupported operand type(s) for op: ('Expr' - 'Selector')"
        raise TypeError(msg)

    @overload
    def __and__(self, other: SelectorType) -> SelectorType: ...

    @overload
    def __and__(self, other: Any) -> Expr: ...

    def __and__(self, other: Any) -> SelectorType | Expr:
        if is_column(other):
            colname = other.meta.output_name()
            other = by_name(colname)
        if is_selector(other):
            return _selector_proxy_(
                self.meta._as_selector().meta._selector_and(other),
                parameters={"self": self, "other": other},
                name="and",
            )
        else:
            return self.as_expr().__and__(other)

    def __rand__(self, other: Any) -> Expr:
        if is_column(other):
            colname = other.meta.output_name()
            return by_name(colname) & self
        return self.as_expr().__rand__(other)

    @overload
    def __or__(self, other: SelectorType) -> SelectorType: ...

    @overload
    def __or__(self, other: Any) -> Expr: ...

    def __or__(self, other: Any) -> SelectorType | Expr:
        if is_column(other):
            other = by_name(other.meta.output_name())
        if is_selector(other):
            return _selector_proxy_(
                self.meta._as_selector().meta._selector_add(other),
                parameters={"self": self, "other": other},
                name="or",
            )
        else:
            return self.as_expr().__or__(other)

    def __ror__(self, other: Any) -> Expr:
        if is_column(other):
            other = by_name(other.meta.output_name())
        return self.as_expr().__ror__(other)

    @overload
    def __xor__(self, other: SelectorType) -> SelectorType: ...

    @overload
    def __xor__(self, other: Any) -> Expr: ...

    def __xor__(self, other: Any) -> SelectorType | Expr:
        if is_column(other):
            other = by_name(other.meta.output_name())
        if is_selector(other):
            return _selector_proxy_(
                self.meta._as_selector().meta._selector_xor(other),
                parameters={"self": self, "other": other},
                name="xor",
            )
        else:
            return self.as_expr().__or__(other)

    def __rxor__(self, other: Any) -> Expr:
        if is_column(other):
            other = by_name(other.meta.output_name())
        return self.as_expr().__rxor__(other)

    def as_expr(self) -> Expr:
        """
        Materialize the `selector` as a normal expression.

        This ensures that the operators `|`, `&`, `~` and `-`
        are applied on the data and not on the selector sets.

        Examples
        --------
        >>> import polars.selectors as cs
        >>> df = pl.DataFrame(
        ...     {
        ...         "colx": ["aa", "bb", "cc"],
        ...         "coly": [True, False, True],
        ...         "colz": [1, 2, 3],
        ...     }
        ... )

        Inverting the boolean selector will choose the non-boolean columns:

        >>> df.select(~cs.boolean())
        shape: (3, 2)
        ┌──────┬──────┐
        │ colx ┆ colz │
        │ ---  ┆ ---  │
        │ str  ┆ i64  │
        ╞══════╪══════╡
        │ aa   ┆ 1    │
        │ bb   ┆ 2    │
        │ cc   ┆ 3    │
        └──────┴──────┘

        To invert the *values* in the selected boolean columns, we need to
        materialize the selector as a standard expression instead:

        >>> df.select(~cs.boolean().as_expr())
        shape: (3, 1)
        ┌───────┐
        │ coly  │
        │ ---   │
        │ bool  │
        ╞═══════╡
        │ false │
        │ true  │
        │ false │
        └───────┘
        """
        return Expr._from_pyexpr(self._pyexpr)


def _re_string(string: str | Collection[str], *, escape: bool = True) -> str:
    """Return escaped regex, potentially representing multiple string fragments."""
    if isinstance(string, str):
        rx = re_escape(string) if escape else string
    else:
        strings: list[str] = []
        for st in string:
            if isinstance(st, Collection) and not isinstance(st, str):  # type: ignore[redundant-expr]
                strings.extend(st)
            else:
                strings.append(st)
        rx = "|".join((re_escape(x) if escape else x) for x in strings)
    return f"({rx})"


def all() -> SelectorType:
    """
    Select all columns.

    See Also
    --------
    first : Select the first column in the current scope.
    last : Select the last column in the current scope.

    Examples
    --------
    >>> from datetime import date
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "dt": [date(1999, 12, 31), date(2024, 1, 1)],
    ...         "value": [1_234_500, 5_000_555],
    ...     },
    ...     schema_overrides={"value": pl.Int32},
    ... )

    Select all columns, casting them to string:

    >>> df.select(cs.all().cast(pl.String))
    shape: (2, 2)
    ┌────────────┬─────────┐
    │ dt         ┆ value   │
    │ ---        ┆ ---     │
    │ str        ┆ str     │
    ╞════════════╪═════════╡
    │ 1999-12-31 ┆ 1234500 │
    │ 2024-01-01 ┆ 5000555 │
    └────────────┴─────────┘

    Select all columns *except* for those matching the given dtypes:

    >>> df.select(cs.all() - cs.numeric())
    shape: (2, 1)
    ┌────────────┐
    │ dt         │
    │ ---        │
    │ date       │
    ╞════════════╡
    │ 1999-12-31 │
    │ 2024-01-01 │
    └────────────┘
    """
    return _selector_proxy_(F.all(), name="all")


def alpha(ascii_only: bool = False, *, ignore_spaces: bool = False) -> SelectorType:  # noqa: FBT001
    r"""
    Select all columns with alphabetic names (eg: only letters).

    Parameters
    ----------
    ascii_only
        Indicate whether to consider only ASCII alphabetic characters, or the full
        Unicode range of valid letters (accented, idiographic, etc).
    ignore_spaces
        Indicate whether to ignore the presence of spaces in column names; if so,
        only the other (non-space) characters are considered.

    Notes
    -----
    Matching column names cannot contain *any* non-alphabetic characters. Note
    that the definition of "alphabetic" consists of all valid Unicode alphabetic
    characters (`\p{Alphabetic}`) by default; this can be changed by setting
    `ascii_only=True`.

    Examples
    --------
    >>> import polars as pl
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "no1": [100, 200, 300],
    ...         "café": ["espresso", "latte", "mocha"],
    ...         "t or f": [True, False, None],
    ...         "hmm": ["aaa", "bbb", "ccc"],
    ...         "都市": ["東京", "大阪", "京都"],
    ...     }
    ... )

    Select columns with alphabetic names; note that accented
    characters and kanji are recognised as alphabetic here:

    >>> df.select(cs.alpha())
    shape: (3, 3)
    ┌──────────┬─────┬──────┐
    │ café     ┆ hmm ┆ 都市 │
    │ ---      ┆ --- ┆ ---  │
    │ str      ┆ str ┆ str  │
    ╞══════════╪═════╪══════╡
    │ espresso ┆ aaa ┆ 東京 │
    │ latte    ┆ bbb ┆ 大阪 │
    │ mocha    ┆ ccc ┆ 京都 │
    └──────────┴─────┴──────┘

    Constrain the definition of "alphabetic" to ASCII characters only:

    >>> df.select(cs.alpha(ascii_only=True))
    shape: (3, 1)
    ┌─────┐
    │ hmm │
    │ --- │
    │ str │
    ╞═════╡
    │ aaa │
    │ bbb │
    │ ccc │
    └─────┘

    >>> df.select(cs.alpha(ascii_only=True, ignore_spaces=True))
    shape: (3, 2)
    ┌────────┬─────┐
    │ t or f ┆ hmm │
    │ ---    ┆ --- │
    │ bool   ┆ str │
    ╞════════╪═════╡
    │ true   ┆ aaa │
    │ false  ┆ bbb │
    │ null   ┆ ccc │
    └────────┴─────┘

    Select all columns *except* for those with alphabetic names:

    >>> df.select(~cs.alpha())
    shape: (3, 2)
    ┌─────┬────────┐
    │ no1 ┆ t or f │
    │ --- ┆ ---    │
    │ i64 ┆ bool   │
    ╞═════╪════════╡
    │ 100 ┆ true   │
    │ 200 ┆ false  │
    │ 300 ┆ null   │
    └─────┴────────┘

    >>> df.select(~cs.alpha(ignore_spaces=True))
    shape: (3, 1)
    ┌─────┐
    │ no1 │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 100 │
    │ 200 │
    │ 300 │
    └─────┘
    """
    # note that we need to supply a pattern compatible with the *rust* regex crate
    re_alpha = r"a-zA-Z" if ascii_only else r"\p{Alphabetic}"
    re_space = " " if ignore_spaces else ""
    return _selector_proxy_(
        F.col(f"^[{re_alpha}{re_space}]+$"),
        name="alpha",
        parameters={"ascii_only": ascii_only, "ignore_spaces": ignore_spaces},
    )


def alphanumeric(
    ascii_only: bool = False,  # noqa: FBT001
    *,
    ignore_spaces: bool = False,
) -> SelectorType:
    r"""
    Select all columns with alphanumeric names (eg: only letters and the digits 0-9).

    Parameters
    ----------
    ascii_only
        Indicate whether to consider only ASCII alphabetic characters, or the full
        Unicode range of valid letters (accented, idiographic, etc).
    ignore_spaces
        Indicate whether to ignore the presence of spaces in column names; if so,
        only the other (non-space) characters are considered.

    Notes
    -----
    Matching column names cannot contain *any* non-alphabetic or integer characters.
    Note that the definition of "alphabetic" consists of all valid Unicode alphabetic
    characters (`\p{Alphabetic}`) and digit characters (`\d`) by default; this
    can be changed by setting `ascii_only=True`.

    Examples
    --------
    >>> import polars as pl
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "1st_col": [100, 200, 300],
    ...         "flagged": [True, False, True],
    ...         "00prefix": ["01:aa", "02:bb", "03:cc"],
    ...         "last col": ["x", "y", "z"],
    ...     }
    ... )

    Select columns with alphanumeric names:

    >>> df.select(cs.alphanumeric())
    shape: (3, 2)
    ┌─────────┬──────────┐
    │ flagged ┆ 00prefix │
    │ ---     ┆ ---      │
    │ bool    ┆ str      │
    ╞═════════╪══════════╡
    │ true    ┆ 01:aa    │
    │ false   ┆ 02:bb    │
    │ true    ┆ 03:cc    │
    └─────────┴──────────┘

    >>> df.select(cs.alphanumeric(ignore_spaces=True))
    shape: (3, 3)
    ┌─────────┬──────────┬──────────┐
    │ flagged ┆ 00prefix ┆ last col │
    │ ---     ┆ ---      ┆ ---      │
    │ bool    ┆ str      ┆ str      │
    ╞═════════╪══════════╪══════════╡
    │ true    ┆ 01:aa    ┆ x        │
    │ false   ┆ 02:bb    ┆ y        │
    │ true    ┆ 03:cc    ┆ z        │
    └─────────┴──────────┴──────────┘

    Select all columns *except* for those with alphanumeric names:

    >>> df.select(~cs.alphanumeric())
    shape: (3, 2)
    ┌─────────┬──────────┐
    │ 1st_col ┆ last col │
    │ ---     ┆ ---      │
    │ i64     ┆ str      │
    ╞═════════╪══════════╡
    │ 100     ┆ x        │
    │ 200     ┆ y        │
    │ 300     ┆ z        │
    └─────────┴──────────┘

    >>> df.select(~cs.alphanumeric(ignore_spaces=True))
    shape: (3, 1)
    ┌─────────┐
    │ 1st_col │
    │ ---     │
    │ i64     │
    ╞═════════╡
    │ 100     │
    │ 200     │
    │ 300     │
    └─────────┘
    """
    # note that we need to supply patterns compatible with the *rust* regex crate
    re_alpha = r"a-zA-Z" if ascii_only else r"\p{Alphabetic}"
    re_digit = "0-9" if ascii_only else r"\d"
    re_space = " " if ignore_spaces else ""
    return _selector_proxy_(
        F.col(f"^[{re_alpha}{re_digit}{re_space}]+$"),
        name="alphanumeric",
        parameters={"ascii_only": ascii_only, "ignore_spaces": ignore_spaces},
    )


def binary() -> SelectorType:
    """
    Select all binary columns.

    See Also
    --------
    by_dtype : Select all columns matching the given dtype(s).
    string : Select all string columns (optionally including categoricals).

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame({"a": [b"hello"], "b": ["world"], "c": [b"!"], "d": [":)"]})
    >>> df
    shape: (1, 4)
    ┌──────────┬───────┬────────┬─────┐
    │ a        ┆ b     ┆ c      ┆ d   │
    │ ---      ┆ ---   ┆ ---    ┆ --- │
    │ binary   ┆ str   ┆ binary ┆ str │
    ╞══════════╪═══════╪════════╪═════╡
    │ b"hello" ┆ world ┆ b"!"   ┆ :)  │
    └──────────┴───────┴────────┴─────┘

    Select binary columns and export as a dict:

    >>> df.select(cs.binary()).to_dict(as_series=False)
    {'a': [b'hello'], 'c': [b'!']}

    Select all columns *except* for those that are binary:

    >>> df.select(~cs.binary()).to_dict(as_series=False)
    {'b': ['world'], 'd': [':)']}
    """
    return _selector_proxy_(F.col(Binary), name="binary")


def boolean() -> SelectorType:
    """
    Select all boolean columns.

    See Also
    --------
    by_dtype : Select all columns matching the given dtype(s).

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame({"n": range(1, 5)}).with_columns(n_even=pl.col("n") % 2 == 0)
    >>> df
    shape: (4, 2)
    ┌─────┬────────┐
    │ n   ┆ n_even │
    │ --- ┆ ---    │
    │ i64 ┆ bool   │
    ╞═════╪════════╡
    │ 1   ┆ false  │
    │ 2   ┆ true   │
    │ 3   ┆ false  │
    │ 4   ┆ true   │
    └─────┴────────┘

    Select and invert boolean columns:

    >>> df.with_columns(is_odd=cs.boolean().not_())
    shape: (4, 3)
    ┌─────┬────────┬────────┐
    │ n   ┆ n_even ┆ is_odd │
    │ --- ┆ ---    ┆ ---    │
    │ i64 ┆ bool   ┆ bool   │
    ╞═════╪════════╪════════╡
    │ 1   ┆ false  ┆ true   │
    │ 2   ┆ true   ┆ false  │
    │ 3   ┆ false  ┆ true   │
    │ 4   ┆ true   ┆ false  │
    └─────┴────────┴────────┘

    Select all columns *except* for those that are boolean:

    >>> df.select(~cs.boolean())
    shape: (4, 1)
    ┌─────┐
    │ n   │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 1   │
    │ 2   │
    │ 3   │
    │ 4   │
    └─────┘
    """
    return _selector_proxy_(F.col(Boolean), name="boolean")


def by_dtype(
    *dtypes: (
        PolarsDataType
        | PythonDataType
        | Iterable[PolarsDataType]
        | Iterable[PythonDataType]
    ),
) -> SelectorType:
    """
    Select all columns matching the given dtypes.

    See Also
    --------
    by_name : Select all columns matching the given names.
    by_index : Select all columns matching the given indices.

    Examples
    --------
    >>> from datetime import date
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "dt": [date(1999, 12, 31), date(2024, 1, 1), date(2010, 7, 5)],
    ...         "value": [1_234_500, 5_000_555, -4_500_000],
    ...         "other": ["foo", "bar", "foo"],
    ...     }
    ... )

    Select all columns with date or string dtypes:

    >>> df.select(cs.by_dtype(pl.Date, pl.String))
    shape: (3, 2)
    ┌────────────┬───────┐
    │ dt         ┆ other │
    │ ---        ┆ ---   │
    │ date       ┆ str   │
    ╞════════════╪═══════╡
    │ 1999-12-31 ┆ foo   │
    │ 2024-01-01 ┆ bar   │
    │ 2010-07-05 ┆ foo   │
    └────────────┴───────┘

    Select all columns that are not of date or string dtype:

    >>> df.select(~cs.by_dtype(pl.Date, pl.String))
    shape: (3, 1)
    ┌──────────┐
    │ value    │
    │ ---      │
    │ i64      │
    ╞══════════╡
    │ 1234500  │
    │ 5000555  │
    │ -4500000 │
    └──────────┘

    Group by string columns and sum the numeric columns:

    >>> df.group_by(cs.string()).agg(cs.numeric().sum()).sort(by="other")
    shape: (2, 2)
    ┌───────┬──────────┐
    │ other ┆ value    │
    │ ---   ┆ ---      │
    │ str   ┆ i64      │
    ╞═══════╪══════════╡
    │ bar   ┆ 5000555  │
    │ foo   ┆ -3265500 │
    └───────┴──────────┘
    """
    all_dtypes: list[PolarsDataType | PythonDataType] = []
    for tp in dtypes:
        if is_polars_dtype(tp) or isinstance(tp, type):
            all_dtypes.append(tp)
        elif isinstance(tp, Collection):
            for t in tp:
                if not (is_polars_dtype(t) or isinstance(t, type)):
                    msg = f"invalid dtype: {t!r}"
                    raise TypeError(msg)
                all_dtypes.append(t)
        else:
            msg = f"invalid dtype: {tp!r}"
            raise TypeError(msg)

    return _selector_proxy_(
        F.col(all_dtypes), name="by_dtype", parameters={"dtypes": all_dtypes}
    )


def by_index(*indices: int | range | Sequence[int | range]) -> SelectorType:
    """
    Select all columns matching the given indices (or range objects).

    Parameters
    ----------
    *indices
        One or more column indices (or range objects).
        Negative indexing is supported.

    Notes
    -----
    Matching columns are returned in the order in which their indexes
    appear in the selector, not the underlying schema order.

    See Also
    --------
    by_dtype : Select all columns matching the given dtypes.
    by_name : Select all columns matching the given names.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "key": ["abc"],
    ...         **{f"c{i:02}": [0.5 * i] for i in range(100)},
    ...     },
    ... )
    >>> print(df)
    shape: (1, 101)
    ┌─────┬─────┬─────┬─────┬───┬──────┬──────┬──────┬──────┐
    │ key ┆ c00 ┆ c01 ┆ c02 ┆ … ┆ c96  ┆ c97  ┆ c98  ┆ c99  │
    │ --- ┆ --- ┆ --- ┆ --- ┆   ┆ ---  ┆ ---  ┆ ---  ┆ ---  │
    │ str ┆ f64 ┆ f64 ┆ f64 ┆   ┆ f64  ┆ f64  ┆ f64  ┆ f64  │
    ╞═════╪═════╪═════╪═════╪═══╪══════╪══════╪══════╪══════╡
    │ abc ┆ 0.0 ┆ 0.5 ┆ 1.0 ┆ … ┆ 48.0 ┆ 48.5 ┆ 49.0 ┆ 49.5 │
    └─────┴─────┴─────┴─────┴───┴──────┴──────┴──────┴──────┘

    Select columns by index ("key" column and the two first/last columns):

    >>> df.select(cs.by_index(0, 1, 2, -2, -1))
    shape: (1, 5)
    ┌─────┬─────┬─────┬──────┬──────┐
    │ key ┆ c00 ┆ c01 ┆ c98  ┆ c99  │
    │ --- ┆ --- ┆ --- ┆ ---  ┆ ---  │
    │ str ┆ f64 ┆ f64 ┆ f64  ┆ f64  │
    ╞═════╪═════╪═════╪══════╪══════╡
    │ abc ┆ 0.0 ┆ 0.5 ┆ 49.0 ┆ 49.5 │
    └─────┴─────┴─────┴──────┴──────┘

    Select the "key" column and use a `range` object to select various columns.
    Note that you can freely mix and match integer indices and `range` objects:

    >>> df.select(cs.by_index(0, range(1, 101, 20)))
    shape: (1, 6)
    ┌─────┬─────┬──────┬──────┬──────┬──────┐
    │ key ┆ c00 ┆ c20  ┆ c40  ┆ c60  ┆ c80  │
    │ --- ┆ --- ┆ ---  ┆ ---  ┆ ---  ┆ ---  │
    │ str ┆ f64 ┆ f64  ┆ f64  ┆ f64  ┆ f64  │
    ╞═════╪═════╪══════╪══════╪══════╪══════╡
    │ abc ┆ 0.0 ┆ 10.0 ┆ 20.0 ┆ 30.0 ┆ 40.0 │
    └─────┴─────┴──────┴──────┴──────┴──────┘

    >>> df.select(cs.by_index(0, range(101, 0, -25)))
    shape: (1, 5)
    ┌─────┬──────┬──────┬──────┬─────┐
    │ key ┆ c75  ┆ c50  ┆ c25  ┆ c00 │
    │ --- ┆ ---  ┆ ---  ┆ ---  ┆ --- │
    │ str ┆ f64  ┆ f64  ┆ f64  ┆ f64 │
    ╞═════╪══════╪══════╪══════╪═════╡
    │ abc ┆ 37.5 ┆ 25.0 ┆ 12.5 ┆ 0.0 │
    └─────┴──────┴──────┴──────┴─────┘

    Select all columns *except* for the even-indexed ones:

    >>> df.select(~cs.by_index(range(1, 100, 2)))
    shape: (1, 51)
    ┌─────┬─────┬─────┬─────┬───┬──────┬──────┬──────┬──────┐
    │ key ┆ c01 ┆ c03 ┆ c05 ┆ … ┆ c93  ┆ c95  ┆ c97  ┆ c99  │
    │ --- ┆ --- ┆ --- ┆ --- ┆   ┆ ---  ┆ ---  ┆ ---  ┆ ---  │
    │ str ┆ f64 ┆ f64 ┆ f64 ┆   ┆ f64  ┆ f64  ┆ f64  ┆ f64  │
    ╞═════╪═════╪═════╪═════╪═══╪══════╪══════╪══════╪══════╡
    │ abc ┆ 0.5 ┆ 1.5 ┆ 2.5 ┆ … ┆ 46.5 ┆ 47.5 ┆ 48.5 ┆ 49.5 │
    └─────┴─────┴─────┴─────┴───┴──────┴──────┴──────┴──────┘
    """
    all_indices: list[int] = []
    for idx in indices:
        if isinstance(idx, (range, Sequence)):
            all_indices.extend(idx)  # type: ignore[arg-type]
        elif isinstance(idx, int):
            all_indices.append(idx)
        else:
            msg = f"invalid index value: {idx!r}"
            raise TypeError(msg)

    return _selector_proxy_(
        F.nth(*all_indices), name="by_index", parameters={"*indices": indices}
    )


def by_name(*names: str | Collection[str], require_all: bool = True) -> SelectorType:
    """
    Select all columns matching the given names.

    .. versionadded:: 0.20.27
      The `require_all` parameter was added.

    Parameters
    ----------
    *names
        One or more names of columns to select.
    require_all
        Whether to match *all* names (the default) or *any* of the names.

    Notes
    -----
    Matching columns are returned in the order in which they are declared in
    the selector, not the underlying schema order.

    See Also
    --------
    by_dtype : Select all columns matching the given dtypes.
    by_index : Select all columns matching the given indices.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [False, True],
    ...     }
    ... )

    Select columns by name:

    >>> df.select(cs.by_name("foo", "bar"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ foo ┆ bar │
    │ --- ┆ --- │
    │ str ┆ i64 │
    ╞═════╪═════╡
    │ x   ┆ 123 │
    │ y   ┆ 456 │
    └─────┴─────┘

    Match *any* of the given columns by name:

    >>> df.select(cs.by_name("baz", "moose", "foo", "bear", require_all=False))
    shape: (2, 2)
    ┌─────┬─────┐
    │ foo ┆ baz │
    │ --- ┆ --- │
    │ str ┆ f64 │
    ╞═════╪═════╡
    │ x   ┆ 2.0 │
    │ y   ┆ 5.5 │
    └─────┴─────┘

    Match all columns *except* for those given:

    >>> df.select(~cs.by_name("foo", "bar"))
    shape: (2, 2)
    ┌─────┬───────┐
    │ baz ┆ zap   │
    │ --- ┆ ---   │
    │ f64 ┆ bool  │
    ╞═════╪═══════╡
    │ 2.0 ┆ false │
    │ 5.5 ┆ true  │
    └─────┴───────┘
    """
    all_names = []
    for nm in names:
        if isinstance(nm, str):
            all_names.append(nm)
        elif isinstance(nm, Collection):
            for n in nm:
                if not isinstance(n, str):
                    msg = f"invalid name: {n!r}"
                    raise TypeError(msg)
                all_names.append(n)
        else:
            msg = f"invalid name: {nm!r}"
            raise TypeError(msg)

    selector_params: dict[str, Any] = {"*names": all_names}
    match_cols: list[str] | str = all_names
    if not require_all:
        match_cols = f"^({'|'.join(re_escape(nm) for nm in all_names)})$"
        selector_params["require_all"] = require_all

    return _selector_proxy_(
        F.col(match_cols),
        name="by_name",
        parameters=selector_params,
    )


def categorical() -> SelectorType:
    """
    Select all categorical columns.

    See Also
    --------
    by_dtype : Select all columns matching the given dtype(s).
    string : Select all string columns (optionally including categoricals).

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["xx", "yy"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...     },
    ...     schema_overrides={"foo": pl.Categorical},
    ... )

    Select all categorical columns:

    >>> df.select(cs.categorical())
    shape: (2, 1)
    ┌─────┐
    │ foo │
    │ --- │
    │ cat │
    ╞═════╡
    │ xx  │
    │ yy  │
    └─────┘

    Select all columns *except* for those that are categorical:

    >>> df.select(~cs.categorical())
    shape: (2, 2)
    ┌─────┬─────┐
    │ bar ┆ baz │
    │ --- ┆ --- │
    │ i64 ┆ f64 │
    ╞═════╪═════╡
    │ 123 ┆ 2.0 │
    │ 456 ┆ 5.5 │
    └─────┴─────┘
    """
    return _selector_proxy_(F.col(Categorical), name="categorical")


def contains(*substring: str) -> SelectorType:
    """
    Select columns whose names contain the given literal substring(s).

    Parameters
    ----------
    substring
        Substring(s) that matching column names should contain.

    See Also
    --------
    matches : Select all columns that match the given regex pattern.
    ends_with : Select columns that end with the given substring(s).
    starts_with : Select columns that start with the given substring(s).

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [False, True],
    ...     }
    ... )

    Select columns that contain the substring 'ba':

    >>> df.select(cs.contains("ba"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ bar ┆ baz │
    │ --- ┆ --- │
    │ i64 ┆ f64 │
    ╞═════╪═════╡
    │ 123 ┆ 2.0 │
    │ 456 ┆ 5.5 │
    └─────┴─────┘

    Select columns that contain the substring 'ba' or the letter 'z':

    >>> df.select(cs.contains("ba", "z"))
    shape: (2, 3)
    ┌─────┬─────┬───────┐
    │ bar ┆ baz ┆ zap   │
    │ --- ┆ --- ┆ ---   │
    │ i64 ┆ f64 ┆ bool  │
    ╞═════╪═════╪═══════╡
    │ 123 ┆ 2.0 ┆ false │
    │ 456 ┆ 5.5 ┆ true  │
    └─────┴─────┴───────┘

    Select all columns *except* for those that contain the substring 'ba':

    >>> df.select(~cs.contains("ba"))
    shape: (2, 2)
    ┌─────┬───────┐
    │ foo ┆ zap   │
    │ --- ┆ ---   │
    │ str ┆ bool  │
    ╞═════╪═══════╡
    │ x   ┆ false │
    │ y   ┆ true  │
    └─────┴───────┘
    """
    escaped_substring = _re_string(substring)
    raw_params = f"^.*{escaped_substring}.*$"

    return _selector_proxy_(
        F.col(raw_params),
        name="contains",
        parameters={"*substring": escaped_substring},
    )


def date() -> SelectorType:
    """
    Select all date columns.

    See Also
    --------
    datetime : Select all datetime columns, optionally filtering by time unit/zone.
    duration : Select all duration columns, optionally filtering by time unit.
    temporal : Select all temporal columns.
    time : Select all time columns.

    Examples
    --------
    >>> from datetime import date, datetime, time
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "dtm": [datetime(2001, 5, 7, 10, 25), datetime(2031, 12, 31, 0, 30)],
    ...         "dt": [date(1999, 12, 31), date(2024, 8, 9)],
    ...         "tm": [time(0, 0, 0), time(23, 59, 59)],
    ...     },
    ... )

    Select all date columns:

    >>> df.select(cs.date())
    shape: (2, 1)
    ┌────────────┐
    │ dt         │
    │ ---        │
    │ date       │
    ╞════════════╡
    │ 1999-12-31 │
    │ 2024-08-09 │
    └────────────┘

    Select all columns *except* for those that are dates:

    >>> df.select(~cs.date())
    shape: (2, 2)
    ┌─────────────────────┬──────────┐
    │ dtm                 ┆ tm       │
    │ ---                 ┆ ---      │
    │ datetime[μs]        ┆ time     │
    ╞═════════════════════╪══════════╡
    │ 2001-05-07 10:25:00 ┆ 00:00:00 │
    │ 2031-12-31 00:30:00 ┆ 23:59:59 │
    └─────────────────────┴──────────┘
    """
    return _selector_proxy_(F.col(Date), name="date")


def datetime(
    time_unit: TimeUnit | Collection[TimeUnit] | None = None,
    time_zone: (str | timezone | Collection[str | timezone | None] | None) = (
        "*",
        None,
    ),
) -> SelectorType:
    """
    Select all datetime columns, optionally filtering by time unit/zone.

    Parameters
    ----------
    time_unit
        One (or more) of the allowed timeunit precision strings, "ms", "us", and "ns".
        Omit to select columns with any valid timeunit.
    time_zone
        * One or more timezone strings, as defined in zoneinfo (to see valid options
          run `import zoneinfo; zoneinfo.available_timezones()` for a full list).
        * Set `None` to select Datetime columns that do not have a timezone.
        * Set "*" to select Datetime columns that have *any* timezone.

    See Also
    --------
    date : Select all date columns.
    duration : Select all duration columns, optionally filtering by time unit.
    temporal : Select all temporal columns.
    time : Select all time columns.

    Examples
    --------
    >>> from datetime import datetime, date, timezone
    >>> import polars.selectors as cs
    >>> from zoneinfo import ZoneInfo
    >>> tokyo_tz = ZoneInfo("Asia/Tokyo")
    >>> utc_tz = timezone.utc
    >>> df = pl.DataFrame(
    ...     {
    ...         "tstamp_tokyo": [
    ...             datetime(1999, 7, 21, 5, 20, 16, 987654, tzinfo=tokyo_tz),
    ...             datetime(2000, 5, 16, 6, 21, 21, 123465, tzinfo=tokyo_tz),
    ...         ],
    ...         "tstamp_utc": [
    ...             datetime(2023, 4, 10, 12, 14, 16, 999000, tzinfo=utc_tz),
    ...             datetime(2025, 8, 25, 14, 18, 22, 666000, tzinfo=utc_tz),
    ...         ],
    ...         "tstamp": [
    ...             datetime(2000, 11, 20, 18, 12, 16, 600000),
    ...             datetime(2020, 10, 30, 10, 20, 25, 123000),
    ...         ],
    ...         "dt": [date(1999, 12, 31), date(2010, 7, 5)],
    ...     },
    ...     schema_overrides={
    ...         "tstamp_tokyo": pl.Datetime("ns", "Asia/Tokyo"),
    ...         "tstamp_utc": pl.Datetime("us", "UTC"),
    ...     },
    ... )

    Select all datetime columns:

    >>> df.select(cs.datetime())
    shape: (2, 3)
    ┌────────────────────────────────┬─────────────────────────────┬─────────────────────────┐
    │ tstamp_tokyo                   ┆ tstamp_utc                  ┆ tstamp                  │
    │ ---                            ┆ ---                         ┆ ---                     │
    │ datetime[ns, Asia/Tokyo]       ┆ datetime[μs, UTC]           ┆ datetime[μs]            │
    ╞════════════════════════════════╪═════════════════════════════╪═════════════════════════╡
    │ 1999-07-21 05:20:16.987654 JST ┆ 2023-04-10 12:14:16.999 UTC ┆ 2000-11-20 18:12:16.600 │
    │ 2000-05-16 06:21:21.123465 JST ┆ 2025-08-25 14:18:22.666 UTC ┆ 2020-10-30 10:20:25.123 │
    └────────────────────────────────┴─────────────────────────────┴─────────────────────────┘

    Select all datetime columns that have 'us' precision:

    >>> df.select(cs.datetime("us"))
    shape: (2, 2)
    ┌─────────────────────────────┬─────────────────────────┐
    │ tstamp_utc                  ┆ tstamp                  │
    │ ---                         ┆ ---                     │
    │ datetime[μs, UTC]           ┆ datetime[μs]            │
    ╞═════════════════════════════╪═════════════════════════╡
    │ 2023-04-10 12:14:16.999 UTC ┆ 2000-11-20 18:12:16.600 │
    │ 2025-08-25 14:18:22.666 UTC ┆ 2020-10-30 10:20:25.123 │
    └─────────────────────────────┴─────────────────────────┘

    Select all datetime columns that have *any* timezone:

    >>> df.select(cs.datetime(time_zone="*"))
    shape: (2, 2)
    ┌────────────────────────────────┬─────────────────────────────┐
    │ tstamp_tokyo                   ┆ tstamp_utc                  │
    │ ---                            ┆ ---                         │
    │ datetime[ns, Asia/Tokyo]       ┆ datetime[μs, UTC]           │
    ╞════════════════════════════════╪═════════════════════════════╡
    │ 1999-07-21 05:20:16.987654 JST ┆ 2023-04-10 12:14:16.999 UTC │
    │ 2000-05-16 06:21:21.123465 JST ┆ 2025-08-25 14:18:22.666 UTC │
    └────────────────────────────────┴─────────────────────────────┘

    Select all datetime columns that have a *specific* timezone:

    >>> df.select(cs.datetime(time_zone="UTC"))
    shape: (2, 1)
    ┌─────────────────────────────┐
    │ tstamp_utc                  │
    │ ---                         │
    │ datetime[μs, UTC]           │
    ╞═════════════════════════════╡
    │ 2023-04-10 12:14:16.999 UTC │
    │ 2025-08-25 14:18:22.666 UTC │
    └─────────────────────────────┘

    Select all datetime columns that have NO timezone:

    >>> df.select(cs.datetime(time_zone=None))
    shape: (2, 1)
    ┌─────────────────────────┐
    │ tstamp                  │
    │ ---                     │
    │ datetime[μs]            │
    ╞═════════════════════════╡
    │ 2000-11-20 18:12:16.600 │
    │ 2020-10-30 10:20:25.123 │
    └─────────────────────────┘

    Select all columns *except* for datetime columns:

    >>> df.select(~cs.datetime())
    shape: (2, 1)
    ┌────────────┐
    │ dt         │
    │ ---        │
    │ date       │
    ╞════════════╡
    │ 1999-12-31 │
    │ 2010-07-05 │
    └────────────┘
    """  # noqa: W505
    if time_unit is None:
        time_unit = ["ms", "us", "ns"]
    else:
        time_unit = [time_unit] if isinstance(time_unit, str) else list(time_unit)

    if time_zone is None:
        time_zone = [None]
    elif time_zone:
        time_zone = (
            [time_zone] if isinstance(time_zone, (str, timezone)) else list(time_zone)
        )

    datetime_dtypes = [Datetime(tu, tz) for tu in time_unit for tz in time_zone]

    return _selector_proxy_(
        F.col(datetime_dtypes),
        name="datetime",
        parameters={"time_unit": time_unit, "time_zone": time_zone},
    )


def decimal() -> SelectorType:
    """
    Select all decimal columns.

    See Also
    --------
    float : Select all float columns.
    integer : Select all integer columns.
    numeric : Select all numeric columns.

    Examples
    --------
    >>> from decimal import Decimal as D
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [D(123), D(456)],
    ...         "baz": [D("2.0005"), D("-50.5555")],
    ...     },
    ...     schema_overrides={"baz": pl.Decimal(scale=5, precision=10)},
    ... )

    Select all decimal columns:

    >>> df.select(cs.decimal())
    shape: (2, 2)
    ┌──────────────┬───────────────┐
    │ bar          ┆ baz           │
    │ ---          ┆ ---           │
    │ decimal[*,0] ┆ decimal[10,5] │
    ╞══════════════╪═══════════════╡
    │ 123          ┆ 2.00050       │
    │ 456          ┆ -50.55550     │
    └──────────────┴───────────────┘

    Select all columns *except* the decimal ones:

    >>> df.select(~cs.decimal())
    shape: (2, 1)
    ┌─────┐
    │ foo │
    │ --- │
    │ str │
    ╞═════╡
    │ x   │
    │ y   │
    └─────┘
    """
    # TODO: allow explicit selection by scale/precision?
    return _selector_proxy_(F.col(Decimal), name="decimal")


def digit(ascii_only: bool = False) -> SelectorType:  # noqa: FBT001
    r"""
    Select all columns having names consisting only of digits.

    Notes
    -----
    Matching column names cannot contain *any* non-digit characters. Note that the
    definition of "digit" consists of all valid Unicode digit characters (`\d`)
    by default; this can be changed by setting `ascii_only=True`.

    Examples
    --------
    >>> import polars as pl
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "key": ["aaa", "bbb", "aaa", "bbb", "bbb"],
    ...         "year": [2001, 2001, 2025, 2025, 2001],
    ...         "value": [-25, 100, 75, -15, -5],
    ...     }
    ... ).pivot(
    ...     values="value",
    ...     index="key",
    ...     on="year",
    ...     aggregate_function="sum",
    ... )
    >>> print(df)
    shape: (2, 3)
    ┌─────┬──────┬──────┐
    │ key ┆ 2001 ┆ 2025 │
    │ --- ┆ ---  ┆ ---  │
    │ str ┆ i64  ┆ i64  │
    ╞═════╪══════╪══════╡
    │ aaa ┆ -25  ┆ 75   │
    │ bbb ┆ 95   ┆ -15  │
    └─────┴──────┴──────┘

    Select columns with digit names:

    >>> df.select(cs.digit())
    shape: (2, 2)
    ┌──────┬──────┐
    │ 2001 ┆ 2025 │
    │ ---  ┆ ---  │
    │ i64  ┆ i64  │
    ╞══════╪══════╡
    │ -25  ┆ 75   │
    │ 95   ┆ -15  │
    └──────┴──────┘

    Select all columns *except* for those with digit names:

    >>> df.select(~cs.digit())
    shape: (2, 1)
    ┌─────┐
    │ key │
    │ --- │
    │ str │
    ╞═════╡
    │ aaa │
    │ bbb │
    └─────┘

    Demonstrate use of `ascii_only` flag (by default all valid unicode digits
    are considered, but this can be constrained to ascii 0-9):

    >>> df = pl.DataFrame({"१९९९": [1999], "२०७७": [2077], "3000": [3000]})
    >>> df.select(cs.digit())
    shape: (1, 3)
    ┌──────┬──────┬──────┐
    │ १९९९ ┆ २०७७ ┆ 3000 │
    │ ---  ┆ ---  ┆ ---  │
    │ i64  ┆ i64  ┆ i64  │
    ╞══════╪══════╪══════╡
    │ 1999 ┆ 2077 ┆ 3000 │
    └──────┴──────┴──────┘

    >>> df.select(cs.digit(ascii_only=True))
    shape: (1, 1)
    ┌──────┐
    │ 3000 │
    │ ---  │
    │ i64  │
    ╞══════╡
    │ 3000 │
    └──────┘
    """
    re_digit = r"[0-9]" if ascii_only else r"\d"
    return _selector_proxy_(F.col(rf"^{re_digit}+$"), name="digit")


def duration(
    time_unit: TimeUnit | Collection[TimeUnit] | None = None,
) -> SelectorType:
    """
    Select all duration columns, optionally filtering by time unit.

    Parameters
    ----------
    time_unit
        One (or more) of the allowed timeunit precision strings, "ms", "us", and "ns".
        Omit to select columns with any valid timeunit.

    See Also
    --------
    date : Select all date columns.
    datetime : Select all datetime columns, optionally filtering by time unit/zone.
    temporal : Select all temporal columns.
    time : Select all time columns.

    Examples
    --------
    >>> from datetime import date, timedelta
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "dt": [date(2022, 1, 31), date(2025, 7, 5)],
    ...         "td1": [
    ...             timedelta(days=1, milliseconds=123456),
    ...             timedelta(days=1, hours=23, microseconds=987000),
    ...         ],
    ...         "td2": [
    ...             timedelta(days=7, microseconds=456789),
    ...             timedelta(days=14, minutes=999, seconds=59),
    ...         ],
    ...         "td3": [
    ...             timedelta(weeks=4, days=-10, microseconds=999999),
    ...             timedelta(weeks=3, milliseconds=123456, microseconds=1),
    ...         ],
    ...     },
    ...     schema_overrides={
    ...         "td1": pl.Duration("ms"),
    ...         "td2": pl.Duration("us"),
    ...         "td3": pl.Duration("ns"),
    ...     },
    ... )

    Select all duration columns:

    >>> df.select(cs.duration())
    shape: (2, 3)
    ┌────────────────┬─────────────────┬────────────────────┐
    │ td1            ┆ td2             ┆ td3                │
    │ ---            ┆ ---             ┆ ---                │
    │ duration[ms]   ┆ duration[μs]    ┆ duration[ns]       │
    ╞════════════════╪═════════════════╪════════════════════╡
    │ 1d 2m 3s 456ms ┆ 7d 456789µs     ┆ 18d 999999µs       │
    │ 1d 23h 987ms   ┆ 14d 16h 39m 59s ┆ 21d 2m 3s 456001µs │
    └────────────────┴─────────────────┴────────────────────┘

    Select all duration columns that have 'ms' precision:

    >>> df.select(cs.duration("ms"))
    shape: (2, 1)
    ┌────────────────┐
    │ td1            │
    │ ---            │
    │ duration[ms]   │
    ╞════════════════╡
    │ 1d 2m 3s 456ms │
    │ 1d 23h 987ms   │
    └────────────────┘

    Select all duration columns that have 'ms' OR 'ns' precision:

    >>> df.select(cs.duration(["ms", "ns"]))
    shape: (2, 2)
    ┌────────────────┬────────────────────┐
    │ td1            ┆ td3                │
    │ ---            ┆ ---                │
    │ duration[ms]   ┆ duration[ns]       │
    ╞════════════════╪════════════════════╡
    │ 1d 2m 3s 456ms ┆ 18d 999999µs       │
    │ 1d 23h 987ms   ┆ 21d 2m 3s 456001µs │
    └────────────────┴────────────────────┘

    Select all columns *except* for duration columns:

    >>> df.select(~cs.duration())
    shape: (2, 1)
    ┌────────────┐
    │ dt         │
    │ ---        │
    │ date       │
    ╞════════════╡
    │ 2022-01-31 │
    │ 2025-07-05 │
    └────────────┘
    """
    if time_unit is None:
        time_unit = ["ms", "us", "ns"]
    else:
        time_unit = [time_unit] if isinstance(time_unit, str) else list(time_unit)

    duration_dtypes = [Duration(tu) for tu in time_unit]
    return _selector_proxy_(
        F.col(duration_dtypes),
        name="duration",
        parameters={"time_unit": time_unit},
    )


def ends_with(*suffix: str) -> SelectorType:
    """
    Select columns that end with the given substring(s).

    See Also
    --------
    contains : Select columns that contain the given literal substring(s).
    matches : Select all columns that match the given regex pattern.
    starts_with : Select columns that start with the given substring(s).

    Parameters
    ----------
    suffix
        Substring(s) that matching column names should end with.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [False, True],
    ...     }
    ... )

    Select columns that end with the substring 'z':

    >>> df.select(cs.ends_with("z"))
    shape: (2, 1)
    ┌─────┐
    │ baz │
    │ --- │
    │ f64 │
    ╞═════╡
    │ 2.0 │
    │ 5.5 │
    └─────┘

    Select columns that end with *either* the letter 'z' or 'r':

    >>> df.select(cs.ends_with("z", "r"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ bar ┆ baz │
    │ --- ┆ --- │
    │ i64 ┆ f64 │
    ╞═════╪═════╡
    │ 123 ┆ 2.0 │
    │ 456 ┆ 5.5 │
    └─────┴─────┘

    Select all columns *except* for those that end with the substring 'z':

    >>> df.select(~cs.ends_with("z"))
    shape: (2, 3)
    ┌─────┬─────┬───────┐
    │ foo ┆ bar ┆ zap   │
    │ --- ┆ --- ┆ ---   │
    │ str ┆ i64 ┆ bool  │
    ╞═════╪═════╪═══════╡
    │ x   ┆ 123 ┆ false │
    │ y   ┆ 456 ┆ true  │
    └─────┴─────┴───────┘
    """
    escaped_suffix = _re_string(suffix)
    raw_params = f"^.*{escaped_suffix}$"

    return _selector_proxy_(
        F.col(raw_params),
        name="ends_with",
        parameters={"*suffix": escaped_suffix},
    )


def exclude(
    columns: (
        str
        | PolarsDataType
        | SelectorType
        | Expr
        | Collection[str | PolarsDataType | SelectorType | Expr]
    ),
    *more_columns: str | PolarsDataType | SelectorType | Expr,
) -> Expr:
    """
    Select all columns except those matching the given columns, datatypes, or selectors.

    Parameters
    ----------
    columns
        One or more columns (col or name), datatypes, columns, or selectors representing
        the columns to exclude.
    *more_columns
        Additional columns, datatypes, or selectors to exclude, specified as positional
        arguments.

    Notes
    -----
    If excluding a single selector it is simpler to write as `~selector` instead.

    Examples
    --------
    Exclude by column name(s):

    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "aa": [1, 2, 3],
    ...         "ba": ["a", "b", None],
    ...         "cc": [None, 2.5, 1.5],
    ...     }
    ... )
    >>> df.select(cs.exclude("ba", "xx"))
    shape: (3, 2)
    ┌─────┬──────┐
    │ aa  ┆ cc   │
    │ --- ┆ ---  │
    │ i64 ┆ f64  │
    ╞═════╪══════╡
    │ 1   ┆ null │
    │ 2   ┆ 2.5  │
    │ 3   ┆ 1.5  │
    └─────┴──────┘

    Exclude using a column name, a selector, and a dtype:

    >>> df.select(cs.exclude("aa", cs.string(), pl.UInt32))
    shape: (3, 1)
    ┌──────┐
    │ cc   │
    │ ---  │
    │ f64  │
    ╞══════╡
    │ null │
    │ 2.5  │
    │ 1.5  │
    └──────┘
    """
    return ~_combine_as_selector(columns, *more_columns)


def first() -> SelectorType:
    """
    Select the first column in the current scope.

    See Also
    --------
    all : Select all columns.
    last : Select the last column in the current scope.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [0, 1],
    ...     }
    ... )

    Select the first column:

    >>> df.select(cs.first())
    shape: (2, 1)
    ┌─────┐
    │ foo │
    │ --- │
    │ str │
    ╞═════╡
    │ x   │
    │ y   │
    └─────┘

    Select everything  *except* for the first column:

    >>> df.select(~cs.first())
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ bar ┆ baz ┆ zap │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ f64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 123 ┆ 2.0 ┆ 0   │
    │ 456 ┆ 5.5 ┆ 1   │
    └─────┴─────┴─────┘
    """
    return _selector_proxy_(F.first(), name="first")


def float() -> SelectorType:
    """
    Select all float columns.

    See Also
    --------
    integer : Select all integer columns.
    numeric : Select all numeric columns.
    signed_integer : Select all signed integer columns.
    unsigned_integer : Select all unsigned integer columns.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [0.0, 1.0],
    ...     },
    ...     schema_overrides={"baz": pl.Float32, "zap": pl.Float64},
    ... )

    Select all float columns:

    >>> df.select(cs.float())
    shape: (2, 2)
    ┌─────┬─────┐
    │ baz ┆ zap │
    │ --- ┆ --- │
    │ f32 ┆ f64 │
    ╞═════╪═════╡
    │ 2.0 ┆ 0.0 │
    │ 5.5 ┆ 1.0 │
    └─────┴─────┘

    Select all columns *except* for those that are float:

    >>> df.select(~cs.float())
    shape: (2, 2)
    ┌─────┬─────┐
    │ foo ┆ bar │
    │ --- ┆ --- │
    │ str ┆ i64 │
    ╞═════╪═════╡
    │ x   ┆ 123 │
    │ y   ┆ 456 │
    └─────┴─────┘
    """
    return _selector_proxy_(F.col(FLOAT_DTYPES), name="float")


def integer() -> SelectorType:
    """
    Select all integer columns.

    See Also
    --------
    by_dtype : Select columns by dtype.
    float : Select all float columns.
    numeric : Select all numeric columns.
    signed_integer : Select all signed integer columns.
    unsigned_integer : Select all unsigned integer columns.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [0, 1],
    ...     }
    ... )

    Select all integer columns:

    >>> df.select(cs.integer())
    shape: (2, 2)
    ┌─────┬─────┐
    │ bar ┆ zap │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 123 ┆ 0   │
    │ 456 ┆ 1   │
    └─────┴─────┘

    Select all columns *except* for those that are integer :

    >>> df.select(~cs.integer())
    shape: (2, 2)
    ┌─────┬─────┐
    │ foo ┆ baz │
    │ --- ┆ --- │
    │ str ┆ f64 │
    ╞═════╪═════╡
    │ x   ┆ 2.0 │
    │ y   ┆ 5.5 │
    └─────┴─────┘
    """
    return _selector_proxy_(F.col(INTEGER_DTYPES), name="integer")


def signed_integer() -> SelectorType:
    """
    Select all signed integer columns.

    See Also
    --------
    by_dtype : Select columns by dtype.
    float : Select all float columns.
    integer : Select all integer columns.
    numeric : Select all numeric columns.
    unsigned_integer : Select all unsigned integer columns.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": [-123, -456],
    ...         "bar": [3456, 6789],
    ...         "baz": [7654, 4321],
    ...         "zap": ["ab", "cd"],
    ...     },
    ...     schema_overrides={"bar": pl.UInt32, "baz": pl.UInt64},
    ... )

    Select all signed integer columns:

    >>> df.select(cs.signed_integer())
    shape: (2, 1)
    ┌──────┐
    │ foo  │
    │ ---  │
    │ i64  │
    ╞══════╡
    │ -123 │
    │ -456 │
    └──────┘

    >>> df.select(~cs.signed_integer())
    shape: (2, 3)
    ┌──────┬──────┬─────┐
    │ bar  ┆ baz  ┆ zap │
    │ ---  ┆ ---  ┆ --- │
    │ u32  ┆ u64  ┆ str │
    ╞══════╪══════╪═════╡
    │ 3456 ┆ 7654 ┆ ab  │
    │ 6789 ┆ 4321 ┆ cd  │
    └──────┴──────┴─────┘

    Select all integer columns (both signed and unsigned):

    >>> df.select(cs.integer())
    shape: (2, 3)
    ┌──────┬──────┬──────┐
    │ foo  ┆ bar  ┆ baz  │
    │ ---  ┆ ---  ┆ ---  │
    │ i64  ┆ u32  ┆ u64  │
    ╞══════╪══════╪══════╡
    │ -123 ┆ 3456 ┆ 7654 │
    │ -456 ┆ 6789 ┆ 4321 │
    └──────┴──────┴──────┘
    """
    return _selector_proxy_(F.col(SIGNED_INTEGER_DTYPES), name="signed_integer")


def unsigned_integer() -> SelectorType:
    """
    Select all unsigned integer columns.

    See Also
    --------
    by_dtype : Select columns by dtype.
    float : Select all float columns.
    integer : Select all integer columns.
    numeric : Select all numeric columns.
    signed_integer : Select all signed integer columns.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": [-123, -456],
    ...         "bar": [3456, 6789],
    ...         "baz": [7654, 4321],
    ...         "zap": ["ab", "cd"],
    ...     },
    ...     schema_overrides={"bar": pl.UInt32, "baz": pl.UInt64},
    ... )

    Select all unsigned integer columns:

    >>> df.select(cs.unsigned_integer())
    shape: (2, 2)
    ┌──────┬──────┐
    │ bar  ┆ baz  │
    │ ---  ┆ ---  │
    │ u32  ┆ u64  │
    ╞══════╪══════╡
    │ 3456 ┆ 7654 │
    │ 6789 ┆ 4321 │
    └──────┴──────┘

    Select all columns *except* for those that are unsigned integers:

    >>> df.select(~cs.unsigned_integer())
    shape: (2, 2)
    ┌──────┬─────┐
    │ foo  ┆ zap │
    │ ---  ┆ --- │
    │ i64  ┆ str │
    ╞══════╪═════╡
    │ -123 ┆ ab  │
    │ -456 ┆ cd  │
    └──────┴─────┘

    Select all integer columns (both signed and unsigned):

    >>> df.select(cs.integer())
    shape: (2, 3)
    ┌──────┬──────┬──────┐
    │ foo  ┆ bar  ┆ baz  │
    │ ---  ┆ ---  ┆ ---  │
    │ i64  ┆ u32  ┆ u64  │
    ╞══════╪══════╪══════╡
    │ -123 ┆ 3456 ┆ 7654 │
    │ -456 ┆ 6789 ┆ 4321 │
    └──────┴──────┴──────┘
    """
    return _selector_proxy_(F.col(UNSIGNED_INTEGER_DTYPES), name="unsigned_integer")


def last() -> SelectorType:
    """
    Select the last column in the current scope.

    See Also
    --------
    all : Select all columns.
    first : Select the first column in the current scope.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [0, 1],
    ...     }
    ... )

    Select the last column:

    >>> df.select(cs.last())
    shape: (2, 1)
    ┌─────┐
    │ zap │
    │ --- │
    │ i64 │
    ╞═════╡
    │ 0   │
    │ 1   │
    └─────┘

    Select everything  *except* for the last column:

    >>> df.select(~cs.last())
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ foo ┆ bar ┆ baz │
    │ --- ┆ --- ┆ --- │
    │ str ┆ i64 ┆ f64 │
    ╞═════╪═════╪═════╡
    │ x   ┆ 123 ┆ 2.0 │
    │ y   ┆ 456 ┆ 5.5 │
    └─────┴─────┴─────┘
    """
    return _selector_proxy_(F.last(), name="last")


def matches(pattern: str) -> SelectorType:
    """
    Select all columns that match the given regex pattern.

    See Also
    --------
    contains : Select all columns that contain the given substring.
    ends_with : Select all columns that end with the given substring(s).
    starts_with : Select all columns that start with the given substring(s).

    Parameters
    ----------
    pattern
        A valid regular expression pattern, compatible with the `regex crate
        <https://docs.rs/regex/latest/regex/>`_.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [0, 1],
    ...     }
    ... )

    Match column names containing an 'a', preceded by a character that is not 'z':

    >>> df.select(cs.matches("[^z]a"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ bar ┆ baz │
    │ --- ┆ --- │
    │ i64 ┆ f64 │
    ╞═════╪═════╡
    │ 123 ┆ 2.0 │
    │ 456 ┆ 5.5 │
    └─────┴─────┘

    Do not match column names ending in 'R' or 'z' (case-insensitively):

    >>> df.select(~cs.matches(r"(?i)R|z$"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ foo ┆ zap │
    │ --- ┆ --- │
    │ str ┆ i64 │
    ╞═════╪═════╡
    │ x   ┆ 0   │
    │ y   ┆ 1   │
    └─────┴─────┘
    """
    if pattern == ".*":
        return all()
    else:
        if pattern.startswith(".*"):
            pattern = pattern[2:]
        elif pattern.endswith(".*"):
            pattern = pattern[:-2]

        pfx = "^.*" if not pattern.startswith("^") else ""
        sfx = ".*$" if not pattern.endswith("$") else ""
        raw_params = f"{pfx}{pattern}{sfx}"

        return _selector_proxy_(
            F.col(raw_params),
            name="matches",
            parameters={"pattern": pattern},
        )


def numeric() -> SelectorType:
    """
    Select all numeric columns.

    See Also
    --------
    by_dtype : Select columns by dtype.
    float : Select all float columns.
    integer : Select all integer columns.
    signed_integer : Select all signed integer columns.
    unsigned_integer : Select all unsigned integer columns.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": ["x", "y"],
    ...         "bar": [123, 456],
    ...         "baz": [2.0, 5.5],
    ...         "zap": [0, 0],
    ...     },
    ...     schema_overrides={"bar": pl.Int16, "baz": pl.Float32, "zap": pl.UInt8},
    ... )

    Match all numeric columns:

    >>> df.select(cs.numeric())
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ bar ┆ baz ┆ zap │
    │ --- ┆ --- ┆ --- │
    │ i16 ┆ f32 ┆ u8  │
    ╞═════╪═════╪═════╡
    │ 123 ┆ 2.0 ┆ 0   │
    │ 456 ┆ 5.5 ┆ 0   │
    └─────┴─────┴─────┘

    Match all columns *except* for those that are numeric:

    >>> df.select(~cs.numeric())
    shape: (2, 1)
    ┌─────┐
    │ foo │
    │ --- │
    │ str │
    ╞═════╡
    │ x   │
    │ y   │
    └─────┘
    """
    return _selector_proxy_(F.col(NUMERIC_DTYPES), name="numeric")


def object() -> SelectorType:
    """
    Select all object columns.

    See Also
    --------
    by_dtype : Select columns by dtype.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> from uuid import uuid4
    >>> with pl.Config(fmt_str_lengths=36):
    ...     df = pl.DataFrame(
    ...         {
    ...             "idx": [0, 1],
    ...             "uuid_obj": [uuid4(), uuid4()],
    ...             "uuid_str": [str(uuid4()), str(uuid4())],
    ...         },
    ...         schema_overrides={"idx": pl.Int32},
    ...     )
    ...     print(df)  # doctest: +IGNORE_RESULT
    shape: (2, 3)
    ┌─────┬──────────────────────────────────────┬──────────────────────────────────────┐
    │ idx ┆ uuid_obj                             ┆ uuid_str                             │
    │ --- ┆ ---                                  ┆ ---                                  │
    │ i32 ┆ object                               ┆ str                                  │
    ╞═════╪══════════════════════════════════════╪══════════════════════════════════════╡
    │ 0   ┆ 6be063cf-c9c6-43be-878e-e446cfd42981 ┆ acab9fea-c05d-4b91-b639-418004a63f33 │
    │ 1   ┆ 7849d8f9-2cac-48e7-96d3-63cf81c14869 ┆ 28c65415-8b7d-4857-a4ce-300dca14b12b │
    └─────┴──────────────────────────────────────┴──────────────────────────────────────┘

    Select object columns and export as a dict:

    >>> df.select(cs.object()).to_dict(as_series=False)  # doctest: +IGNORE_RESULT
    {
        "uuid_obj": [
            UUID("6be063cf-c9c6-43be-878e-e446cfd42981"),
            UUID("7849d8f9-2cac-48e7-96d3-63cf81c14869"),
        ]
    }

    Select all columns *except* for those that are object and export as dict:

    >>> df.select(~cs.object())  # doctest: +IGNORE_RESULT
    {
        "idx": [0, 1],
        "uuid_str": [
            "acab9fea-c05d-4b91-b639-418004a63f33",
            "28c65415-8b7d-4857-a4ce-300dca14b12b",
        ],
    }
    """  # noqa: W505
    return _selector_proxy_(F.col(Object), name="object")


def starts_with(*prefix: str) -> SelectorType:
    """
    Select columns that start with the given substring(s).

    Parameters
    ----------
    prefix
        Substring(s) that matching column names should start with.

    See Also
    --------
    contains : Select all columns that contain the given substring.
    ends_with : Select all columns that end with the given substring(s).
    matches : Select all columns that match the given regex pattern.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "foo": [1.0, 2.0],
    ...         "bar": [3.0, 4.0],
    ...         "baz": [5, 6],
    ...         "zap": [7, 8],
    ...     }
    ... )

    Match columns starting with a 'b':

    >>> df.select(cs.starts_with("b"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ bar ┆ baz │
    │ --- ┆ --- │
    │ f64 ┆ i64 │
    ╞═════╪═════╡
    │ 3.0 ┆ 5   │
    │ 4.0 ┆ 6   │
    └─────┴─────┘

    Match columns starting with *either* the letter 'b' or 'z':

    >>> df.select(cs.starts_with("b", "z"))
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ bar ┆ baz ┆ zap │
    │ --- ┆ --- ┆ --- │
    │ f64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 3.0 ┆ 5   ┆ 7   │
    │ 4.0 ┆ 6   ┆ 8   │
    └─────┴─────┴─────┘

    Match all columns *except* for those starting with 'b':

    >>> df.select(~cs.starts_with("b"))
    shape: (2, 2)
    ┌─────┬─────┐
    │ foo ┆ zap │
    │ --- ┆ --- │
    │ f64 ┆ i64 │
    ╞═════╪═════╡
    │ 1.0 ┆ 7   │
    │ 2.0 ┆ 8   │
    └─────┴─────┘
    """
    escaped_prefix = _re_string(prefix)
    raw_params = f"^{escaped_prefix}.*$"

    return _selector_proxy_(
        F.col(raw_params),
        name="starts_with",
        parameters={"*prefix": prefix},
    )


def string(*, include_categorical: bool = False) -> SelectorType:
    """
    Select all String (and, optionally, Categorical) string columns.

    See Also
    --------
    binary : Select all binary columns.
    by_dtype : Select all columns matching the given dtype(s).
    categorical: Select all categorical columns.

    Examples
    --------
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "w": ["xx", "yy", "xx", "yy", "xx"],
    ...         "x": [1, 2, 1, 4, -2],
    ...         "y": [3.0, 4.5, 1.0, 2.5, -2.0],
    ...         "z": ["a", "b", "a", "b", "b"],
    ...     },
    ... ).with_columns(
    ...     z=pl.col("z").cast(pl.Categorical("lexical")),
    ... )

    Group by all string columns, sum the numeric columns, then sort by the string cols:

    >>> df.group_by(cs.string()).agg(cs.numeric().sum()).sort(by=cs.string())
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ w   ┆ x   ┆ y   │
    │ --- ┆ --- ┆ --- │
    │ str ┆ i64 ┆ f64 │
    ╞═════╪═════╪═════╡
    │ xx  ┆ 0   ┆ 2.0 │
    │ yy  ┆ 6   ┆ 7.0 │
    └─────┴─────┴─────┘

    Group by all string *and* categorical columns:

    >>> df.group_by(cs.string(include_categorical=True)).agg(cs.numeric().sum()).sort(
    ...     by=cs.string(include_categorical=True)
    ... )
    shape: (3, 4)
    ┌─────┬─────┬─────┬──────┐
    │ w   ┆ z   ┆ x   ┆ y    │
    │ --- ┆ --- ┆ --- ┆ ---  │
    │ str ┆ cat ┆ i64 ┆ f64  │
    ╞═════╪═════╪═════╪══════╡
    │ xx  ┆ a   ┆ 2   ┆ 4.0  │
    │ xx  ┆ b   ┆ -2  ┆ -2.0 │
    │ yy  ┆ b   ┆ 6   ┆ 7.0  │
    └─────┴─────┴─────┴──────┘
    """
    string_dtypes: list[PolarsDataType] = [String]
    if include_categorical:
        string_dtypes.append(Categorical)

    return _selector_proxy_(
        F.col(string_dtypes),
        name="string",
        parameters={"include_categorical": include_categorical},
    )


def temporal() -> SelectorType:
    """
    Select all temporal columns.

    See Also
    --------
    by_dtype : Select all columns matching the given dtype(s).
    date : Select all date columns.
    datetime : Select all datetime columns, optionally filtering by time unit/zone.
    duration : Select all duration columns, optionally filtering by time unit.
    time : Select all time columns.

    Examples
    --------
    >>> from datetime import date, time
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "dt": [date(2021, 1, 1), date(2021, 1, 2)],
    ...         "tm": [time(12, 0, 0), time(20, 30, 45)],
    ...         "value": [1.2345, 2.3456],
    ...     }
    ... )

    Match all temporal columns:

    >>> df.select(cs.temporal())
    shape: (2, 2)
    ┌────────────┬──────────┐
    │ dt         ┆ tm       │
    │ ---        ┆ ---      │
    │ date       ┆ time     │
    ╞════════════╪══════════╡
    │ 2021-01-01 ┆ 12:00:00 │
    │ 2021-01-02 ┆ 20:30:45 │
    └────────────┴──────────┘

    Match all temporal columns *except* for time columns:

    >>> df.select(cs.temporal() - cs.time())
    shape: (2, 1)
    ┌────────────┐
    │ dt         │
    │ ---        │
    │ date       │
    ╞════════════╡
    │ 2021-01-01 │
    │ 2021-01-02 │
    └────────────┘

    Match all columns *except* for temporal columns:

    >>> df.select(~cs.temporal())
    shape: (2, 1)
    ┌────────┐
    │ value  │
    │ ---    │
    │ f64    │
    ╞════════╡
    │ 1.2345 │
    │ 2.3456 │
    └────────┘
    """
    return _selector_proxy_(F.col(TEMPORAL_DTYPES), name="temporal")


def time() -> SelectorType:
    """
    Select all time columns.

    See Also
    --------
    date : Select all date columns.
    datetime : Select all datetime columns, optionally filtering by time unit/zone.
    duration : Select all duration columns, optionally filtering by time unit.
    temporal : Select all temporal columns.

    Examples
    --------
    >>> from datetime import date, datetime, time
    >>> import polars.selectors as cs
    >>> df = pl.DataFrame(
    ...     {
    ...         "dtm": [datetime(2001, 5, 7, 10, 25), datetime(2031, 12, 31, 0, 30)],
    ...         "dt": [date(1999, 12, 31), date(2024, 8, 9)],
    ...         "tm": [time(0, 0, 0), time(23, 59, 59)],
    ...     },
    ... )

    Select all time columns:

    >>> df.select(cs.time())
    shape: (2, 1)
    ┌──────────┐
    │ tm       │
    │ ---      │
    │ time     │
    ╞══════════╡
    │ 00:00:00 │
    │ 23:59:59 │
    └──────────┘

    Select all columns *except* for those that are times:

    >>> df.select(~cs.time())
    shape: (2, 2)
    ┌─────────────────────┬────────────┐
    │ dtm                 ┆ dt         │
    │ ---                 ┆ ---        │
    │ datetime[μs]        ┆ date       │
    ╞═════════════════════╪════════════╡
    │ 2001-05-07 10:25:00 ┆ 1999-12-31 │
    │ 2031-12-31 00:30:00 ┆ 2024-08-09 │
    └─────────────────────┴────────────┘
    """
    return _selector_proxy_(F.col(Time), name="time")
