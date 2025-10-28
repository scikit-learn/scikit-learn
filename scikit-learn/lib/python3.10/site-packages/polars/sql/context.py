from __future__ import annotations

import contextlib
import re
from typing import (
    TYPE_CHECKING,
    Callable,
    Generic,
    Union,
    overload,
)

from polars._dependencies import _check_for_pandas, _check_for_pyarrow
from polars._dependencies import pandas as pd
from polars._dependencies import pyarrow as pa
from polars._typing import FrameType
from polars._utils.deprecation import deprecate_renamed_parameter
from polars._utils.pycapsule import is_pycapsule
from polars._utils.unstable import issue_unstable_warning
from polars._utils.various import _get_stack_locals, qualified_type_name
from polars._utils.wrap import wrap_ldf
from polars.convert import from_arrow, from_pandas
from polars.dataframe import DataFrame
from polars.lazyframe import LazyFrame
from polars.series import Series

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PySQLContext

if TYPE_CHECKING:
    import sys
    from collections.abc import Collection, Mapping
    from types import TracebackType
    from typing import Any, Final, Literal

    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

    CompatibleFrameType: TypeAlias = Union[
        DataFrame,
        LazyFrame,
        Series,
        pd.DataFrame,
        pd.Series[Any],
        pa.Table,
        pa.RecordBatch,
    ]

__all__ = ["SQLContext"]


def _compatible_frame(obj: Any) -> bool:
    """Check if the object can be converted to DataFrame."""
    return (
        is_pycapsule(obj)
        or isinstance(obj, LazyFrame)
        or (_check_for_pandas(obj) and isinstance(obj, (pd.DataFrame, pd.Series)))
        or (_check_for_pyarrow(obj) and isinstance(obj, (pa.Table, pa.RecordBatch)))
    )


def _ensure_lazyframe(obj: Any) -> LazyFrame:
    """Return LazyFrame from compatible input."""
    if isinstance(obj, (DataFrame, LazyFrame)):
        return obj.lazy()
    elif isinstance(obj, Series):
        return obj.to_frame().lazy()
    elif _check_for_pandas(obj) and isinstance(obj, (pd.DataFrame, pd.Series)):
        if isinstance(frame := from_pandas(obj), Series):
            frame = frame.to_frame()
        return frame.lazy()
    elif is_pycapsule(obj) or (
        _check_for_pyarrow(obj) and isinstance(obj, (pa.Table, pa.RecordBatch))
    ):
        return from_arrow(obj).lazy()  # type: ignore[union-attr]
    else:
        msg = f"unrecognised frame type: {qualified_type_name(obj)}"
        raise ValueError(msg)


def _get_frame_locals(
    *,
    all_compatible: bool,
    n_objects: int | None = None,
    named: str | Collection[str] | Callable[[str], bool] | None = None,
) -> dict[str, Any]:
    """Return compatible frame objects from the local stack."""
    of_type = _compatible_frame if all_compatible else (DataFrame, LazyFrame, Series)
    return _get_stack_locals(of_type=of_type, n_objects=n_objects, named=named)  # type: ignore[arg-type]


class SQLContext(Generic[FrameType]):
    """
    Run SQL queries against DataFrame, LazyFrame, and Series data.

    .. warning::
        This functionality is considered **unstable**, although it is close to being
        considered stable. It may be changed at any point without it being considered
        a breaking change.
    """

    _ctxt: PySQLContext
    _eager_execution: Final[bool]
    _tables_scope_stack: list[set[str]]

    # note: the type-overloaded methods are required to support accurate typing
    # of the frame return from "execute" (which may be DataFrame or LazyFrame),
    # as that is influenced by both the "eager" flag at init-time AND the "eager"
    # flag at query-time (if anyone can find a lighter-weight set of annotations
    # that successfully resolves this, please go for it... ;)

    @overload
    def __init__(
        self: SQLContext[LazyFrame],
        frames: Mapping[str, CompatibleFrameType | None] | None = ...,
        *,
        register_globals: bool | int = ...,
        all_compatible: bool = ...,
        eager: Literal[False] = False,
        **named_frames: CompatibleFrameType | None,
    ) -> None: ...

    @overload
    def __init__(
        self: SQLContext[DataFrame],
        frames: Mapping[str, CompatibleFrameType | None] | None = ...,
        *,
        register_globals: bool | int = ...,
        all_compatible: bool = ...,
        eager: Literal[True],
        **named_frames: CompatibleFrameType | None,
    ) -> None: ...

    @overload
    def __init__(
        self: SQLContext[DataFrame],
        frames: Mapping[str, CompatibleFrameType | None] | None = ...,
        *,
        register_globals: bool | int = ...,
        all_compatible: bool = ...,
        eager: bool,
        **named_frames: CompatibleFrameType | None,
    ) -> None: ...

    @deprecate_renamed_parameter("eager_execution", "eager", version="0.20.31")
    def __init__(
        self,
        frames: Mapping[str, CompatibleFrameType | None] | None = None,
        *,
        register_globals: bool | int = False,
        eager: bool = False,
        **named_frames: CompatibleFrameType | None,
    ) -> None:
        """
        Initialize a new `SQLContext`.

        .. versionchanged:: 0.20.31
            The `eager_execution` parameter was renamed `eager`.

        Parameters
        ----------
        frames
            A `{name:frame, ...}` mapping which can include Polars frames *and*
            pandas DataFrames, Series and pyarrow Table and RecordBatch objects.
        register_globals
            Register compatible objects (polars DataFrame, LazyFrame, and Series) found
            in the globals, automatically mapping their variable name to a table name.
            To register other objects (pandas/pyarrow data) pass them explicitly, or
            call the `execute_global` classmethod. If given an integer then only the
            most recent "n" objects found will be registered.
        eager
            If True, returns execution results as `DataFrame` instead of `LazyFrame`.
            (Note that the query itself is always executed in lazy-mode; this parameter
            impacts whether :meth:`execute` returns an eager or lazy result frame).
        **named_frames
            Named eager/lazy frames, provided as kwargs.

        Examples
        --------
        >>> lf = pl.LazyFrame({"a": [1, 2, 3], "b": ["x", None, "z"]})
        >>> res = pl.SQLContext(frame=lf).execute(
        ...     "SELECT b, a*2 AS two_a FROM frame WHERE b IS NOT NULL"
        ... )
        >>> res.collect()
        shape: (2, 2)
        ┌─────┬───────┐
        │ b   ┆ two_a │
        │ --- ┆ ---   │
        │ str ┆ i64   │
        ╞═════╪═══════╡
        │ x   ┆ 2     │
        │ z   ┆ 6     │
        └─────┴───────┘
        """
        issue_unstable_warning(
            "`SQLContext` is considered **unstable**, although it is close to being considered stable."
        )
        self._ctxt = PySQLContext.new()
        self._eager_execution = eager

        frames = dict(frames or {})
        if register_globals:
            for name, obj in _get_frame_locals(
                all_compatible=False,
            ).items():
                if name not in frames and name not in named_frames:
                    named_frames[name] = obj

        if frames or named_frames:
            frames.update(named_frames)
            self.register_many(frames)

    @overload
    @classmethod
    def execute_global(
        cls, query: str, *, eager: Literal[False] = False
    ) -> LazyFrame: ...

    @overload
    @classmethod
    def execute_global(cls, query: str, *, eager: Literal[True]) -> DataFrame: ...

    @overload
    @classmethod
    def execute_global(cls, query: str, *, eager: bool) -> DataFrame | LazyFrame: ...

    @classmethod
    def execute_global(
        cls, query: str, *, eager: bool = False
    ) -> DataFrame | LazyFrame:
        """
        Immediately execute a SQL query, automatically registering frame globals.

        Notes
        -----
        * This convenience method automatically registers all compatible objects in
          the local stack that are referenced in the query, mapping their variable name
          to a table name. Note that in addition to polars DataFrame, LazyFrame, and
          Series this method *also* registers pandas DataFrame, Series, and pyarrow
          Table and RecordBatch objects.
        * Instead of calling this classmethod you should consider using `pl.sql`,
          which will use this code internally.

        Parameters
        ----------
        query
            A valid SQL query string.
        eager
            If True, returns execution results as `DataFrame` instead of `LazyFrame`.
            (Note that the query itself is always executed in lazy-mode).

        Examples
        --------
        >>> import pandas as pd
        >>> df = pl.LazyFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        >>> df_pandas = pd.DataFrame({"a": [2, 3, 4], "c": [7, 8, 9]})

        Join a polars LazyFrame with a pandas DataFrame (note use of the preferred
        `pl.sql` method, which is equivalent to `SQLContext.execute_global`):

        >>> pl.sql("SELECT df.*, c FROM df JOIN df_pandas USING(a)").collect()
        shape: (2, 3)
        ┌─────┬─────┬─────┐
        │ a   ┆ b   ┆ c   │
        │ --- ┆ --- ┆ --- │
        │ i64 ┆ i64 ┆ i64 │
        ╞═════╪═════╪═════╡
        │ 2   ┆ 5   ┆ 7   │
        │ 3   ┆ 6   ┆ 8   │
        └─────┴─────┴─────┘
        """
        # basic extraction of possible table names from the query, so we don't register
        # unnecessary objects from the globals (ideally we shuoold look to make the
        # underlying `sqlparser-rs` lib parse the query to identify table names)
        q = re.split(r"\bFROM\b", query, maxsplit=1, flags=re.I)
        possible_names = (
            {
                nm.strip('"')
                for nm in re.split(r"\b", q[1])
                if re.match(r'^("[^"]+")$', nm) or nm.isidentifier()
            }
            if len(q) > 1
            else set()
        )
        # get compatible frame objects from the globals, constraining by possible names
        named_frames = _get_frame_locals(all_compatible=True, named=possible_names)
        with cls(frames=named_frames, register_globals=False) as ctx:
            return ctx.execute(query=query, eager=eager)

    def __enter__(self) -> SQLContext[FrameType]:
        """Track currently registered tables on scope entry; supports nested scopes."""
        self._tables_scope_stack = getattr(self, "_tables_scope_stack", [])
        self._tables_scope_stack.append(set(self.tables()))
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        """
        Unregister any tables created within the given scope on context exit.

        See Also
        --------
        unregister
        """
        self.unregister(
            names=(set(self.tables()) - self._tables_scope_stack.pop()),
        )

    def __repr__(self) -> str:
        n_tables = len(self.tables())
        return f"<SQLContext [tables:{n_tables}] at 0x{id(self):x}>"

    # these overloads are necessary to cover the possible permutations
    # of the init-time "eager" param, and the local "eager" param.

    @overload
    def execute(
        self: SQLContext[DataFrame], query: str, *, eager: None = ...
    ) -> DataFrame: ...

    @overload
    def execute(
        self: SQLContext[DataFrame], query: str, *, eager: Literal[False]
    ) -> LazyFrame: ...

    @overload
    def execute(
        self: SQLContext[DataFrame], query: str, *, eager: Literal[True]
    ) -> DataFrame: ...

    @overload
    def execute(
        self: SQLContext[LazyFrame], query: str, *, eager: None = ...
    ) -> LazyFrame: ...

    @overload
    def execute(
        self: SQLContext[LazyFrame], query: str, *, eager: Literal[False]
    ) -> LazyFrame: ...

    @overload
    def execute(
        self: SQLContext[LazyFrame], query: str, *, eager: Literal[True]
    ) -> DataFrame: ...

    @overload
    def execute(
        self, query: str, *, eager: bool | None = ...
    ) -> LazyFrame | DataFrame: ...

    def execute(
        self, query: str, *, eager: bool | None = None
    ) -> LazyFrame | DataFrame:
        """
        Parse the given SQL query and execute it against the registered frame data.

        Parameters
        ----------
        query
            A valid string SQL query.
        eager
            Apply the query eagerly, returning `DataFrame` instead of `LazyFrame`.
            If unset, the value of the init-time "eager" parameter will be used.
            Note that the query itself is always executed in lazy-mode; this
            parameter only impacts the type of the returned frame.

        Examples
        --------
        Declare frame data and register with a SQLContext:

        >>> df = pl.DataFrame(
        ...     data=[
        ...         ("The Godfather", 1972, 6_000_000, 134_821_952, 9.2),
        ...         ("The Dark Knight", 2008, 185_000_000, 533_316_061, 9.0),
        ...         ("Schindler's List", 1993, 22_000_000, 96_067_179, 8.9),
        ...         ("Pulp Fiction", 1994, 8_000_000, 107_930_000, 8.9),
        ...         ("The Shawshank Redemption", 1994, 25_000_000, 28_341_469, 9.3),
        ...     ],
        ...     schema=["title", "release_year", "budget", "gross", "imdb_score"],
        ...     orient="row",
        ... )
        >>> ctx = pl.SQLContext(films=df)

        Execute a SQL query against the registered frame data:

        >>> ctx.execute(
        ...     '''
        ...     SELECT title, release_year, imdb_score
        ...     FROM films
        ...     WHERE release_year > 1990
        ...     ORDER BY imdb_score DESC
        ...     ''',
        ...     eager=True,
        ... )
        shape: (4, 3)
        ┌──────────────────────────┬──────────────┬────────────┐
        │ title                    ┆ release_year ┆ imdb_score │
        │ ---                      ┆ ---          ┆ ---        │
        │ str                      ┆ i64          ┆ f64        │
        ╞══════════════════════════╪══════════════╪════════════╡
        │ The Shawshank Redemption ┆ 1994         ┆ 9.3        │
        │ The Dark Knight          ┆ 2008         ┆ 9.0        │
        │ Schindler's List         ┆ 1993         ┆ 8.9        │
        │ Pulp Fiction             ┆ 1994         ┆ 8.9        │
        └──────────────────────────┴──────────────┴────────────┘

        Execute a GROUP BY query:

        >>> ctx.execute(
        ...     '''
        ...     SELECT
        ...         MAX(release_year / 10) * 10 AS decade,
        ...         SUM(gross) AS total_gross,
        ...         COUNT(title) AS n_films,
        ...     FROM films
        ...     GROUP BY (release_year / 10) -- decade
        ...     ORDER BY total_gross DESC
        ...     ''',
        ...     eager=True,
        ... )
        shape: (3, 3)
        ┌────────┬─────────────┬─────────┐
        │ decade ┆ total_gross ┆ n_films │
        │ ---    ┆ ---         ┆ ---     │
        │ i64    ┆ i64         ┆ u32     │
        ╞════════╪═════════════╪═════════╡
        │ 2000   ┆ 533316061   ┆ 1       │
        │ 1990   ┆ 232338648   ┆ 3       │
        │ 1970   ┆ 134821952   ┆ 1       │
        └────────┴─────────────┴─────────┘
        """
        res = wrap_ldf(self._ctxt.execute(query))
        return res.collect() if (eager or self._eager_execution) else res

    def register(self, name: str, frame: CompatibleFrameType | None) -> Self:
        """
        Register a single frame as a table, using the given name.

        Parameters
        ----------
        name
            Name of the table.
        frame
            eager/lazy frame to associate with this table name.

        See Also
        --------
        register_globals
        register_many
        unregister

        Examples
        --------
        >>> df = pl.DataFrame({"hello": ["world"]})
        >>> ctx = pl.SQLContext()
        >>> ctx.register("frame_data", df).execute("SELECT * FROM frame_data").collect()
        shape: (1, 1)
        ┌───────┐
        │ hello │
        │ ---   │
        │ str   │
        ╞═══════╡
        │ world │
        └───────┘
        """
        frame = LazyFrame() if frame is None else _ensure_lazyframe(frame)
        self._ctxt.register(name, frame._ldf)
        return self

    def register_globals(
        self, n: int | None = None, *, all_compatible: bool = True
    ) -> Self:
        """
        Register all frames (lazy or eager) found in the current globals scope.

        Automatically maps variable names to table names.

        See Also
        --------
        register
        register_many
        unregister

        Parameters
        ----------
        n
            Register only the most recent "n" frames.
        all_compatible
            Control whether we *also* register pandas DataFrame, Series, and
            pyarrow Table and RecordBatch objects. If False, only Polars
            classes are registered with the SQL engine.

        Examples
        --------
        >>> df1 = pl.DataFrame({"a": [1, 2, 3], "b": ["x", None, "z"]})
        >>> df2 = pl.DataFrame({"a": [2, 3, 4], "c": ["t", "w", "v"]})

        Register frames directly from variables found in the current globals scope:

        >>> ctx = pl.SQLContext(register_globals=True)
        >>> ctx.tables()
        ['df1', 'df2']

        Query using the register variable/frame names

        >>> ctx.execute(
        ...     "SELECT a, b, c FROM df1 LEFT JOIN df2 USING (a) ORDER BY a DESC"
        ... ).collect()
        shape: (3, 3)
        ┌─────┬──────┬──────┐
        │ a   ┆ b    ┆ c    │
        │ --- ┆ ---  ┆ ---  │
        │ i64 ┆ str  ┆ str  │
        ╞═════╪══════╪══════╡
        │ 3   ┆ z    ┆ w    │
        │ 2   ┆ null ┆ t    │
        │ 1   ┆ x    ┆ null │
        └─────┴──────┴──────┘
        """
        frames = _get_frame_locals(all_compatible=all_compatible, n_objects=n)
        return self.register_many(frames=frames)

    def register_many(
        self,
        frames: Mapping[str, CompatibleFrameType | None] | None = None,
        **named_frames: CompatibleFrameType | None,
    ) -> Self:
        """
        Register multiple eager/lazy frames as tables, using the associated names.

        Parameters
        ----------
        frames
            A `{name:frame, ...}` mapping.
        **named_frames
            Named eager/lazy frames, provided as kwargs.

        See Also
        --------
        register
        register_globals
        unregister

        Examples
        --------
        >>> lf1 = pl.LazyFrame({"a": [1, 2, 3], "b": ["m", "n", "o"]})
        >>> lf2 = pl.LazyFrame({"a": [2, 3, 4], "c": ["p", "q", "r"]})
        >>> lf3 = pl.LazyFrame({"a": [3, 4, 5], "b": ["s", "t", "u"]})
        >>> lf4 = pl.LazyFrame({"a": [4, 5, 6], "c": ["v", "w", "x"]})

        Register multiple frames at once, either by passing in as a dict...

        >>> ctx = pl.SQLContext().register_many({"tbl1": lf1, "tbl2": lf2})
        >>> ctx.tables()
        ['tbl1', 'tbl2']

        ...or using keyword args:

        >>> ctx.register_many(tbl3=lf3, tbl4=lf4).tables()
        ['tbl1', 'tbl2', 'tbl3', 'tbl4']
        """
        frames = dict(frames or {})
        frames.update(named_frames)
        for name, frame in frames.items():
            self.register(name, frame)
        return self

    def unregister(self, names: str | Collection[str]) -> Self:
        """
        Unregister one or more eager/lazy frames by name.

        Parameters
        ----------
        names
            Names of the tables to unregister.

        Notes
        -----
        You can also control table registration lifetime by using `SQLContext` as a
        context manager; this can often be more useful when such control is wanted:

        >>> df0 = pl.DataFrame({"colx": [0, 1, 2]})
        >>> df1 = pl.DataFrame({"colx": [1, 2, 3]})
        >>> df2 = pl.DataFrame({"colx": [2, 3, 4]})

        Frames registered in-scope are automatically unregistered on scope-exit. Note
        that frames registered on construction will persist through subsequent scopes.

        >>> # register one frame at construction time, and the other two in-scope
        >>> with pl.SQLContext(tbl0=df0) as ctx:
        ...     ctx.register_many(tbl1=df1, tbl2=df2).tables()
        ['tbl0', 'tbl1', 'tbl2']

        After scope exit, none of the tables registered in-scope remain:

        >>> ctx.tables()
        ['tbl0']

        See Also
        --------
        register
        register_globals
        register_many

        Examples
        --------
        >>> df0 = pl.DataFrame({"ints": [9, 8, 7, 6, 5]})
        >>> lf1 = pl.LazyFrame({"text": ["a", "b", "c"]})
        >>> lf2 = pl.LazyFrame({"misc": ["testing1234"]})

        Register with a SQLContext object:

        >>> ctx = pl.SQLContext(test1=df0, test2=lf1, test3=lf2)
        >>> ctx.tables()
        ['test1', 'test2', 'test3']

        Unregister one or more of the tables:

        >>> ctx.unregister(["test1", "test3"]).tables()
        ['test2']
        >>> ctx.unregister("test2").tables()
        []
        """
        if isinstance(names, str):
            names = [names]
        for nm in names:
            self._ctxt.unregister(nm)
        return self

    def tables(self) -> list[str]:
        """
        Return a list of the registered table names.

        Notes
        -----
        The :meth:`tables` method will return the same values as the
        "SHOW TABLES" SQL statement, but as a list instead of a frame.

        Executing as SQL:

        >>> frame_data = pl.DataFrame({"hello": ["world"]})
        >>> ctx = pl.SQLContext(hello_world=frame_data)
        >>> ctx.execute("SHOW TABLES", eager=True)
        shape: (1, 1)
        ┌─────────────┐
        │ name        │
        │ ---         │
        │ str         │
        ╞═════════════╡
        │ hello_world │
        └─────────────┘

        Calling the method:

        >>> ctx.tables()
        ['hello_world']

        Examples
        --------
        >>> df1 = pl.DataFrame({"hello": ["world"]})
        >>> df2 = pl.DataFrame({"foo": ["bar", "baz"]})
        >>> ctx = pl.SQLContext(hello_data=df1, foo_bar=df2)
        >>> ctx.tables()
        ['foo_bar', 'hello_data']
        """
        return sorted(self._ctxt.get_tables())
