from __future__ import annotations

import inspect
import os
import re
import sys
import warnings
from collections import Counter
from collections.abc import (
    Collection,
    Generator,
    Iterable,
    MappingView,
    Sequence,
    Sized,
)
from enum import Enum
from io import BytesIO
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    TypeVar,
    overload,
)

import polars as pl
from polars import functions as F
from polars._dependencies import _check_for_numpy, import_optional, subprocess
from polars._dependencies import numpy as np
from polars.datatypes import (
    Boolean,
    Date,
    Datetime,
    Decimal,
    Duration,
    Int64,
    String,
    Time,
)
from polars.datatypes.group import FLOAT_DTYPES, INTEGER_DTYPES

if TYPE_CHECKING:
    from collections.abc import Iterator, MutableMapping, Reversible

    from polars import DataFrame, Expr
    from polars._typing import PolarsDataType, SizeUnit

    if sys.version_info >= (3, 13):
        from typing import TypeIs
    else:
        from typing_extensions import TypeIs

    if sys.version_info >= (3, 10):
        from typing import ParamSpec, TypeGuard
    else:
        from typing_extensions import ParamSpec, TypeGuard

    P = ParamSpec("P")
    T = TypeVar("T")

# note: reversed views don't match as instances of MappingView
if sys.version_info >= (3, 11):
    _views: list[Reversible[Any]] = [{}.keys(), {}.values(), {}.items()]
    _reverse_mapping_views = tuple(type(reversed(view)) for view in _views)


def _process_null_values(
    null_values: None | str | Sequence[str] | dict[str, str] = None,
) -> None | str | Sequence[str] | list[tuple[str, str]]:
    if isinstance(null_values, dict):
        return list(null_values.items())
    else:
        return null_values


def _is_generator(val: object | Iterator[T]) -> TypeIs[Iterator[T]]:
    return (
        (isinstance(val, (Generator, Iterable)) and not isinstance(val, Sized))
        or isinstance(val, MappingView)
        or (sys.version_info >= (3, 11) and isinstance(val, _reverse_mapping_views))
    )


def _is_iterable_of(val: Iterable[object], eltype: type | tuple[type, ...]) -> bool:
    """Check whether the given iterable is of the given type(s)."""
    return all(isinstance(x, eltype) for x in val)


def is_path_or_str_sequence(
    val: object, *, allow_str: bool = False, include_series: bool = False
) -> TypeGuard[Sequence[str | Path]]:
    """
    Check that `val` is a sequence of strings or paths.

    Note that a single string is a sequence of strings by definition, use
    `allow_str=False` to return False on a single string.
    """
    if allow_str is False and isinstance(val, str):
        return False
    elif _check_for_numpy(val) and isinstance(val, np.ndarray):
        return np.issubdtype(val.dtype, np.str_)
    elif include_series and isinstance(val, pl.Series):
        return val.dtype == pl.String
    return (
        not isinstance(val, bytes)
        and isinstance(val, Sequence)
        and _is_iterable_of(val, (Path, str))
    )


def is_bool_sequence(
    val: object, *, include_series: bool = False
) -> TypeGuard[Sequence[bool]]:
    """Check whether the given sequence is a sequence of booleans."""
    if _check_for_numpy(val) and isinstance(val, np.ndarray):
        return val.dtype == np.bool_
    elif include_series and isinstance(val, pl.Series):
        return val.dtype == pl.Boolean
    return isinstance(val, Sequence) and _is_iterable_of(val, bool)


def is_int_sequence(
    val: object, *, include_series: bool = False
) -> TypeGuard[Sequence[int]]:
    """Check whether the given sequence is a sequence of integers."""
    if _check_for_numpy(val) and isinstance(val, np.ndarray):
        return np.issubdtype(val.dtype, np.integer)
    elif include_series and isinstance(val, pl.Series):
        return val.dtype.is_integer()
    return isinstance(val, Sequence) and _is_iterable_of(val, int)


def is_sequence(
    val: object, *, include_series: bool = False
) -> TypeGuard[Sequence[Any]]:
    """Check whether the given input is a numpy array or python sequence."""
    return (_check_for_numpy(val) and isinstance(val, np.ndarray)) or (
        isinstance(val, (pl.Series, Sequence) if include_series else Sequence)
        and not isinstance(val, str)
    )


def is_str_sequence(
    val: object, *, allow_str: bool = False, include_series: bool = False
) -> TypeGuard[Sequence[str]]:
    """
    Check that `val` is a sequence of strings.

    Note that a single string is a sequence of strings by definition, use
    `allow_str=False` to return False on a single string.
    """
    if allow_str is False and isinstance(val, str):
        return False
    elif _check_for_numpy(val) and isinstance(val, np.ndarray):
        return np.issubdtype(val.dtype, np.str_)
    elif include_series and isinstance(val, pl.Series):
        return val.dtype == pl.String
    return isinstance(val, Sequence) and _is_iterable_of(val, str)


def is_column(obj: Any) -> bool:
    """Indicate if the given object is a basic/unaliased column."""
    from polars.expr import Expr

    return isinstance(obj, Expr) and obj.meta.is_column()


def warn_null_comparison(obj: Any) -> None:
    """Warn for possibly unintentional comparisons with None."""
    if obj is None:
        warnings.warn(
            "Comparisons with None always result in null. Consider using `.is_null()` or `.is_not_null()`.",
            UserWarning,
            stacklevel=find_stacklevel(),
        )


def range_to_series(
    name: str, rng: range, dtype: PolarsDataType | None = None
) -> pl.Series:
    """Fast conversion of the given range to a Series."""
    dtype = dtype or Int64
    if dtype.is_integer():
        range = F.int_range(  # type: ignore[call-overload]
            start=rng.start, end=rng.stop, step=rng.step, dtype=dtype, eager=True
        )
    else:
        range = F.int_range(
            start=rng.start, end=rng.stop, step=rng.step, eager=True
        ).cast(dtype)
    return range.alias(name)


def range_to_slice(rng: range) -> slice:
    """Return the given range as an equivalent slice."""
    return slice(rng.start, rng.stop, rng.step)


def _in_notebook() -> bool:
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


def _in_marimo_notebook() -> bool:
    try:
        import marimo as mo

        return mo.running_in_notebook()  # pragma: no cover
    except ImportError:
        return False


def arrlen(obj: Any) -> int | None:
    """Return length of (non-string/dict) sequence; returns None for non-sequences."""
    try:
        return None if isinstance(obj, (str, dict)) else len(obj)
    except TypeError:
        return None


def normalize_filepath(path: str | Path, *, check_not_directory: bool = True) -> str:
    """Create a string path, expanding the home directory if present."""
    # don't use pathlib here as it modifies slashes (s3:// -> s3:/)
    path = os.path.expanduser(path)  # noqa: PTH111
    if (
        check_not_directory
        and os.path.exists(path)  # noqa: PTH110
        and os.path.isdir(path)  # noqa: PTH112
    ):
        msg = f"expected a file path; {path!r} is a directory"
        raise IsADirectoryError(msg)
    return path


def parse_version(version: Sequence[str | int]) -> tuple[int, ...]:
    """Simple version parser; split into a tuple of ints for comparison."""
    if isinstance(version, str):
        version = version.split(".")
    return tuple(int(re.sub(r"\D", "", str(v))) for v in version)


def ordered_unique(values: Sequence[Any]) -> list[Any]:
    """Return unique list of sequence values, maintaining their order of appearance."""
    seen: set[Any] = set()
    add_ = seen.add
    return [v for v in values if not (v in seen or add_(v))]


def deduplicate_names(names: Iterable[str]) -> list[str]:
    """Ensure name uniqueness by appending a counter to subsequent duplicates."""
    seen: MutableMapping[str, int] = Counter()
    deduped = []
    for nm in names:
        deduped.append(f"{nm}{seen[nm] - 1}" if nm in seen else nm)
        seen[nm] += 1
    return deduped


@overload
def scale_bytes(sz: int, unit: SizeUnit) -> int | float: ...


@overload
def scale_bytes(sz: Expr, unit: SizeUnit) -> Expr: ...


def scale_bytes(sz: int | Expr, unit: SizeUnit) -> int | float | Expr:
    """Scale size in bytes to other size units (eg: "kb", "mb", "gb", "tb")."""
    if unit in {"b", "bytes"}:
        return sz
    elif unit in {"kb", "kilobytes"}:
        return sz / 1024
    elif unit in {"mb", "megabytes"}:
        return sz / 1024**2
    elif unit in {"gb", "gigabytes"}:
        return sz / 1024**3
    elif unit in {"tb", "terabytes"}:
        return sz / 1024**4
    else:
        msg = f"`unit` must be one of {{'b', 'kb', 'mb', 'gb', 'tb'}}, got {unit!r}"
        raise ValueError(msg)


def _cast_repr_strings_with_schema(
    df: DataFrame, schema: dict[str, PolarsDataType | None]
) -> DataFrame:
    """
    Utility function to cast table repr/string values into frame-native types.

    Parameters
    ----------
    df
        Dataframe containing string-repr column data.
    schema
        DataFrame schema containing the desired end-state types.

    Notes
    -----
    Table repr strings are less strict (or different) than equivalent CSV data, so need
    special handling; as this function is only used for reprs, parsing is flexible.
    """
    tp: PolarsDataType | None
    if not df.is_empty():
        for tp in df.schema.values():
            if tp != String:
                msg = f"DataFrame should contain only String repr data; found {tp!r}"
                raise TypeError(msg)

    special_floats = {"-inf", "+inf", "inf", "nan"}

    # duration string scaling
    ns_sec = 1_000_000_000
    duration_scaling = {
        "ns": 1,
        "us": 1_000,
        "µs": 1_000,
        "ms": 1_000_000,
        "s": ns_sec,
        "m": ns_sec * 60,
        "h": ns_sec * 60 * 60,
        "d": ns_sec * 3_600 * 24,
        "w": ns_sec * 3_600 * 24 * 7,
    }

    # identify duration units and convert to nanoseconds
    def str_duration_(td: str | None) -> int | None:
        return (
            None
            if td is None
            else sum(
                int(value) * duration_scaling[unit.strip()]
                for value, unit in re.findall(r"([+-]?\d+)(\D+)", td)
            )
        )

    cast_cols = {}
    for c, tp in schema.items():
        if tp is not None:
            if tp.base_type() == Datetime:
                tp_base = Datetime(tp.time_unit)  # type: ignore[union-attr]
                d = F.col(c).str.replace(r"[A-Z ]+$", "")
                cast_cols[c] = (
                    F.when(d.str.len_bytes() == 19)
                    .then(d + ".000000000")
                    .otherwise(d + "000000000")
                    .str.slice(0, 29)
                    .str.strptime(tp_base, "%Y-%m-%d %H:%M:%S.%9f")
                )
                if getattr(tp, "time_zone", None) is not None:
                    cast_cols[c] = cast_cols[c].dt.replace_time_zone(tp.time_zone)  # type: ignore[union-attr]
            elif tp == Date:
                cast_cols[c] = F.col(c).str.strptime(tp, "%Y-%m-%d")  # type: ignore[arg-type]
            elif tp == Time:
                cast_cols[c] = (
                    F.when(F.col(c).str.len_bytes() == 8)
                    .then(F.col(c) + ".000000000")
                    .otherwise(F.col(c) + "000000000")
                    .str.slice(0, 18)
                    .str.strptime(tp, "%H:%M:%S.%9f")  # type: ignore[arg-type]
                )
            elif tp == Duration:
                cast_cols[c] = (
                    F.col(c)
                    .map_elements(str_duration_, return_dtype=Int64)
                    .cast(Duration("ns"))
                    .cast(tp)
                )
            elif tp == Boolean:
                cast_cols[c] = F.col(c).replace_strict({"true": True, "false": False})
            elif tp in INTEGER_DTYPES:
                int_string = F.col(c).str.replace_all(r"[^\d+-]", "")
                cast_cols[c] = (
                    pl.when(int_string.str.len_bytes() > 0).then(int_string).cast(tp)
                )
            elif tp in FLOAT_DTYPES or tp.base_type() == Decimal:
                # identify integer/fractional parts
                integer_part = F.col(c).str.replace(r"^(.*)\D(\d*)$", "$1")
                fractional_part = F.col(c).str.replace(r"^(.*)\D(\d*)$", "$2")
                cast_cols[c] = (
                    # check for empty string, special floats, or integer format
                    pl.when(
                        F.col(c).str.contains(r"^[+-]?\d*$")
                        | F.col(c).str.to_lowercase().is_in(special_floats)
                    )
                    .then(pl.when(F.col(c).str.len_bytes() > 0).then(F.col(c)))
                    # check for scientific notation
                    .when(F.col(c).str.contains("[eE]"))
                    .then(F.col(c).str.replace(r"[^eE\d]", "."))
                    .otherwise(
                        # recombine sanitised integer/fractional components
                        pl.concat_str(
                            integer_part.str.replace_all(r"[^\d+-]", ""),
                            fractional_part,
                            separator=".",
                        )
                    )
                    .cast(String)
                    .cast(tp)
                )
            elif tp != df.schema[c]:
                cast_cols[c] = F.col(c).cast(tp)

    return df.with_columns(**cast_cols) if cast_cols else df


# when building docs (with Sphinx) we need access to the functions
# associated with the namespaces from the class, as we don't have
# an instance; @sphinx_accessor is a @property that allows this.
NS = TypeVar("NS")


class sphinx_accessor(property):
    def __get__(  # type: ignore[override]
        self,
        instance: Any,
        cls: type[NS],
    ) -> NS:
        try:
            return self.fget(  # type: ignore[misc]
                instance if isinstance(instance, cls) else cls
            )
        except (AttributeError, ImportError):
            return self  # type: ignore[return-value]


BUILDING_SPHINX_DOCS = os.getenv("BUILDING_SPHINX_DOCS")


class _NoDefault(Enum):
    # "borrowed" from
    # https://github.com/pandas-dev/pandas/blob/e7859983a814b1823cf26e3b491ae2fa3be47c53/pandas/_libs/lib.pyx#L2736-L2748
    no_default = "NO_DEFAULT"

    def __repr__(self) -> str:
        return "<no_default>"


# the "no_default" sentinel should typically be used when one of the valid parameter
# values is None, as otherwise we cannot determine if the caller has set that value.
no_default = _NoDefault.no_default
NoDefault = Literal[_NoDefault.no_default]


def find_stacklevel() -> int:
    """
    Find the first place in the stack that is not inside Polars.

    Taken from:
    https://github.com/pandas-dev/pandas/blob/ab89c53f48df67709a533b6a95ce3d911871a0a8/pandas/util/_exceptions.py#L30-L51
    """
    pkg_dir = str(Path(pl.__file__).parent)

    # https://stackoverflow.com/questions/17407119/python-inspect-stack-is-slow
    frame = inspect.currentframe()
    n = 0
    try:
        while frame:
            fname = inspect.getfile(frame)
            if fname.startswith(pkg_dir) or (
                (qualname := getattr(frame.f_code, "co_qualname", None))
                # ignore @singledispatch wrappers
                and qualname.startswith("singledispatch.")
            ):
                frame = frame.f_back
                n += 1
            else:
                break
    finally:
        # https://docs.python.org/3/library/inspect.html
        # > Though the cycle detector will catch these, destruction of the frames
        # > (and local variables) can be made deterministic by removing the cycle
        # > in a 'finally' clause.
        del frame
    return n


def issue_warning(message: str, category: type[Warning], **kwargs: Any) -> None:
    """
    Issue a warning.

    Parameters
    ----------
    message
        The message associated with the warning.
    category
        The warning category.
    **kwargs
        Additional arguments for `warnings.warn`. Note that the `stacklevel` is
        determined automatically.
    """
    warnings.warn(
        message=message, category=category, stacklevel=find_stacklevel(), **kwargs
    )


def _get_stack_locals(
    of_type: type | Collection[type] | Callable[[Any], bool] | None = None,
    *,
    named: str | Collection[str] | None = None,
    n_objects: int | None = None,
    n_frames: int | None = None,
) -> dict[str, Any]:
    """
    Retrieve f_locals from all (or the last 'n') stack frames from the calling location.

    Parameters
    ----------
    of_type
        Only return objects of this type; can be a single class, tuple of
        classes, or a callable that returns True/False if the object being
        tested is considered a match.
    n_objects
        If specified, return only the most recent `n` matching objects.
    n_frames
        If specified, look at objects in the last `n` stack frames only.
    named
        If specified, only return objects matching the given name(s).
    """
    objects = {}
    examined_frames = 0

    if isinstance(named, str):
        named = (named,)
    if n_frames is None:
        n_frames = sys.maxsize

    if inspect.isfunction(of_type):
        matches_type = of_type
    else:
        if isinstance(of_type, Collection):
            of_type = tuple(of_type)

        def matches_type(obj: Any) -> bool:  # type: ignore[misc]
            return isinstance(obj, of_type)  # type: ignore[arg-type]

    if named is not None:
        if isinstance(named, str):
            named = (named,)
        elif not isinstance(named, set):
            named = set(named)

    stack_frame = inspect.currentframe()
    stack_frame = getattr(stack_frame, "f_back", None)
    try:
        while stack_frame and examined_frames < n_frames:
            local_items = list(stack_frame.f_locals.items())
            for nm, obj in reversed(local_items):
                if (
                    nm not in objects
                    and (named is None or nm in named)
                    and (of_type is None or matches_type(obj))
                ):
                    objects[nm] = obj
                    if n_objects is not None and len(objects) >= n_objects:
                        return objects

            stack_frame = stack_frame.f_back
            examined_frames += 1
    finally:
        # https://docs.python.org/3/library/inspect.html
        # > Though the cycle detector will catch these, destruction of the frames
        # > (and local variables) can be made deterministic by removing the cycle
        # > in a finally clause.
        del stack_frame

    return objects


# this is called from rust
def _polars_warn(msg: str, category: type[Warning] = UserWarning) -> None:
    warnings.warn(
        msg,
        category=category,
        stacklevel=find_stacklevel(),
    )


def extend_bool(
    value: bool | Sequence[bool],
    n_match: int,
    value_name: str,
    match_name: str,
) -> Sequence[bool]:
    """Ensure the given bool or sequence of bools is the correct length."""
    values = [value] * n_match if isinstance(value, bool) else value
    if n_match != len(values):
        msg = (
            f"the length of `{value_name}` ({len(values)}) "
            f"does not match the length of `{match_name}` ({n_match})"
        )
        raise ValueError(msg)
    return values


def in_terminal_that_supports_colour() -> bool:
    """
    Determine (within reason) if we are in an interactive terminal that supports color.

    Note: this is not exhaustive, but it covers a lot (most?) of the common cases.
    """
    if hasattr(sys.stdout, "isatty"):
        # can enhance as necessary, but this is a reasonable start
        return (
            sys.stdout.isatty()
            and (
                sys.platform != "win32"
                or "ANSICON" in os.environ
                or "WT_SESSION" in os.environ
                or os.environ.get("TERM_PROGRAM") == "vscode"
                or os.environ.get("TERM") == "xterm-256color"
            )
        ) or os.environ.get("PYCHARM_HOSTED") == "1"
    return False


def parse_percentiles(
    percentiles: Sequence[float] | float | None, *, inject_median: bool = False
) -> Sequence[float]:
    """
    Transforms raw percentiles into our preferred format, adding the 50th percentile.

    Raises a ValueError if the percentile sequence is invalid
    (e.g. outside the range [0, 1])
    """
    if isinstance(percentiles, float):
        percentiles = [percentiles]
    elif percentiles is None:
        percentiles = []
    if not all((0 <= p <= 1) for p in percentiles):
        msg = "`percentiles` must all be in the range [0, 1]"
        raise ValueError(msg)

    sub_50_percentiles = sorted(p for p in percentiles if p < 0.5)
    at_or_above_50_percentiles = sorted(p for p in percentiles if p >= 0.5)

    if inject_median and (
        not at_or_above_50_percentiles or at_or_above_50_percentiles[0] != 0.5
    ):
        at_or_above_50_percentiles = [0.5, *at_or_above_50_percentiles]

    return [*sub_50_percentiles, *at_or_above_50_percentiles]


def re_escape(s: str) -> str:
    """Escape a string for use in a Polars (Rust) regex."""
    # note: almost the same as the standard python 're.escape' function, but
    # escapes _only_ those metachars with meaning to the rust regex crate
    re_rust_metachars = r"\\?()|\[\]{}^$#&~.+*-"
    return re.sub(f"([{re_rust_metachars}])", r"\\\1", s)


# Don't rename or move. This is used by polars cloud
def display_dot_graph(
    *,
    dot: str,
    show: bool = True,
    output_path: str | Path | None = None,
    raw_output: bool = False,
    figsize: tuple[float, float] = (16.0, 12.0),
) -> str | None:
    if raw_output:
        # we do not show a graph, nor save a graph to disk
        return dot

    output_type = (
        "svg"
        if _in_notebook()
        or _in_marimo_notebook()
        or "POLARS_DOT_SVG_VIEWER" in os.environ
        else "png"
    )

    try:
        graph = subprocess.check_output(
            ["dot", "-Nshape=box", "-T" + output_type], input=f"{dot}".encode()
        )
    except (ImportError, FileNotFoundError):
        msg = (
            "the graphviz `dot` binary should be on your PATH."
            "(If not installed you can download here: https://graphviz.org/download/)"
        )
        raise ImportError(msg) from None

    if output_path:
        Path(output_path).write_bytes(graph)

    if not show:
        return None

    if _in_notebook():
        from IPython.display import SVG, display

        return display(SVG(graph))
    elif _in_marimo_notebook():
        import marimo as mo

        return mo.Html(f"{graph.decode()}")
    else:
        if (cmd := os.environ.get("POLARS_DOT_SVG_VIEWER", None)) is not None:
            import tempfile

            with tempfile.NamedTemporaryFile(suffix=".svg") as file:
                file.write(graph)
                file.flush()
                cmd = cmd.replace("%file%", file.name)
                subprocess.run(cmd, shell=True)
            return None

        import_optional(
            "matplotlib",
            err_prefix="",
            err_suffix="should be installed to show graphs",
        )
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt

        plt.figure(figsize=figsize)
        img = mpimg.imread(BytesIO(graph))
        plt.axis("off")
        plt.imshow(img)
        plt.show()
        return None


def qualified_type_name(obj: Any, *, qualify_polars: bool = False) -> str:
    """
    Return the module-qualified name of the given object as a string.

    Parameters
    ----------
    obj
        The object to get the qualified name for.
    qualify_polars
        If False (default), omit the module path for our own (Polars) objects.
    """
    if isinstance(obj, type):
        module = obj.__module__
        name = obj.__name__
    else:
        module = obj.__class__.__module__
        name = obj.__class__.__name__

    if (
        not module
        or module == "builtins"
        or (not qualify_polars and module.startswith("polars."))
    ):
        return name

    return f"{module}.{name}"


def require_same_type(current: Any, other: Any) -> None:
    """
    Raise an error if the two arguments are not of the same type.

    The check will not raise an error if one object is of a subclass of the other.

    Parameters
    ----------
    current
        The object the type of which is being checked against.
    other
        An object that has to be of the same type.
    """
    if not isinstance(other, type(current)) and not isinstance(current, type(other)):
        msg = (
            f"expected `other` to be a {qualified_type_name(current)!r}, "
            f"not {qualified_type_name(other)!r}"
        )
        raise TypeError(msg)
