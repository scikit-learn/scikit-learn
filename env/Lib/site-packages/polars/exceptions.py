try:
    from polars.polars import (
        CategoricalRemappingWarning,
        ColumnNotFoundError,
        ComputeError,
        DuplicateError,
        InvalidOperationError,
        MapWithoutReturnDtypeWarning,
        NoDataError,
        OutOfBoundsError,
        PanicException,
        PerformanceWarning,
        PolarsError,
        PolarsWarning,
        SchemaError,
        SchemaFieldNotFoundError,
        ShapeError,
        SQLInterfaceError,
        SQLSyntaxError,
        StringCacheMismatchError,
        StructFieldNotFoundError,
    )
except ImportError:
    # redefined for documentation purposes when there is no binary

    class PolarsError(Exception):  # type: ignore[no-redef]
        """Base class for all Polars errors."""

    class ColumnNotFoundError(PolarsError):  # type: ignore[no-redef, misc]
        """
        Exception raised when a specified column is not found.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 2, 3]})
        >>> df.select("b")
        polars.exceptions.ColumnNotFoundError: b
        """

    class ComputeError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when Polars could not perform an underlying computation."""

    class DuplicateError(PolarsError):  # type: ignore[no-redef, misc]
        """
        Exception raised when a column name is duplicated.

        Examples
        --------
        >>> df = pl.DataFrame({"a": [1, 1, 1]})
        >>> pl.concat([df, df], how="horizontal")
        polars.exceptions.DuplicateError: unable to hstack, column with name "a" already exists
        """  # noqa: W505

    class InvalidOperationError(PolarsError):  # type: ignore[no-redef, misc]
        """
        Exception raised when an operation is not allowed (or possible) against a given object or data structure.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.is_in(["x", "y"])
        polars.exceptions.InvalidOperationError: `is_in` cannot check for String values in Int64 data
        """  # noqa: W505

    class NoDataError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when an operation cannot be performed on an empty data structure."""  # noqa: W505

    class OutOfBoundsError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when the given index is out of bounds."""

    class PanicException(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when an unexpected state causes a panic in the underlying Rust library."""  # noqa: W505

    class SchemaError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when an unexpected schema mismatch causes an error."""

    class SchemaFieldNotFoundError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when a specified schema field is not found."""

    class ShapeError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when trying to perform operations on data structures with incompatible shapes."""  # noqa: W505

    class SQLInterfaceError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when an error occurs in the SQL interface."""

    class SQLSyntaxError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised from the SQL interface when encountering invalid syntax."""

    class StringCacheMismatchError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when string caches come from different sources."""

    class StructFieldNotFoundError(PolarsError):  # type: ignore[no-redef, misc]
        """Exception raised when a specified Struct field is not found."""

    class PolarsWarning(Exception):  # type: ignore[no-redef]
        """Base class for all Polars warnings."""

    class PerformanceWarning(PolarsWarning):  # type: ignore[no-redef, misc]
        """Warning issued to indicate potential performance pitfalls."""

    class CategoricalRemappingWarning(PerformanceWarning):  # type: ignore[no-redef, misc]
        """Warning issued when a categorical needs to be remapped to be compatible with another categorical."""  # noqa: W505

    class MapWithoutReturnDtypeWarning(PolarsWarning):  # type: ignore[no-redef, misc]
        """Warning issued when `map_elements` is performed without specifying the return dtype."""  # noqa: W505


class RowsError(PolarsError):  # type: ignore[misc]
    """Exception raised when the number of returned rows does not match expectation."""


class NoRowsReturnedError(RowsError):
    """Exception raised when no rows are returned, but at least one row is expected."""


class TooManyRowsReturnedError(RowsError):
    """Exception raised when more rows than expected are returned."""


class ModuleUpgradeRequiredError(ModuleNotFoundError):
    """Exception raised when a module is installed but needs to be upgraded."""


class ParameterCollisionError(PolarsError):  # type: ignore[misc]
    """Exception raised when the same parameter occurs multiple times."""


class UnsuitableSQLError(PolarsError):  # type: ignore[misc]
    """Exception raised when unsuitable SQL is given to a database method."""


class ChronoFormatWarning(PolarsWarning):  # type: ignore[misc]
    """
    Warning issued when a chrono format string contains dubious patterns.

    Polars uses Rust's chrono crate to convert between string data and temporal data.
    The patterns used by chrono differ slightly from Python's built-in datetime module.
    Refer to the `chrono strftime documentation
    <https://docs.rs/chrono/latest/chrono/format/strftime/index.html>`_ for the full
    specification.
    """


class CustomUFuncWarning(PolarsWarning):  # type: ignore[misc]
    """Warning issued when a custom ufunc is handled differently than numpy ufunc would."""  # noqa: W505


class DataOrientationWarning(PolarsWarning):  # type: ignore[misc]
    """
    Warning issued to indicate row orientation was inferred from the inputs.

    Occurs when constructing a DataFrame from a list of rows without explicitly
    specifying row orientation. Polars is usually able to infer the data orientation
    from the data and schema, but there are cases where this is not possible. This is a
    common source of confusion. Use the `orient` parameter to be explicit about the
    data orientation.

    Examples
    --------
    >>> pl.DataFrame([(1, 2, 3), (4, 5, 6)], schema=["a", "b", "c"])  # doctest: +SKIP
    DataOrientationWarning: Row orientation inferred during DataFrame construction.
    Explicitly specify the orientation by passing `orient="row"` to silence this warning.
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ a   ┆ b   ┆ c   │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 2   ┆ 3   │
    │ 4   ┆ 5   ┆ 6   │
    └─────┴─────┴─────┘

    Pass `orient="row"` to silence the warning.

    >>> pl.DataFrame([[1, 2, 3], [4, 5, 6]], schema=["a", "b", "c"], orient="row")
    shape: (2, 3)
    ┌─────┬─────┬─────┐
    │ a   ┆ b   ┆ c   │
    │ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╡
    │ 1   ┆ 2   ┆ 3   │
    │ 4   ┆ 5   ┆ 6   │
    └─────┴─────┴─────┘
    """  # noqa: W505


class PolarsInefficientMapWarning(PerformanceWarning):  # type: ignore[misc]
    """Warning issued when a potentially slow `map_*` operation is performed."""


class UnstableWarning(PolarsWarning):  # type: ignore[misc]
    """Warning issued when unstable functionality is used."""


__all__ = [
    # Errors
    "PolarsError",
    "ColumnNotFoundError",
    "ComputeError",
    "DuplicateError",
    "InvalidOperationError",
    "ModuleUpgradeRequiredError",
    "NoDataError",
    "NoRowsReturnedError",
    "OutOfBoundsError",
    "ParameterCollisionError",
    "RowsError",
    "SQLInterfaceError",
    "SQLSyntaxError",
    "SchemaError",
    "SchemaFieldNotFoundError",
    "ShapeError",
    "StringCacheMismatchError",
    "StructFieldNotFoundError",
    "TooManyRowsReturnedError",
    "UnsuitableSQLError",
    # Warnings
    "PolarsWarning",
    "CategoricalRemappingWarning",
    "ChronoFormatWarning",
    "CustomUFuncWarning",
    "DataOrientationWarning",
    "MapWithoutReturnDtypeWarning",
    "PerformanceWarning",
    "PolarsInefficientMapWarning",
    "UnstableWarning",
    # Panic
    "PanicException",
]
