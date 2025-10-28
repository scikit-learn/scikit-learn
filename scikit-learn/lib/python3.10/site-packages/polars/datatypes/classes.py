from __future__ import annotations

import contextlib
import enum
from collections import OrderedDict
from collections.abc import Mapping
from datetime import tzinfo
from inspect import isclass
from typing import TYPE_CHECKING, Any, Callable, Generic, TypeVar, overload

import polars._reexport as pl
import polars.datatypes
import polars.functions as F

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr
    from polars._plr import PyCategories
    from polars._plr import dtype_str_repr as _dtype_str_repr

import polars.datatypes.classes as pldt

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

    from polars import Series
    from polars._typing import (
        CategoricalOrdering,
        PolarsDataType,
        PythonDataType,
        SchemaDict,
        TimeUnit,
    )


T = TypeVar("T")
R = TypeVar("R")


class classinstmethod(Generic[R]):
    """Decorator that allows a method to be called from the class OR instance."""

    def __init__(self, func: Callable[..., R]) -> None:
        self.func = func

    def __get__(self, instance: Any, type_: Any) -> Callable[..., R]:
        if instance is not None:
            return self.func.__get__(instance, type_)
        return self.func.__get__(type_, type_)


class DataTypeClass(type):
    """Metaclass for nicely printing DataType classes."""

    def __repr__(cls) -> str:
        return cls.__name__

    def _string_repr(cls) -> str:
        return _dtype_str_repr(cls)

    # Methods below defined here in signature only to satisfy mypy

    @classmethod
    def base_type(cls) -> DataTypeClass:  # noqa: D102
        ...

    @classmethod
    def is_(cls, other: PolarsDataType) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_numeric(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_decimal(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_integer(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_object(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_signed_integer(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_unsigned_integer(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_float(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_temporal(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def is_nested(cls) -> bool:  # noqa: D102
        ...

    @classmethod
    def from_python(cls, py_type: PythonDataType) -> PolarsDataType:  # noqa: D102
        ...

    @classmethod
    def to_python(cls) -> PythonDataType:  # noqa: D102
        ...

    @classmethod
    def to_dtype_expr(cls) -> pl.DataTypeExpr:  # noqa: D102
        ...


class DataType(metaclass=DataTypeClass):
    """Base class for all Polars data types."""

    def _string_repr(self) -> str:
        return _dtype_str_repr(self)

    @overload  # type: ignore[override]
    def __eq__(self, other: pl.DataTypeExpr) -> pl.Expr: ...

    @overload
    def __eq__(self, other: PolarsDataType) -> bool: ...

    def __eq__(self, other: pl.DataTypeExpr | PolarsDataType) -> pl.Expr | bool:
        if isinstance(other, pl.DataTypeExpr):
            return self.to_dtype_expr() == other
        elif type(other) is DataTypeClass:
            return issubclass(other, type(self))
        else:
            return isinstance(other, type(self))

    def __hash__(self) -> int:
        return hash(self.__class__)

    def __repr__(self) -> str:
        return self.__class__.__name__

    @classmethod
    def base_type(cls) -> DataTypeClass:
        """
        Return this DataType's fundamental/root type class.

        Examples
        --------
        >>> pl.Datetime("ns").base_type()
        Datetime
        >>> pl.List(pl.Int32).base_type()
        List
        >>> pl.Struct([pl.Field("a", pl.Int64), pl.Field("b", pl.Boolean)]).base_type()
        Struct
        """
        return cls

    @classinstmethod
    def is_(self, other: PolarsDataType) -> bool:
        """
        Check if this DataType is the same as another DataType.

        This is a stricter check than `self == other`, as it enforces an exact
        match of all dtype attributes for nested and/or uninitialised dtypes.

        Parameters
        ----------
        other
            the other Polars dtype to compare with.

        Examples
        --------
        >>> pl.List == pl.List(pl.Int32)
        True
        >>> pl.List.is_(pl.List(pl.Int32))
        False
        """
        return self == other and hash(self) == hash(other)

    @classmethod
    def is_numeric(cls) -> bool:
        """Check whether the data type is a numeric type."""
        return issubclass(cls, NumericType)

    @classmethod
    def is_decimal(cls) -> bool:
        """Check whether the data type is a decimal type."""
        return issubclass(cls, Decimal)

    @classmethod
    def is_integer(cls) -> bool:
        """Check whether the data type is an integer type."""
        return issubclass(cls, IntegerType)

    @classmethod
    def is_object(cls) -> bool:
        """Check whether the data type is an object type."""
        return issubclass(cls, ObjectType)

    @classmethod
    def is_signed_integer(cls) -> bool:
        """Check whether the data type is a signed integer type."""
        return issubclass(cls, SignedIntegerType)

    @classmethod
    def is_unsigned_integer(cls) -> bool:
        """Check whether the data type is an unsigned integer type."""
        return issubclass(cls, UnsignedIntegerType)

    @classmethod
    def is_float(cls) -> bool:
        """Check whether the data type is a floating point type."""
        return issubclass(cls, FloatType)

    @classmethod
    def is_temporal(cls) -> bool:
        """Check whether the data type is a temporal type."""
        return issubclass(cls, TemporalType)

    @classmethod
    def is_nested(cls) -> bool:
        """Check whether the data type is a nested type."""
        return issubclass(cls, NestedType)

    @classmethod
    def from_python(cls, py_type: PythonDataType) -> PolarsDataType:
        """
        Return the Polars data type corresponding to a given Python type.

        Notes
        -----
        Not every Python type has a corresponding Polars data type; in general
        you should declare Polars data types explicitly to exactly specify
        the desired type and its properties (such as scale/unit).

        Examples
        --------
        >>> pl.DataType.from_python(int)
        Int64
        >>> pl.DataType.from_python(float)
        Float64
        >>> from datetime import tzinfo
        >>> pl.DataType.from_python(tzinfo)  # doctest: +SKIP
        TypeError: cannot parse input <class 'datetime.tzinfo'> into Polars data type
        """
        from polars.datatypes._parse import parse_into_dtype

        return parse_into_dtype(py_type)

    @classinstmethod
    def to_python(self) -> PythonDataType:
        """
        Return the Python type corresponding to this Polars data type.

        Examples
        --------
        >>> pl.Int16().to_python()
        <class 'int'>
        >>> pl.Float32().to_python()
        <class 'float'>
        >>> pl.Array(pl.Date(), 10).to_python()
        <class 'list'>
        """
        from polars.datatypes import dtype_to_py_type

        return dtype_to_py_type(self)

    @classinstmethod
    def to_dtype_expr(self) -> pl.DataTypeExpr:
        """
        Return a :class:`DataTypeExpr` with a static :class:`DataType`.

        Examples
        --------
        >>> pl.Int16().to_dtype_expr().collect_dtype({})
        Int16
        """
        from polars._plr import PyDataTypeExpr

        return pl.DataTypeExpr._from_pydatatype_expr(PyDataTypeExpr.from_dtype(self))


class NumericType(DataType):
    """Base class for numeric data types."""

    @classmethod
    def max(cls) -> pl.Expr:
        """
        Return a literal expression representing the maximum value of this data type.

        Examples
        --------
        >>> pl.select(pl.Int8.max() == 127)
        shape: (1, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ bool    │
        ╞═════════╡
        │ true    │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(plr._get_dtype_max(cls))

    @classmethod
    def min(cls) -> pl.Expr:
        """
        Return a literal expression representing the minimum value of this data type.

        Examples
        --------
        >>> pl.select(pl.Int8.min() == -128)
        shape: (1, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ bool    │
        ╞═════════╡
        │ true    │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(plr._get_dtype_min(cls))


class IntegerType(NumericType):
    """Base class for integer data types."""


class SignedIntegerType(IntegerType):
    """Base class for signed integer data types."""


class UnsignedIntegerType(IntegerType):
    """Base class for unsigned integer data types."""


class FloatType(NumericType):
    """Base class for float data types."""


class TemporalType(DataType):
    """Base class for temporal data types."""


class NestedType(DataType):
    """Base class for nested data types."""


class ObjectType(DataType):
    """Base class for object data types."""


class Int8(SignedIntegerType):
    """8-bit signed integer type."""


class Int16(SignedIntegerType):
    """16-bit signed integer type."""


class Int32(SignedIntegerType):
    """32-bit signed integer type."""


class Int64(SignedIntegerType):
    """64-bit signed integer type."""


class Int128(SignedIntegerType):
    """
    128-bit signed integer type.

    .. warning::
        This functionality is considered **unstable**.
        It is a work-in-progress feature and may not always work as expected.
        It may be changed at any point without it being considered a breaking change.
    """


class UInt8(UnsignedIntegerType):
    """8-bit unsigned integer type."""


class UInt16(UnsignedIntegerType):
    """16-bit unsigned integer type."""


class UInt32(UnsignedIntegerType):
    """32-bit unsigned integer type."""


class UInt64(UnsignedIntegerType):
    """64-bit unsigned integer type."""


class UInt128(UnsignedIntegerType):
    """128-bit unsigned integer type.

    .. warning::
        This functionality is considered **unstable**.
        It is a work-in-progress feature and may not always work as expected.
        It may be changed at any point without it being considered a breaking change.
    """


class Float32(FloatType):
    """32-bit floating point type."""


class Float64(FloatType):
    """64-bit floating point type."""


class Decimal(NumericType):
    """
    Decimal 128-bit type with an optional precision and non-negative scale.

    Parameters
    ----------
    precision
        Maximum number of digits in each number.
        If set to `None` (default), the precision is set to 38 (the maximum
        supported by Polars).
    scale
        Number of digits to the right of the decimal point in each number.
    """

    precision: int | None
    scale: int

    def __init__(
        self,
        precision: int | None = None,
        scale: int = 0,
    ) -> None:
        if precision is None:
            precision = 38

        self.precision = precision
        self.scale = scale

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(precision={self.precision}, scale={self.scale})"
        )

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, Decimal):
            return True
        elif isinstance(other, Decimal):
            return self.precision == other.precision and self.scale == other.scale
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, self.precision, self.scale))


class Boolean(DataType):
    """Boolean type."""


class String(DataType):
    """UTF-8 encoded string type."""


# Allow Utf8 as an alias for String
Utf8 = String


class Binary(DataType):
    """Binary type."""


class Date(TemporalType):
    """
    Data type representing a calendar date.

    Notes
    -----
    The underlying representation of this type is a 32-bit signed integer.
    The integer indicates the number of days since the Unix epoch (1970-01-01).
    The number can be negative to indicate dates before the epoch.
    """


class Time(TemporalType):
    """
    Data type representing the time of day.

    Notes
    -----
    The underlying representation of this type is a 64-bit signed integer.
    The integer indicates the number of nanoseconds since midnight.
    """

    @classmethod
    def max(cls) -> pl.Expr:
        """
        Return a literal expression representing the maximum value of this data type.

        Examples
        --------
        >>> pl.select(pl.Time.max() == 86_399_999_999_999)
        shape: (1, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ bool    │
        ╞═════════╡
        │ true    │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(plr._get_dtype_max(cls))

    @classmethod
    def min(cls) -> pl.Expr:
        """
        Return a literal expression representing the minimum value of this data type.

        Examples
        --------
        >>> pl.select(pl.Time.min() == 0)
        shape: (1, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ bool    │
        ╞═════════╡
        │ true    │
        └─────────┘
        """
        return pl.Expr._from_pyexpr(plr._get_dtype_min(cls))


class Datetime(TemporalType):
    """
    Data type representing a calendar date and time of day.

    Parameters
    ----------
    time_unit : {'us', 'ns', 'ms'}
        Unit of time. Defaults to `'us'` (microseconds).
    time_zone
        Time zone string, as defined in zoneinfo (to see valid strings run
        `import zoneinfo; zoneinfo.available_timezones()` for a full list).
        When used to match dtypes, can set this to "*" to check for Datetime
        columns that have any (non-null) timezone.

    Notes
    -----
    The underlying representation of this type is a 64-bit signed integer.
    The integer indicates the number of time units since the Unix epoch
    (1970-01-01 00:00:00). The number can be negative to indicate datetimes before the
    epoch.
    """

    time_unit: TimeUnit
    time_zone: str | None

    def __init__(
        self, time_unit: TimeUnit = "us", time_zone: str | tzinfo | None = None
    ) -> None:
        if time_unit not in ("ms", "us", "ns"):
            msg = (
                "invalid `time_unit`"
                f"\n\nExpected one of {{'ns','us','ms'}}, got {time_unit!r}."
            )
            raise ValueError(msg)

        if isinstance(time_zone, tzinfo):
            time_zone = str(time_zone)

        self.time_unit = time_unit
        self.time_zone = time_zone

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, Datetime):
            return True
        elif isinstance(other, Datetime):
            return (
                self.time_unit == other.time_unit and self.time_zone == other.time_zone
            )
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, self.time_unit, self.time_zone))

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return (
            f"{class_name}(time_unit={self.time_unit!r}, time_zone={self.time_zone!r})"
        )


class Duration(TemporalType):
    """
    Data type representing a time duration.

    Parameters
    ----------
    time_unit : {'us', 'ns', 'ms'}
        Unit of time. Defaults to `'us'` (microseconds).

    Notes
    -----
    The underlying representation of this type is a 64-bit signed integer.
    The integer indicates an amount of time units and can be negative to indicate
    negative time offsets.
    """

    time_unit: TimeUnit

    def __init__(self, time_unit: TimeUnit = "us") -> None:
        if time_unit not in ("ms", "us", "ns"):
            msg = (
                "invalid `time_unit`"
                f"\n\nExpected one of {{'ns','us','ms'}}, got {time_unit!r}."
            )
            raise ValueError(msg)

        self.time_unit = time_unit

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, Duration):
            return True
        elif isinstance(other, Duration):
            return self.time_unit == other.time_unit
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, self.time_unit))

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}(time_unit={self.time_unit!r})"


class Categories:
    """
    A named collection of categories for `Categorical`.

    Two categories are considered equal (and will use the same physical mapping of
    categories to strings) if they have the same name, namespace and physical backing
    type, even if they are created in separate calls to `Categories`.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.
    """

    _categories: PyCategories

    def __init__(
        self,
        name: str | None = None,
        namespace: str = "",
        physical: PolarsDataType = pldt.UInt32,
    ) -> None:
        if name is None or name == "":
            assert namespace == "", "global categories may not specify a namespace"
            assert physical == pldt.UInt32, (
                "global categories may not specify a physical type"
            )
            self._categories = PyCategories.global_categories()
            return

        if physical == pldt.UInt32:
            internal_phys = "u32"
        elif physical == pldt.UInt16:
            internal_phys = "u16"
        elif physical == pldt.UInt8:
            internal_phys = "u8"
        else:
            msg = "Categorical physical must be one of pl.UInt(8|16|32)"
            raise TypeError(msg)

        self._categories = PyCategories(name, namespace, internal_phys)

    @staticmethod
    def _from_py_categories(py_categories: PyCategories) -> Categories:
        self = Categories.__new__(Categories)
        self._categories = py_categories
        return self

    @staticmethod
    def random(
        namespace: str = "", physical: PolarsDataType = pldt.UInt32
    ) -> Categories:
        """Creates a new Categories with a random name."""
        if physical == pldt.UInt32:
            internal_phys = "u32"
        elif physical == pldt.UInt16:
            internal_phys = "u16"
        elif physical == pldt.UInt8:
            internal_phys = "u8"
        else:
            msg = "Categorical physical must be one of pl.UInt(8|16|32)"
            raise TypeError(msg)

        return Categories._from_py_categories(
            PyCategories.random(namespace, internal_phys)
        )

    def name(self) -> str:
        """The name of this `Categories`."""
        return self._categories.name()

    def namespace(self) -> str:
        """The namespace of this `Categories`."""
        return self._categories.namespace()

    def physical(self) -> PolarsDataType:
        """The physical type used to represent the categories."""
        phys = self._categories.physical()
        if phys == "u8":
            return pldt.UInt8
        elif phys == "u16":
            return pldt.UInt16
        elif phys == "u32":
            return pldt.UInt32
        else:
            msg = "unknown physical dtype"
            raise RuntimeError(msg)

    def is_global(self) -> bool:
        """Returns whether this refers to the global categories."""
        return self._categories.is_global()

    def __getitem__(self, key: str | int | None) -> str | int | None:
        if key is None:
            return key
        elif isinstance(key, str):
            return self._categories.get_cat(key)
        else:
            return self._categories.cat_to_str(key)

    def __repr__(self) -> str:
        name = self.name()
        namespace = self.namespace()
        phys = self.physical()
        if self._categories.is_global():
            return "Categories()"
        elif namespace == "" and phys == pldt.UInt32:
            return f'Categories("{name}")'
        else:
            return f'Categories(name="{name}", namespace="{namespace}", physical=pl.{phys})'

    def __hash__(self) -> int:
        return hash(self._categories)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Categories) and self._categories == other._categories

    def __getstate__(self) -> tuple[str, str, PolarsDataType]:
        return self.name(), self.namespace(), self.physical()

    def __setstate__(self, state: tuple[str, str, PolarsDataType]) -> None:
        self.__dict__ = Categories(*state).__dict__


class Categorical(DataType):
    """
    A categorical encoding of a set of strings.

    Parameters
    ----------
    ordering : {'lexical', 'physical'}
        Ordering by order of appearance (`'physical'`, default)
        or string value (`'lexical'`).

        .. deprecated:: 1.32.0
            Parameter is now ignored. Always behaves as if `'lexical'` was passed.
    """

    ordering: CategoricalOrdering | None
    categories: Categories

    def __init__(
        self,
        ordering: CategoricalOrdering | Categories | None = "lexical",
        **kwargs: Any,
    ) -> None:
        # Future API will be this, already support it for backwards compat.
        if isinstance(ordering, Categories):
            self.ordering = "lexical"
            self.categories = ordering
            assert len(kwargs) == 0
            return

        if ordering == "physical":
            from polars._utils.deprecation import issue_deprecation_warning

            issue_deprecation_warning(
                "the physical Categorical ordering is deprecated. The ordering is now always lexical.",
                version="1.32.0",
            )

        self.ordering = "lexical"
        if kwargs.get("categories") is not None:
            assert len(kwargs) == 1
            self.categories = kwargs["categories"]
        else:
            assert len(kwargs) == 0
            self.categories = Categories()

    def __repr__(self) -> str:
        if self.categories.is_global():
            return f"{self.__class__.__name__}"
        else:
            return f"{self.__class__.__name__}({self.categories!r})"

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, Categorical):
            return self.categories.is_global()
        elif isinstance(other, Categorical):
            return self.categories == other.categories
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, self.categories))


class Enum(DataType):
    """
    A fixed categorical encoding of a unique set of strings.

    Parameters
    ----------
    categories
        The categories in the dataset; must be a unique set of strings, or an
        existing Python string-valued enum.

    Examples
    --------
    Explicitly define enumeration categories:

    >>> pl.Enum(["north", "south", "east", "west"])
    Enum(categories=['north', 'south', 'east', 'west'])

    Initialise from an existing Python enumeration:

    >>> from http import HTTPMethod
    >>> pl.Enum(HTTPMethod)
    Enum(categories=['CONNECT', 'DELETE', 'GET', 'HEAD', 'OPTIONS', 'PATCH', 'POST', 'PUT', 'TRACE'])
    """  # noqa: W505

    categories: Series

    def __init__(self, categories: Series | Iterable[str] | type[enum.Enum]) -> None:
        if isclass(categories) and issubclass(categories, enum.Enum):
            for enum_subclass in (enum.Flag, enum.IntEnum):
                if issubclass(categories, enum_subclass):
                    enum_type_name = categories.__name__
                    msg = f"Enum categories must be strings; `{enum_type_name}` values are integers"
                    raise TypeError(msg)

            enum_values = [
                getattr(v, "value", v) for v in categories.__members__.values()
            ]
            categories = pl.Series(values=enum_values)
        elif not isinstance(categories, pl.Series):
            categories = pl.Series(values=categories)

        if categories.is_empty():
            self.categories = pl.Series(name="category", dtype=String)
            return

        if categories.has_nulls():
            msg = "Enum categories must not contain null values"
            raise TypeError(msg)

        if (dtype := categories.dtype) != String:
            msg = f"Enum categories must be strings; found data of type {dtype}"
            raise TypeError(msg)

        if categories.n_unique() != categories.len():
            duplicate = categories.filter(categories.is_duplicated())[0]
            msg = f"Enum categories must be unique; found duplicate {duplicate!r}"
            raise ValueError(msg)

        self.categories = categories.rechunk().alias("category")

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, Enum):
            return True
        elif isinstance(other, Enum):
            return self.categories.equals(other.categories)
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, tuple(self.categories)))

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}(categories={self.categories.to_list()!r})"

    def union(self, other: Enum) -> Enum:
        """Union of two Enums."""
        return Enum(
            F.concat((self.categories, other.categories)).unique(maintain_order=True)
        )

    __or__ = union


class Object(ObjectType):
    """Data type for wrapping arbitrary Python objects."""


class Null(DataType):
    """Data type representing null values."""


class Unknown(DataType):
    """Type representing DataType values that could not be determined statically."""


class List(NestedType):
    """
    Variable length list type.

    Parameters
    ----------
    inner
        The `DataType` of the values within each list.

    Examples
    --------
    >>> df = pl.DataFrame(
    ...     {
    ...         "integer_lists": [[1, 2], [3, 4]],
    ...         "float_lists": [[1.0, 2.0], [3.0, 4.0]],
    ...     }
    ... )
    >>> df
    shape: (2, 2)
    ┌───────────────┬─────────────┐
    │ integer_lists ┆ float_lists │
    │ ---           ┆ ---         │
    │ list[i64]     ┆ list[f64]   │
    ╞═══════════════╪═════════════╡
    │ [1, 2]        ┆ [1.0, 2.0]  │
    │ [3, 4]        ┆ [3.0, 4.0]  │
    └───────────────┴─────────────┘
    """

    inner: PolarsDataType

    def __init__(self, inner: PolarsDataType | PythonDataType) -> None:
        self.inner = polars.datatypes.parse_into_dtype(inner)

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # This equality check allows comparison of type classes and type instances.
        # If a parent type is not specific about its inner type, we infer it as equal:
        # > list[i64] == list[i64] -> True
        # > list[i64] == list[f32] -> False
        # > list[i64] == list      -> True

        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, List):
            return True
        elif isinstance(other, List):
            return self.inner == other.inner
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, self.inner))

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}({self.inner!r})"


class Array(NestedType):
    """
    Fixed length list type.

    Parameters
    ----------
    inner
        The `DataType` of the values within each array.
    shape
        The shape of the arrays.
    width
        The length of the arrays.

        .. deprecated:: 0.20.31
            The `width` parameter for `Array` is deprecated. Use `shape` instead.

    Examples
    --------
    >>> s = pl.Series("a", [[1, 2], [4, 3]], dtype=pl.Array(pl.Int64, 2))
    >>> s
    shape: (2,)
    Series: 'a' [array[i64, 2]]
    [
            [1, 2]
            [4, 3]
    ]
    """

    inner: PolarsDataType
    size: int
    shape: tuple[int, ...]

    def __init__(
        self,
        inner: PolarsDataType | PythonDataType,
        shape: int | tuple[int, ...] | None = None,
        *,
        width: int | None = None,
    ) -> None:
        if width is not None:
            from polars._utils.deprecation import issue_deprecation_warning

            issue_deprecation_warning(
                "the `width` parameter for `Array` is deprecated. Use `shape` instead.",
                version="0.20.31",
            )
            shape = width
        elif shape is None:
            msg = "Array constructor is missing the required argument `shape`"
            raise TypeError(msg)

        inner_parsed = polars.datatypes.parse_into_dtype(inner)
        inner_shape = inner_parsed.shape if isinstance(inner_parsed, Array) else ()

        if isinstance(shape, int):
            self.inner = inner_parsed
            self.size = shape
            self.shape = (shape,) + inner_shape

        elif isinstance(shape, tuple) and isinstance(shape[0], int):  # type: ignore[redundant-expr]
            if len(shape) > 1:
                inner_parsed = Array(inner_parsed, shape[1:])

            self.inner = inner_parsed
            self.size = shape[0]
            self.shape = shape + inner_shape

        else:
            msg = f"invalid input for shape: {shape!r}"
            raise TypeError(msg)

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # This equality check allows comparison of type classes and type instances.
        # If a parent type is not specific about its inner type, we infer it as equal:
        # > array[i64] == array[i64] -> True
        # > array[i64] == array[f32] -> False
        # > array[i64] == array      -> True

        # allow comparing object instances to class
        if type(other) is DataTypeClass and issubclass(other, Array):
            return True
        elif isinstance(other, Array):
            if self.shape != other.shape:
                return False
            else:
                return self.inner == other.inner
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, self.inner, self.size))

    def __repr__(self) -> str:
        # Get leaf type
        dtype = self.inner
        while isinstance(dtype, Array):
            dtype = dtype.inner

        class_name = self.__class__.__name__
        return f"{class_name}({dtype!r}, shape={self.shape})"

    @property
    def width(self) -> int:
        """The size of the Array."""
        from polars._utils.deprecation import issue_deprecation_warning

        issue_deprecation_warning(
            "the `width` attribute for `Array` is deprecated. Use `size` instead.",
            version="0.20.31",
        )
        return self.size


class Field:
    """
    Definition of a single field within a `Struct` DataType.

    Parameters
    ----------
    name
        The name of the field within its parent `Struct`.
    dtype
        The `DataType` of the field's values.
    """

    name: str
    dtype: PolarsDataType

    def __init__(self, name: str, dtype: PolarsDataType) -> None:
        self.name = name
        self.dtype = polars.datatypes.parse_into_dtype(dtype)

    def __eq__(self, other: Field) -> bool:  # type: ignore[override]
        return (self.name == other.name) & (self.dtype == other.dtype)

    def __hash__(self) -> int:
        return hash((self.name, self.dtype))

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}({self.name!r}, {self.dtype})"


class Struct(NestedType):
    """
    Struct composite type.

    Parameters
    ----------
    fields
        The fields that make up the struct. Can be either a sequence of Field
        objects or a mapping of column names to data types.

    Examples
    --------
    Initialize using a dictionary:

    >>> dtype = pl.Struct({"a": pl.Int8, "b": pl.List(pl.String)})
    >>> dtype
    Struct({'a': Int8, 'b': List(String)})

    Initialize using a list of Field objects:

    >>> dtype = pl.Struct([pl.Field("a", pl.Int8), pl.Field("b", pl.List(pl.String))])
    >>> dtype
    Struct({'a': Int8, 'b': List(String)})

    When initializing a Series, Polars can infer a struct data type from the data.

    >>> s = pl.Series([{"a": 1, "b": ["x", "y"]}, {"a": 2, "b": ["z"]}])
    >>> s
    shape: (2,)
    Series: '' [struct[2]]
    [
            {1,["x", "y"]}
            {2,["z"]}
    ]
    >>> s.dtype
    Struct({'a': Int64, 'b': List(String)})
    """

    fields: list[Field]

    def __init__(self, fields: Sequence[Field] | SchemaDict) -> None:
        if isinstance(fields, Mapping):
            self.fields = [Field(name, dtype) for name, dtype in fields.items()]
        else:
            self.fields = list(fields)

    def __eq__(self, other: PolarsDataType) -> bool:  # type: ignore[override]
        # The comparison allows comparing objects to classes, and specific
        # inner types to those without (eg: inner=None). if one of the
        # arguments is not specific about its inner type we infer it
        # as being equal. (See the List type for more info).
        if isclass(other) and issubclass(other, Struct):
            return True
        elif isinstance(other, Struct):
            return self.fields == other.fields
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.__class__, tuple(self.fields)))

    def __iter__(self) -> Iterator[tuple[str, PolarsDataType]]:
        for fld in self.fields:
            yield fld.name, fld.dtype

    def __reversed__(self) -> Iterator[tuple[str, PolarsDataType]]:
        for fld in reversed(self.fields):
            yield fld.name, fld.dtype

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}({dict(self)})"

    def to_schema(self) -> OrderedDict[str, PolarsDataType]:
        """Return Struct dtype as a schema dict."""
        return OrderedDict(self)
