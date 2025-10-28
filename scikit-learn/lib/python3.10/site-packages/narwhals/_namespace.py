"""Narwhals-level equivalent of `CompliantNamespace`."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar, Generic, overload

from narwhals._compliant.typing import CompliantNamespaceAny, CompliantNamespaceT_co
from narwhals._native import (
    NativeAny,
    NativeArrow,
    NativeCuDF,
    NativeDask,
    NativeDuckDB,
    NativeIbis,
    NativeModin,
    NativePandas,
    NativePandasLike,
    NativePolars,
    NativeSparkLike,
    NativeUnknown,
    _CuDFDataFrame,
    _CuDFSeries,
    _ModinDataFrame,
    _ModinSeries,
    is_native_arrow,
    is_native_cudf,
    is_native_dask,
    is_native_duckdb,
    is_native_ibis,
    is_native_modin,
    is_native_pandas,
    is_native_polars,
    is_native_pyspark_connect,
    is_native_spark_like,
    is_native_sqlframe,
)
from narwhals._utils import Implementation, Version

if TYPE_CHECKING:
    import pandas as pd
    from typing_extensions import TypeAlias

    from narwhals._arrow.namespace import ArrowNamespace
    from narwhals._dask.namespace import DaskNamespace
    from narwhals._duckdb.namespace import DuckDBNamespace
    from narwhals._ibis.namespace import IbisNamespace
    from narwhals._pandas_like.namespace import PandasLikeNamespace
    from narwhals._polars.namespace import PolarsNamespace
    from narwhals._spark_like.namespace import SparkLikeNamespace
    from narwhals._typing import (
        Arrow,
        Backend,
        Dask,
        DuckDB,
        EagerAllowed,
        Ibis,
        IntoBackend,
        PandasLike,
        Polars,
        SparkLike,
    )

    EagerAllowedNamespace: TypeAlias = "Namespace[PandasLikeNamespace] | Namespace[ArrowNamespace] | Namespace[PolarsNamespace]"

__all__ = ["Namespace"]


class Namespace(Generic[CompliantNamespaceT_co]):
    _compliant_namespace: CompliantNamespaceT_co
    _version: ClassVar[Version] = Version.MAIN

    def __init__(self, namespace: CompliantNamespaceT_co, /) -> None:
        self._compliant_namespace = namespace

    def __init_subclass__(cls, *args: Any, version: Version, **kwds: Any) -> None:
        super().__init_subclass__(*args, **kwds)

        if isinstance(version, Version):
            cls._version = version
        else:
            msg = f"Expected {Version} but got {type(version).__name__!r}"
            raise TypeError(msg)

    def __repr__(self) -> str:
        return f"Namespace[{type(self.compliant).__name__}]"

    @property
    def compliant(self) -> CompliantNamespaceT_co:
        return self._compliant_namespace

    @property
    def implementation(self) -> Implementation:
        return self.compliant._implementation

    @property
    def version(self) -> Version:
        return self._version

    @overload
    @classmethod
    def from_backend(cls, backend: PandasLike, /) -> Namespace[PandasLikeNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: Polars, /) -> Namespace[PolarsNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: Arrow, /) -> Namespace[ArrowNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: SparkLike, /) -> Namespace[SparkLikeNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: DuckDB, /) -> Namespace[DuckDBNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: Dask, /) -> Namespace[DaskNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: Ibis, /) -> Namespace[IbisNamespace]: ...

    @overload
    @classmethod
    def from_backend(cls, backend: EagerAllowed, /) -> EagerAllowedNamespace: ...

    @overload
    @classmethod
    def from_backend(
        cls, backend: IntoBackend[Backend], /
    ) -> Namespace[CompliantNamespaceAny]: ...

    @classmethod
    def from_backend(
        cls: type[Namespace[Any]], backend: IntoBackend[Backend], /
    ) -> Namespace[Any]:
        """Instantiate from native namespace module, string, or Implementation.

        Arguments:
            backend: native namespace module, string, or Implementation.

        Examples:
            >>> from narwhals._namespace import Namespace
            >>> Namespace.from_backend("polars")
            Namespace[PolarsNamespace]
        """
        impl = Implementation.from_backend(backend)
        backend_version = impl._backend_version()  # noqa: F841
        version = cls._version
        ns: CompliantNamespaceAny
        if impl.is_pandas_like():
            from narwhals._pandas_like.namespace import PandasLikeNamespace

            ns = PandasLikeNamespace(implementation=impl, version=version)

        elif impl.is_polars():
            from narwhals._polars.namespace import PolarsNamespace

            ns = PolarsNamespace(version=version)
        elif impl.is_pyarrow():
            from narwhals._arrow.namespace import ArrowNamespace

            ns = ArrowNamespace(version=version)
        elif impl.is_spark_like():
            from narwhals._spark_like.namespace import SparkLikeNamespace

            ns = SparkLikeNamespace(implementation=impl, version=version)
        elif impl.is_duckdb():
            from narwhals._duckdb.namespace import DuckDBNamespace

            ns = DuckDBNamespace(version=version)
        elif impl.is_dask():
            from narwhals._dask.namespace import DaskNamespace

            ns = DaskNamespace(version=version)
        elif impl.is_ibis():
            from narwhals._ibis.namespace import IbisNamespace

            ns = IbisNamespace(version=version)
        else:
            msg = "Not supported Implementation"  # pragma: no cover
            raise AssertionError(msg)
        return cls(ns)

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativePolars, /
    ) -> Namespace[PolarsNamespace]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativePandas, /
    ) -> Namespace[PandasLikeNamespace[pd.DataFrame, pd.Series[Any]]]: ...

    @overload
    @classmethod
    def from_native_object(cls, native: NativeArrow, /) -> Namespace[ArrowNamespace]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativeSparkLike, /
    ) -> Namespace[SparkLikeNamespace]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativeDuckDB, /
    ) -> Namespace[DuckDBNamespace]: ...

    @overload
    @classmethod
    def from_native_object(cls, native: NativeDask, /) -> Namespace[DaskNamespace]: ...

    @overload
    @classmethod
    def from_native_object(cls, native: NativeIbis, /) -> Namespace[IbisNamespace]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativeModin, /
    ) -> Namespace[PandasLikeNamespace[_ModinDataFrame, _ModinSeries]]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativeCuDF, /
    ) -> Namespace[PandasLikeNamespace[_CuDFDataFrame, _CuDFSeries]]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativePandasLike, /
    ) -> Namespace[PandasLikeNamespace[Any, Any]]: ...

    @overload
    @classmethod
    def from_native_object(
        cls, native: NativeUnknown, /
    ) -> Namespace[CompliantNamespaceAny]: ...

    @classmethod
    def from_native_object(
        cls: type[Namespace[Any]], native: NativeAny, /
    ) -> Namespace[Any]:
        impl: Backend
        if is_native_polars(native):
            impl = Implementation.POLARS
        elif is_native_pandas(native):
            impl = Implementation.PANDAS
        elif is_native_arrow(native):
            impl = Implementation.PYARROW
        elif is_native_spark_like(native):
            impl = (
                Implementation.SQLFRAME
                if is_native_sqlframe(native)
                else Implementation.PYSPARK_CONNECT
                if is_native_pyspark_connect(native)
                else Implementation.PYSPARK
            )
        elif is_native_dask(native):  # pragma: no cover
            impl = Implementation.DASK
        elif is_native_duckdb(native):
            impl = Implementation.DUCKDB
        elif is_native_cudf(native):  # pragma: no cover
            impl = Implementation.CUDF
        elif is_native_modin(native):  # pragma: no cover
            impl = Implementation.MODIN
        elif is_native_ibis(native):
            impl = Implementation.IBIS
        else:
            msg = f"Unsupported type: {type(native).__qualname__!r}"
            raise TypeError(msg)
        return cls.from_backend(impl)
