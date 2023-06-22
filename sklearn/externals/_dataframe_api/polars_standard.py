# The files in this folder is vendoered from
# https://github.com/MarcoGorelli/impl-dataframe-api

from __future__ import annotations
import collections

from typing import Sequence, NoReturn, Any, Mapping
import polars as pl
import polars


def dataframe_standard(df: pl.DataFrame) -> PolarsDataFrame:
    return PolarsDataFrame(df)


DTYPE_MAPPING = {  # todo, expand
    "bool": pl.Boolean,
    "int64": pl.Int64,
    "float64": pl.Float64,
}


class PolarsNamespace:
    @classmethod
    def concat(cls, dataframes: Sequence[PolarsDataFrame]) -> PolarsDataFrame:
        dfs = []
        for _df in dataframes:
            dfs.append(_df.dataframe)
        return PolarsDataFrame(pl.concat(dfs))

    @classmethod
    def dataframe_from_dict(cls, data: dict[str, PolarsColumn]) -> PolarsDataFrame:
        return PolarsDataFrame(
            pl.DataFrame({label: column._series for label, column in data.items()})
        )

    @classmethod
    def column_from_sequence(
        cls, sequence: Sequence[object], dtype: str
    ) -> PolarsColumn:
        return PolarsColumn(pl.Series(sequence, dtype=DTYPE_MAPPING[dtype]))


class PolarsColumn:
    def __init__(self, column: pl.Series) -> None:
        self._series = column

    # In the standard
    def __column_namespace__(self, *, api_version: str | None = None) -> Any:
        return PolarsNamespace

    def __len__(self) -> int:
        return len(self._series)

    @property
    def dtype(self) -> object:
        # todo change
        return self._series.dtype

    def get_rows(self, indices: PolarsColumn) -> PolarsColumn:
        return PolarsColumn(self._series.take(indices._series))

    def __getitem__(self, row: int) -> object:
        return self._series[row]

    def __iter__(self) -> NoReturn:
        raise NotImplementedError()

    def is_in(self, values: PolarsColumn) -> PolarsColumn:
        if values.dtype != self.dtype:
            raise ValueError(f"`value` has dtype {values.dtype}, expected {self.dtype}")
        return PolarsColumn(self._series.is_in(values._series))

    def unique(self) -> PolarsColumn:
        return PolarsColumn(self._series.unique())

    def mean(self) -> object:
        return self._series.mean()

    def isnull(self) -> PolarsColumn:
        return PolarsColumn(self._series.is_null())

    def isnan(self) -> PolarsColumn:
        return PolarsColumn(self._series.is_nan())

    def any(self) -> bool:
        return self._series.any()

    def all(self) -> bool:
        return self._series.all()

    def __eq__(  # type: ignore[override]
        self, other: PolarsColumn | object
    ) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series == other._series)
        return PolarsColumn(self._series == other)

    def __ne__(  # type: ignore[override]
        self, other: PolarsColumn | object
    ) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series != other._series)
        return PolarsColumn(self._series != other)

    def __ge__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series >= other._series)
        return PolarsColumn(self._series >= other)

    def __gt__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series > other._series)
        return PolarsColumn(self._series > other)

    def __le__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series <= other._series)
        return PolarsColumn(self._series <= other)

    def __lt__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series < other._series)
        return PolarsColumn(self._series < other)

    def __mul__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series * other._series)
        return PolarsColumn(self._series * other)

    def __floordiv__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series // other._series)
        return PolarsColumn(self._series // other)

    def __truediv__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series / other._series)
        return PolarsColumn(self._series / other)

    def __pow__(self, other: PolarsColumn | float) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series.pow(other._series))
        return PolarsColumn(self._series.pow(other))

    def __mod__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series % other._series)
        return PolarsColumn(self._series % other)

    def __and__(self, other: PolarsColumn | object) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series & other._series)
        return PolarsColumn(self._series & other)  # type: ignore[operator]

    def __invert__(self) -> PolarsColumn:
        return PolarsColumn(~self._series)

    def max(self) -> object:
        return self._series.max()

    def std(self) -> object:
        return self._series.std()

    def __add__(self, other: PolarsColumn) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series + other._series)
        return PolarsColumn(self._series + other)

    def __sub__(self, other: PolarsColumn) -> PolarsColumn:
        if isinstance(other, PolarsColumn):
            return PolarsColumn(self._series - other._series)
        return PolarsColumn(self._series - other)

    def sorted_indices(self) -> PolarsColumn:
        df = self._series.to_frame()
        keys = df.columns
        return PolarsColumn(df.with_row_count().sort(keys, descending=False)["row_nr"])

    def fill_nan(self, value: float) -> PolarsColumn:
        return PolarsColumn(self._series.fill_nan(value))


class PolarsGroupBy:
    def __init__(self, df: pl.DataFrame, keys: Sequence[str]) -> None:
        for key in keys:
            if key not in df.columns:
                raise KeyError(f"key {key} not present in DataFrame's columns")
        self.df = df
        self.keys = keys

    def size(self) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).count().rename({"count": "size"})
        return PolarsDataFrame(result)

    def any(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").any())
        return PolarsDataFrame(result)

    def all(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").all())
        return PolarsDataFrame(result)

    def min(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").min())
        return PolarsDataFrame(result)

    def max(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").max())
        return PolarsDataFrame(result)

    def sum(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").sum())
        return PolarsDataFrame(result)

    def prod(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").product())
        return PolarsDataFrame(result)

    def median(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").median())
        return PolarsDataFrame(result)

    def mean(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").mean())
        return PolarsDataFrame(result)

    def std(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").std())
        return PolarsDataFrame(result)

    def var(self, skipna: bool = True) -> PolarsDataFrame:
        result = self.df.groupby(self.keys).agg(pl.col("*").var())
        return PolarsDataFrame(result)


class PolarsDataFrame:
    def __init__(self, df: pl.DataFrame) -> None:
        # columns already have to be strings, and duplicates aren't
        # allowed, so no validation required
        self.df = df

    def __dataframe_namespace__(self, *, api_version: str | None = None) -> Any:
        return PolarsNamespace

    def __len__(self) -> int:
        return len(self.df)

    @property
    def dataframe(self) -> pl.DataFrame:
        return self.df

    def shape(self) -> tuple[int, int]:
        return self.df.shape

    def groupby(self, keys: Sequence[str]) -> PolarsGroupBy:
        return PolarsGroupBy(self.df, keys)

    def get_column_by_name(self, name: str) -> PolarsColumn:
        return PolarsColumn(self.df[name])

    def get_columns_by_name(self, names: Sequence[str]) -> PolarsDataFrame:
        if isinstance(names, str):
            raise TypeError(f"Expected sequence of str, got {type(names)}")
        return PolarsDataFrame(self.df.select(names))

    def get_rows(self, indices: PolarsColumn) -> PolarsDataFrame:
        return PolarsDataFrame(self.df[indices._series])

    def slice_rows(self, start: int, stop: int, step: int) -> PolarsDataFrame:
        return PolarsDataFrame(
            self.df.with_row_count("idx")
            .filter(pl.col("idx").is_in(range(start, stop, step)))
            .drop("idx")
        )

    def get_rows_by_mask(self, mask: PolarsColumn) -> PolarsDataFrame:
        return PolarsDataFrame(self.df.filter(mask._series))

    def insert(self, loc: int, label: str, value: PolarsColumn) -> PolarsDataFrame:
        df = self.df.clone()
        if len(df) > 0:
            df.insert_at_idx(loc, pl.Series(label, value._series))
            return PolarsDataFrame(df)
        return PolarsDataFrame(pl.DataFrame({label: value._series}))

    def drop_column(self, label: str) -> PolarsDataFrame:
        if not isinstance(label, str):
            raise TypeError(f"Expected str, got: {type(label)}")
        return PolarsDataFrame(self.dataframe.drop(label))

    def rename_columns(self, mapping: Mapping[str, str]) -> PolarsDataFrame:
        if not isinstance(mapping, collections.abc.Mapping):
            raise TypeError(f"Expected Mapping, got: {type(mapping)}")
        return PolarsDataFrame(self.dataframe.rename(dict(mapping)))

    def get_column_names(self) -> Sequence[str]:
        return self.dataframe.columns

    def __eq__(self, other: PolarsDataFrame) -> PolarsDataFrame:  # type: ignore[override]
        return PolarsDataFrame(self.dataframe.__eq__(other.dataframe))

    def __ne__(self, other: PolarsDataFrame) -> PolarsDataFrame:  # type: ignore[override]
        return PolarsDataFrame(self.dataframe.__ne__(other.dataframe))

    def __ge__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__ge__(other.dataframe))

    def __gt__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__gt__(other.dataframe))

    def __le__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__le__(other.dataframe))

    def __lt__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__lt__(other.dataframe))

    def __add__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__add__(other.dataframe))

    def __sub__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__sub__(other.dataframe))

    def __mul__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__mul__(other.dataframe))

    def __truediv__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__truediv__(other.dataframe))

    def __floordiv__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__floordiv__(other.dataframe))

    def __pow__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(
            self.dataframe.select(
                [
                    pl.col(col).pow(other.dataframe[col])
                    for col in self.get_column_names()
                ]
            )
        )

    def __mod__(self, other: PolarsDataFrame) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.__mod__(other.dataframe))

    def __divmod__(
        self,
        other: PolarsDataFrame,
    ) -> tuple[PolarsDataFrame, PolarsDataFrame]:
        quotient = self // other
        remainder = self - quotient * other
        return quotient, remainder

    def isnull(self) -> PolarsDataFrame:
        result = {}
        for column in self.dataframe.columns:
            result[column] = self.dataframe[column].is_null()
        return PolarsDataFrame(pl.DataFrame(result))

    def isnan(self) -> PolarsDataFrame:
        result = {}
        for column in self.dataframe.columns:
            result[column] = self.dataframe[column].is_nan()
        return PolarsDataFrame(pl.DataFrame(result))

    def any(self) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.select(pl.col("*").any()))

    def all(self) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.select(pl.col("*").all()))

    def any_rowwise(self) -> PolarsColumn:
        return PolarsColumn(self.dataframe.select(pl.any(pl.col("*")))["any"])

    def all_rowwise(self) -> PolarsColumn:
        return PolarsColumn(self.dataframe.select(pl.all(pl.col("*")))["all"])

    def sorted_indices(self, keys: Sequence[str]) -> PolarsColumn:
        df = self.dataframe.select(keys)
        return PolarsColumn(df.with_row_count().sort(keys, descending=False)["row_nr"])

    def fill_nan(self, value: float) -> PolarsDataFrame:
        return PolarsDataFrame(self.dataframe.fill_nan(value))
