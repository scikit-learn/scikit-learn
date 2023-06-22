# The files in this folder is vendoered from
# https://github.com/MarcoGorelli/impl-dataframe-api

from __future__ import annotations
import numpy as np

import pandas as pd
from pandas.api.types import is_extension_array_dtype
import collections
from typing import Any, Sequence, Mapping, NoReturn, cast

import pandas


def dataframe_standard(df: pd.DataFrame) -> PandasDataFrame:
    return PandasDataFrame(df)


class PandasNamespace:
    @classmethod
    def concat(cls, dataframes: Sequence[PandasDataFrame]) -> PandasDataFrame:
        dtypes = dataframes[0].dataframe.dtypes
        dfs = []
        for _df in dataframes:
            try:
                pd.testing.assert_series_equal(_df.dataframe.dtypes, dtypes)
            except Exception as exc:
                raise ValueError("Expected matching columns") from exc
            else:
                dfs.append(_df.dataframe)
        return PandasDataFrame(
            pd.concat(
                dfs,
                axis=0,
                ignore_index=True,
            )
        )

    @classmethod
    def column_from_sequence(
        cls, sequence: Sequence[object], dtype: str
    ) -> PandasColumn:
        return PandasColumn(pd.Series(sequence, dtype=dtype))

    @classmethod
    def dataframe_from_dict(cls, data: dict[str, PandasColumn]) -> PandasDataFrame:
        return PandasDataFrame(
            pd.DataFrame({label: column._series for label, column in data.items()})
        )


class PandasColumn:
    # private, not technically part of the standard
    def __init__(self, column: pd.Series) -> None:  # type: ignore[type-arg]
        if (
            isinstance(column.index, pd.RangeIndex)
            and column.index.start == 0  # type: ignore[comparison-overlap]
            and column.index.step == 1  # type: ignore[comparison-overlap]
            and (column.index.stop == len(column))  # type: ignore[comparison-overlap]
        ):
            self._series = column
        else:
            self._series = column.reset_index(drop=True)

    # In the standard
    def __column_namespace__(self, *, api_version: str | None = None) -> Any:
        return PandasNamespace

    def __len__(self) -> int:
        return len(self._series)

    def __iter__(self) -> NoReturn:
        raise NotImplementedError()

    @property
    def dtype(self) -> object:
        return self._series.dtype

    def get_rows(self, indices: PandasColumn) -> PandasColumn:
        return PandasColumn(self._series.iloc[indices._series.to_numpy()])

    def __getitem__(self, row: int) -> object:
        return self._series.iloc[row]

    def sorted_indices(self) -> PandasColumn:
        return PandasColumn(pd.Series(self._series.argsort()))

    def is_in(self, values: PandasColumn) -> PandasColumn:
        if values.dtype != self.dtype:
            raise ValueError(f"`value` has dtype {values.dtype}, expected {self.dtype}")
        return PandasColumn(self._series.isin(values._series))

    def unique(self) -> PandasColumn:
        return PandasColumn(pd.Series(self._series.unique()))

    def mean(self) -> float:
        return self._series.mean()

    def std(self) -> float:
        return self._series.std()

    def isnull(self) -> PandasColumn:
        if is_extension_array_dtype(self._series.dtype):
            return PandasColumn(self._series.isnull())
        else:
            return PandasColumn(pd.Series(np.array([False] * len(self))))

    def isnan(self) -> PandasColumn:
        return PandasColumn(self._series.isna())

    def any(self) -> bool:
        return self._series.any()

    def all(self) -> bool:
        return self._series.all()

    def __eq__(  # type: ignore[override]
        self, other: PandasColumn | object
    ) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series == other._series)
        return PandasColumn(self._series == other)

    def __ne__(self, other: PandasColumn) -> PandasColumn:  # type: ignore[override]
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series != other._series)
        return PandasColumn(self._series != other)

    def __ge__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series >= other._series)
        return PandasColumn(self._series >= other)

    def __gt__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series > other._series)
        return PandasColumn(self._series > other)

    def __le__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series <= other._series)
        return PandasColumn(self._series <= other)

    def __lt__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series < other._series)
        return PandasColumn(self._series < other)

    def __add__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series + other._series)
        return PandasColumn(self._series + other)

    def __sub__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series - other._series)
        return PandasColumn(self._series - other)

    def __mul__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series * other._series)
        return PandasColumn(self._series * other)

    def __truediv__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series / other._series)
        return PandasColumn(self._series / other)

    def __floordiv__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series // other._series)
        return PandasColumn(self._series // other)

    def __pow__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series**other._series)
        return PandasColumn(self._series**other)

    def __mod__(self, other: PandasColumn) -> PandasColumn:
        if isinstance(other, PandasColumn):
            return PandasColumn(self._series % other._series)
        return PandasColumn(self._series % other)

    def __invert__(self) -> PandasColumn:
        return PandasColumn(~self._series)

    def __and__(self, other: PandasColumn) -> PandasColumn:
        return PandasColumn(self._series & other._series)

    def max(self) -> object:
        return self._series.max()

    def fill_nan(
        self, value: float | pd.NAType  # type: ignore[name-defined]
    ) -> PandasColumn:
        ser = self._series.copy()
        ser[cast("pd.Series[bool]", np.isnan(ser)).fillna(False).to_numpy(bool)] = value
        return PandasColumn(ser)


class PandasGroupBy:
    def __init__(self, df: pd.DataFrame, keys: Sequence[str]) -> None:
        self.df = df
        self.grouped = df.groupby(list(keys), sort=False, as_index=False)
        self.keys = list(keys)

    def _validate_result(self, result: pd.DataFrame) -> None:
        failed_columns = self.df.columns.difference(result.columns)
        if len(failed_columns) > 0:
            # defensive check
            raise AssertionError(
                "Groupby operation could not be performed on columns "
                f"{failed_columns}. Please drop them before calling groupby."
            )

    def size(self) -> PandasDataFrame:
        # pandas-stubs is wrong
        return PandasDataFrame(self.grouped.size())  # type: ignore[arg-type]

    def any(self, skipna: bool = True) -> PandasDataFrame:
        if not (self.df.drop(columns=self.keys).dtypes == "bool").all():
            raise ValueError("Expected boolean types")
        result = self.grouped.any()
        self._validate_result(result)
        return PandasDataFrame(result)

    def all(self, skipna: bool = True) -> PandasDataFrame:
        if not (self.df.drop(columns=self.keys).dtypes == "bool").all():
            raise ValueError("Expected boolean types")
        result = self.grouped.all()
        self._validate_result(result)
        return PandasDataFrame(result)

    def min(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.min()
        self._validate_result(result)
        return PandasDataFrame(result)

    def max(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.max()
        self._validate_result(result)
        return PandasDataFrame(result)

    def sum(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.sum()
        self._validate_result(result)
        return PandasDataFrame(result)

    def prod(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.prod()
        self._validate_result(result)
        return PandasDataFrame(result)

    def median(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.median()
        self._validate_result(result)
        return PandasDataFrame(result)

    def mean(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.mean()
        self._validate_result(result)
        return PandasDataFrame(result)

    def std(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.std()
        self._validate_result(result)
        return PandasDataFrame(result)

    def var(self, skipna: bool = True) -> PandasDataFrame:
        result = self.grouped.var()
        self._validate_result(result)
        return PandasDataFrame(result)


class PandasDataFrame:
    # Not technically part of the standard

    def __init__(self, dataframe: pd.DataFrame) -> None:
        self._validate_columns(dataframe.columns)  # type: ignore[arg-type]
        if (
            isinstance(dataframe.index, pd.RangeIndex)
            and dataframe.index.start == 0  # type: ignore[comparison-overlap]
            and dataframe.index.step == 1  # type: ignore[comparison-overlap]
            and (
                dataframe.index.stop == len(dataframe)  # type: ignore[comparison-overlap]
            )
        ):
            self._dataframe = dataframe
        else:
            self._dataframe = dataframe.reset_index(drop=True)

    def __len__(self) -> int:
        return self.shape()[0]

    def _validate_columns(self, columns: Sequence[str]) -> None:
        counter = collections.Counter(columns)
        for col, count in counter.items():
            if count > 1:
                raise ValueError(
                    f"Expected unique column names, got {col} {count} time(s)"
                )
        for col in columns:
            if not isinstance(col, str):
                raise TypeError(
                    f"Expected column names to be of type str, got {col} "
                    f"of type {type(col)}"
                )

    def _validate_index(self, index: pd.Index) -> None:
        pd.testing.assert_index_equal(self.dataframe.index, index)

    def _validate_comparand(self, other: PandasDataFrame) -> None:
        if isinstance(other, PandasDataFrame) and not (
            self.dataframe.index.equals(other.dataframe.index)
            and self.dataframe.shape == other.dataframe.shape
            and self.dataframe.columns.equals(other.dataframe.columns)
        ):
            raise ValueError(
                "Expected DataFrame with same length, matching columns, "
                "and matching index."
            )

    def _validate_booleanness(self) -> None:
        if not (self.dataframe.dtypes == "bool").all():
            raise NotImplementedError(
                "'any' can only be called on DataFrame where all dtypes are 'bool'"
            )

    # In the standard
    def __dataframe_namespace__(self, *, api_version: str | None = None) -> Any:
        return PandasNamespace

    @property
    def dataframe(self) -> pd.DataFrame:
        return self._dataframe

    def shape(self) -> tuple[int, int]:
        return self.dataframe.shape

    def groupby(self, keys: Sequence[str]) -> PandasGroupBy:
        if not isinstance(keys, collections.abc.Sequence):
            raise TypeError(f"Expected sequence of strings, got: {type(keys)}")
        if isinstance(keys, str):
            raise TypeError("Expected sequence of strings, got: str")
        for key in keys:
            if key not in self.get_column_names():
                raise KeyError(f"key {key} not present in DataFrame's columns")
        return PandasGroupBy(self.dataframe, keys)

    def get_column_by_name(self, name: str) -> PandasColumn:
        if not isinstance(name, str):
            raise ValueError(f"Expected str, got: {type(name)}")
        return PandasColumn(self.dataframe.loc[:, name])

    def get_columns_by_name(self, names: Sequence[str]) -> PandasDataFrame:
        if isinstance(names, str):
            raise TypeError(f"Expected sequence of str, got {type(names)}")
        self._validate_columns(names)
        return PandasDataFrame(self.dataframe.loc[:, list(names)])

    def get_rows(self, indices: PandasColumn) -> PandasDataFrame:
        return PandasDataFrame(self.dataframe.iloc[indices._series, :])

    def slice_rows(self, start: int, stop: int, step: int) -> PandasDataFrame:
        return PandasDataFrame(self.dataframe.iloc[start:stop:step])

    def get_rows_by_mask(self, mask: PandasColumn) -> PandasDataFrame:
        series = mask._series
        self._validate_index(series.index)
        return PandasDataFrame(self.dataframe.loc[series, :])

    def insert(self, loc: int, label: str, value: PandasColumn) -> PandasDataFrame:
        series = value._series
        self._validate_index(series.index)
        before = self.dataframe.iloc[:, :loc]
        after = self.dataframe.iloc[:, loc:]
        to_insert = value._series.rename(label)
        return PandasDataFrame(pd.concat([before, to_insert, after], axis=1))

    def drop_column(self, label: str) -> PandasDataFrame:
        if not isinstance(label, str):
            raise TypeError(f"Expected str, got: {type(label)}")
        return PandasDataFrame(self.dataframe.drop(label, axis=1))

    def rename_columns(self, mapping: Mapping[str, str]) -> PandasDataFrame:
        if not isinstance(mapping, collections.abc.Mapping):
            raise TypeError(f"Expected Mapping, got: {type(mapping)}")
        return PandasDataFrame(self.dataframe.rename(columns=mapping))

    def get_column_names(self) -> Sequence[str]:
        return self.dataframe.columns.tolist()

    def __iter__(self) -> NoReturn:
        raise NotImplementedError()

    def __eq__(self, other: PandasDataFrame) -> PandasDataFrame:  # type: ignore[override]
        self._validate_comparand(other)
        return PandasDataFrame(self.dataframe.__eq__(other.dataframe))

    def __ne__(self, other: PandasDataFrame) -> PandasDataFrame:  # type: ignore[override]
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__ne__(other.dataframe)))

    def __ge__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__ge__(other.dataframe)))

    def __gt__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__gt__(other.dataframe)))

    def __le__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__le__(other.dataframe)))

    def __lt__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__lt__(other.dataframe)))

    def __add__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__add__(other.dataframe)))

    def __sub__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__sub__(other.dataframe)))

    def __mul__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__mul__(other.dataframe)))

    def __truediv__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__truediv__(other.dataframe)))

    def __floordiv__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__floordiv__(other.dataframe)))

    def __pow__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__pow__(other.dataframe)))

    def __mod__(self, other: PandasDataFrame) -> PandasDataFrame:
        self._validate_comparand(other)
        return PandasDataFrame((self.dataframe.__mod__(other.dataframe)))

    def __divmod__(
        self, other: PandasDataFrame
    ) -> tuple[PandasDataFrame, PandasDataFrame]:
        self._validate_comparand(other)
        quotient, remainder = self.dataframe.__divmod__(other.dataframe)
        return PandasDataFrame(quotient), PandasDataFrame(remainder)

    def any(self) -> PandasDataFrame:
        self._validate_booleanness()
        return PandasDataFrame(self.dataframe.any().to_frame().T)

    def all(self) -> PandasDataFrame:
        self._validate_booleanness()
        return PandasDataFrame(self.dataframe.all().to_frame().T)

    def any_rowwise(self) -> PandasColumn:
        self._validate_booleanness()
        return PandasColumn(self.dataframe.any(axis=1))

    def all_rowwise(self) -> PandasColumn:
        self._validate_booleanness()
        return PandasColumn(self.dataframe.all(axis=1))

    def isnull(self) -> PandasDataFrame:
        result = []
        for column in self.dataframe.columns:
            if is_extension_array_dtype(self.dataframe[column].dtype):
                result.append(self.dataframe[column].isnull())
            else:
                result.append(
                    pd.Series(np.array([False] * self.shape()[0]), name=column)
                )
        return PandasDataFrame(pd.concat(result, axis=1))

    def isnan(self) -> PandasDataFrame:
        result = []
        for column in self.dataframe.columns:
            if is_extension_array_dtype(self.dataframe[column].dtype):
                result.append(
                    np.isnan(self.dataframe[column]).replace(pd.NA, False).astype(bool)
                )
            else:
                result.append(self.dataframe[column].isna())
        return PandasDataFrame(pd.concat(result, axis=1))

    def sorted_indices(self, keys: Sequence[str]) -> PandasColumn:
        df = self.dataframe.loc[:, list(keys)]
        return PandasColumn(df.sort_values(keys).index.to_series())

    def fill_nan(
        self, value: float | pd.NAType  # type: ignore[name-defined]
    ) -> PandasDataFrame:
        df = self.dataframe.copy()
        df[cast(pd.DataFrame, np.isnan(df)).fillna(False).to_numpy(bool)] = value
        return PandasDataFrame(df)
