Polars DataFrame support in permutation_importance
---------------------------------------------------

:class:`~sklearn.inspection.permutation_importance` now supports
Polars DataFrames when used in pipelines with :class:`~sklearn.compose.ColumnTransformer`
and string column names. Previously, Polars DataFrames would be converted to NumPy arrays,
causing string-based column selection to fail.
